import numpy as np
import spiceypy as spice
from scipy.optimize import differential_evolution
import planetary_data as pd
from Propagator import Propagator
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from datetime import datetime
import time


class FlybyOptimizer:

    def __init__(self,
                 target_body,
                 central_body,
                 start_date,
                 end_date,
                 resonant_orbit,
                 perturbations,
                 alt_min_km=100,
                 alt_max_km=5000,
                 e_delta=0.05,
                 steps_coarse=3000,
                 steps_fine=10000,
                 steps_final=100000,
                 popsize_coarse=15,
                 popsize_fine=10,
                 maxiter_coarse=100,
                 maxiter_fine=80,
                 fine_delta_deg=20.0,
                 fine_e_delta=0.01,
                 seed=42):
        """
        Parameters
        ----------
        target_body      : str — SPICE body name e.g. "EUROPA"
        central_body     : planetary_data dict e.g. planetary_data.jupiter
        start_date       : str — propagation start
        end_date         : str — propagation end
        resonant_orbit   : list [a, e, i, M0, argp, raan] in km and radians
                           from get_resonant_orbit()
        perturbations    : list of perturbation dicts
        alt_min_km       : minimum flyby altitude (km)
        alt_max_km       : maximum flyby altitude (km)
        e_delta          : coarse search range around resonant eccentricity (±)
        steps_coarse     : propagation steps during coarse optimization
        steps_fine       : propagation steps during fine optimization
        steps_final      : propagation steps for final single propagation
        popsize_coarse   : differential evolution population size — coarse phase
        popsize_fine     : differential evolution population size — fine phase
        maxiter_coarse   : max iterations — coarse phase
        maxiter_fine     : max iterations — fine phase
        fine_delta_deg   : narrow search window around coarse solution for
                           M0 and argp in fine phase (degrees, ±)
        fine_e_delta     : narrow search window around coarse solution for
                           e in fine phase (±)
        seed             : random seed for reproducibility
        """

        self.target_body     = target_body
        self.central_body    = central_body
        self.start_date      = start_date
        self.end_date        = end_date
        self.resonant_orbit  = resonant_orbit
        self.perturbations   = perturbations
        self.alt_min_km      = alt_min_km
        self.alt_max_km      = alt_max_km
        self.e_delta         = e_delta
        self.steps_coarse    = steps_coarse
        self.steps_fine      = steps_fine
        self.steps_final     = steps_final
        self.popsize_coarse  = popsize_coarse
        self.popsize_fine    = popsize_fine
        self.maxiter_coarse  = maxiter_coarse
        self.maxiter_fine    = maxiter_fine
        self.fine_delta_rad  = np.radians(fine_delta_deg)
        self.fine_e_delta    = fine_e_delta
        self.seed            = seed

        # body radius for altitude computation
        self.body_radius = pd.get_body_radius(target_body)

        # unpack resonant orbit for reference
        self.res_a    = resonant_orbit[0]
        self.res_e    = resonant_orbit[1]
        self.res_i    = resonant_orbit[2]
        self.res_M0   = resonant_orbit[3]
        self.res_argp = resonant_orbit[4]
        self.res_raan = resonant_orbit[5]

        # results — populated by optimize()
        self.coarse_result  = None   # best (M0, argp, e) from phase 1
        self.best_M0        = None
        self.best_argp      = None
        self.best_e         = None
        self.best_count     = 0
        self.best_flybys    = []
        self.final_prop     = None

        # full history across both phases: (eval, M0, argp, e, count, phase)
        self.history        = []

        # internal state used during a phase run
        self._eval_count    = 0
        self._t_start       = None
        self._best_so_far   = 0
        self._current_phase = 1
        self._current_steps = steps_coarse
        self._current_alt_max = None   # set per phase

    # --------------------------------------------------
    def _propagate(self, M0, argp, e, steps):
        """
        Run propagator with given M0, argp, e.
        a, i, raan are fixed from the resonant orbit.
        """
        orbit = self.resonant_orbit.copy()
        orbit[1] = e
        orbit[3] = M0
        orbit[4] = argp

        prop = Propagator(
            orbit, steps,
            start_date=self.start_date,
            end_date=self.end_date,
            coes=True, deg=False,
            cb=self.central_body,
            perts=self.perturbations
        )
        prop.propagate()
        return prop

    # --------------------------------------------------
    def find_all_flybys(self, rs, times, alt_max_km=None):
        """
        Find all local minima in distance to target body within altitude bounds.
        alt_max_km overrides self.alt_max_km if provided — used in coarse phase
        to count approaches within a generous window.
        """
        if alt_max_km is None:
            alt_max_km = self.alt_max_km

        r_body = np.array([
            spice.spkpos(self.target_body, float(t),
                         'J2000', 'NONE',
                         self.central_body['name'].upper())[0]
            for t in times
        ])

        dist = np.linalg.norm(rs - r_body, axis=1) - self.body_radius

        flybys = []
        for i in range(1, len(dist) - 1):
            if (dist[i] < dist[i - 1] and
                dist[i] < dist[i + 1] and
                self.alt_min_km <= dist[i] <= alt_max_km):

                flybys.append({
                    "idx":      i,
                    "time_et":  times[i],
                    "time_utc": spice.et2utc(float(times[i]), "C", 0),
                    "altitude": dist[i],
                })

        return flybys

    # --------------------------------------------------
    def _objective(self, u):
        """
        Objective function used by both phases.
        Phase 1 uses a generous alt_max to guide toward close approaches.
        Phase 2 uses the real alt_max to count actual flybys.
        """
        M0, argp, e = u

        self._eval_count += 1

        try:
            prop   = self._propagate(M0, argp, e, self._current_steps)
            flybys = self.find_all_flybys(
                prop.rs, prop.times,
                alt_max_km=self._current_alt_max
            )
            count  = len(flybys)

            # secondary objective: prefer closer approaches
            if count > 0:
                mean_alt = np.mean([fb["altitude"] for fb in flybys])
                score = count + (1.0 - mean_alt / self._current_alt_max) * 0.1
            else:
                score = 0.0

        except Exception:
            score = 0.0
            count = 0

        if count > self._best_so_far:
            self._best_so_far = count

        self.history.append((
            self._eval_count, M0, argp, e, count, self._current_phase
        ))

        if self._eval_count % 10 == 0:
            elapsed = time.time() - self._t_start
            rate    = self._eval_count / elapsed
            print(f"\r  [phase {self._current_phase}]  "
                  f"eval {self._eval_count:4d}   "
                  f"current: {count:2d}   "
                  f"best: {self._best_so_far:2d}   "
                  f"M0: {np.degrees(M0):.1f}°   "
                  f"argp: {np.degrees(argp):.1f}°   "
                  f"e: {e:.4f}   "
                  f"{elapsed:.0f}s  {rate:.1f} eval/s",
                  end="", flush=True)

        return -score

    # --------------------------------------------------
    def _callback(self, xk, convergence):
        M0, argp, e = xk
        gen = self._eval_count // (
            self.popsize_coarse if self._current_phase == 1
            else self.popsize_fine
        )
        print(f"\n  [phase {self._current_phase}] gen {gen:3d}   "
              f"best: {self._best_so_far:2d}   "
              f"convergence: {convergence:.4f}   "
              f"M0: {np.degrees(M0):.1f}°   "
              f"argp: {np.degrees(argp):.1f}°   "
              f"e: {e:.4f}")

    # --------------------------------------------------
    def _run_phase(self, bounds, steps, alt_max_km, popsize, maxiter, phase):
        """Run one differential evolution phase with the given settings."""
        self._current_phase   = phase
        self._current_steps   = steps
        self._current_alt_max = alt_max_km
        self._best_so_far     = 0
        self._eval_count      = 0
        self._t_start         = time.time()

        result = differential_evolution(
            self._objective,
            bounds=bounds,
            maxiter=maxiter,
            popsize=popsize,
            tol=0.01,
            seed=self.seed,
            callback=self._callback,
            polish=True,
            workers=1,
        )

        elapsed = time.time() - self._t_start
        print(f"\n  Phase {phase} complete — "
              f"{elapsed:.0f}s   "
              f"best: {self._best_so_far} approaches   "
              f"M0: {np.degrees(result.x[0]):.1f}°   "
              f"argp: {np.degrees(result.x[1]):.1f}°   "
              f"e: {result.x[2]:.4f}")

        return result

    # --------------------------------------------------
    def optimize(self):
        """
        Two-phase optimization:

        Phase 1 — coarse:
            Full (M0, argp, e) search space, generous alt_max_km to count
            all close approaches not just true flybys, coarse timestep.
            Finds the right region of the search space quickly.

        Phase 2 — fine:
            Narrow bounds centered on phase 1 solution, real alt_max_km,
            finer timestep. Resolves actual flybys missed by coarse timestep.

        Final:
            Single high-resolution propagation at best solution.
        """
        res_e = self.res_e

        print(f"\nFlyby Optimizer — {self.central_body['name'].capitalize()} "
              f"system → {self.target_body.capitalize()}")
        print(f"Start: {self.start_date}   End: {self.end_date}")
        print(f"Resonant orbit:")
        print(f"  a    = {self.res_a:.0f} km  (fixed)")
        print(f"  e    = {self.res_e:.4f}  (coarse search ±{self.e_delta})")
        print(f"  i    = {self.res_i:.2f}°  (fixed)")
        print(f"  argp = {self.res_argp:.2f}°  (full search 0–360°)")
        print(f"  raan = {self.res_raan:.2f}°  (fixed)")
        print(f"  M0   = {self.res_M0:.2f}°  (full search 0–360°)")
        print(f"Altitude bounds: {self.alt_min_km} – {self.alt_max_km} km")
        print(f"{'─'*70}")

        # --------------------------------------------------
        # Phase 1 — coarse global search
        # use generous alt_max (10× real) so the objective can see
        # the close approaches at 40-50k km and be guided toward them
        # --------------------------------------------------
        print(f"\nPhase 1 — coarse global search")
        print(f"  steps={self.steps_coarse}   "
              f"alt_max={self.alt_max_km * 10:.0f} km (generous)   "
              f"popsize={self.popsize_coarse}   maxiter={self.maxiter_coarse}")

        coarse_bounds = [
            (0,                              2 * np.pi),
            (0,                              2 * np.pi),
            (max(0.0, res_e - self.e_delta), res_e + self.e_delta),
        ]

        coarse_result = self._run_phase(
            bounds   = coarse_bounds,
            steps    = self.steps_coarse,
            alt_max_km = self.alt_max_km * 10,   # generous window
            popsize  = self.popsize_coarse,
            maxiter  = self.maxiter_coarse,
            phase    = 1,
        )

        M0_c, argp_c, e_c = coarse_result.x
        self.coarse_result = coarse_result.x.copy()

        # --------------------------------------------------
        # Phase 2 — fine search around coarse solution
        # narrow bounds, finer timestep, real alt_max
        # --------------------------------------------------
        print(f"\nPhase 2 — fine search around coarse solution")
        print(f"  steps={self.steps_fine}   "
              f"alt_max={self.alt_max_km:.0f} km (real)   "
              f"popsize={self.popsize_fine}   maxiter={self.maxiter_fine}")
        print(f"  search window: M0 ±{np.degrees(self.fine_delta_rad):.1f}°   "
              f"argp ±{np.degrees(self.fine_delta_rad):.1f}°   "
              f"e ±{self.fine_e_delta}")

        fine_bounds = [
            (M0_c   - self.fine_delta_rad, M0_c   + self.fine_delta_rad),
            (argp_c - self.fine_delta_rad, argp_c + self.fine_delta_rad),
            (max(0.0, e_c - self.fine_e_delta), e_c + self.fine_e_delta),
        ]

        fine_result = self._run_phase(
            bounds     = fine_bounds,
            steps      = self.steps_fine,
            alt_max_km = self.alt_max_km*10,
            popsize    = self.popsize_fine,
            maxiter    = self.maxiter_fine,
            phase      = 2,
        )

        self.best_M0   = fine_result.x[0]
        self.best_argp = fine_result.x[1]
        self.best_e    = fine_result.x[2]

        # --------------------------------------------------
        # Final high-resolution propagation
        # --------------------------------------------------
        print(f"\nFinal propagation ({self.steps_final} steps)...")

        self.final_prop  = self._propagate(
            self.best_M0, self.best_argp, self.best_e, self.steps_final
        )
        self.best_flybys = self.find_all_flybys(
            self.final_prop.rs, self.final_prop.times
        )
        self.best_count = len(self.best_flybys)

        print(f"\n{'═'*70}")
        print(f"RESULT: {self.best_count} flybys")
        print(f"  a    = {self.res_a:.2f} km  (fixed)")
        print(f"  e    = {self.best_e:.4f}  "
              f"(resonant: {self.res_e:.4f}  Δe={self.best_e - self.res_e:+.4f})")
        print(f"  argp = {np.degrees(self.best_argp):.2f}°  "
              f"(resonant: {np.degrees(self.res_argp):.2f}°)")
        print(f"  M0   = {np.degrees(self.best_M0):.2f}°  "
              f"(resonant: {np.degrees(self.res_M0):.2f}°)")
        if self.best_flybys:
            print(f"\nFlyby summary:")
            for k, fb in enumerate(self.best_flybys):
                print(f"  [{k+1:2d}]  {fb['time_utc']}   "
                      f"alt = {fb['altitude']:.0f} km")
        print(f"{'═'*70}\n")

        return fine_result

    # --------------------------------------------------
    def plot_results(self, save_html=None):
        """
        Plot distance to target body over time with flyby markers,
        and optimization history split by phase.
        """
        if self.final_prop is None:
            print("Run optimize() first.")
            return

        prop    = self.final_prop
        times   = prop.times
        ts_days = (np.array(times) - times[0]) / 86400

        r_body = np.array([
            spice.spkpos(self.target_body, float(t),
                         'J2000', 'NONE',
                         self.central_body['name'].upper())[0]
            for t in times
        ])
        dist = np.linalg.norm(prop.rs - r_body, axis=1) - self.body_radius

        fig = make_subplots(
            rows=3, cols=1,
            subplot_titles=(
                f"Distance to {self.target_body.capitalize()} — "
                f"{self.best_count} flybys   "
                f"e={self.best_e:.4f}   "
                f"argp={np.degrees(self.best_argp):.1f}°   "
                f"M0={np.degrees(self.best_M0):.1f}°",
                "Optimization history — flyby count per evaluation",
                "Search space — argp vs M0 (color = flyby count)",
            ),
            vertical_spacing=0.10,
            row_heights=[0.50, 0.25, 0.25],
        )

        # --- distance trace ---
        fig.add_trace(go.Scatter(
            x=ts_days, y=dist,
            mode='lines',
            line=dict(color='royalblue', width=1),
            name='altitude',
            hovertemplate="day: %{x:.1f}<br>alt: %{y:.0f} km<extra></extra>",
        ), row=1, col=1)

        fig.add_hline(y=self.alt_max_km, row=1, col=1,
                      line=dict(color='orange', width=1, dash='dash'),
                      annotation_text=f"max alt {self.alt_max_km} km",
                      annotation_position="right",
                      annotation_font=dict(color='orange', size=9))

        fig.add_hline(y=self.alt_min_km, row=1, col=1,
                      line=dict(color='tomato', width=1, dash='dash'),
                      annotation_text=f"min alt {self.alt_min_km} km",
                      annotation_position="right",
                      annotation_font=dict(color='tomato', size=9))

        fig.add_hline(y=0, row=1, col=1,
                      line=dict(color='grey', width=0.8, dash='dot'),
                      annotation_text=f"{self.target_body.capitalize()} surface",
                      annotation_position="right",
                      annotation_font=dict(color='grey', size=9))

        if self.best_flybys:
            fb_days = [(fb["time_et"] - times[0]) / 86400
                       for fb in self.best_flybys]
            fb_alts = [fb["altitude"] for fb in self.best_flybys]
            fb_text = [f"Flyby {k+1}<br>{fb['time_utc']}<br>"
                       f"alt: {fb['altitude']:.0f} km"
                       for k, fb in enumerate(self.best_flybys)]

            fig.add_trace(go.Scatter(
                x=fb_days, y=fb_alts,
                mode='markers+text',
                marker=dict(color='tomato', size=10, symbol='star'),
                text=[f"{k+1}" for k in range(len(self.best_flybys))],
                textposition='top center',
                textfont=dict(color='tomato', size=9),
                name=f"flybys ({self.best_count})",
                hovertext=fb_text,
                hoverinfo='text',
            ), row=1, col=1)

        # --- optimization history split by phase ---
        if self.history:
            for phase, color, name in [(1, 'steelblue', 'phase 1 (coarse)'),
                                       (2, 'seagreen',  'phase 2 (fine)')]:
                ph = [h for h in self.history if h[5] == phase]
                if not ph:
                    continue
                h_evals  = [h[0] for h in ph]
                h_counts = [h[4] for h in ph]

                fig.add_trace(go.Scatter(
                    x=h_evals, y=h_counts,
                    mode='markers',
                    marker=dict(color=color, size=3, opacity=0.5),
                    name=name,
                    hovertemplate="eval: %{x}<br>count: %{y}<extra></extra>",
                ), row=2, col=1)

            # running best across both phases
            all_counts = [h[4] for h in self.history]
            all_evals  = [h[0] for h in self.history]
            running_best = []
            best = 0
            for c in all_counts:
                best = max(best, c)
                running_best.append(best)

            fig.add_trace(go.Scatter(
                x=all_evals, y=running_best,
                mode='lines',
                line=dict(color='tomato', width=2),
                name='running best',
                hovertemplate="eval: %{x}<br>best: %{y}<extra></extra>",
            ), row=2, col=1)

            # --- argp vs M0 scatter colored by count ---
            for phase, symbol in [(1, 'circle'), (2, 'diamond')]:
                ph = [h for h in self.history if h[5] == phase]
                if not ph:
                    continue
                fig.add_trace(go.Scatter(
                    x=[np.degrees(h[2]) for h in ph],
                    y=[np.degrees(h[1]) for h in ph],
                    mode='markers',
                    marker=dict(
                        color=[h[4] for h in ph],
                        colorscale='Viridis',
                        size=4 if phase == 1 else 6,
                        opacity=0.5,
                        symbol=symbol,
                        colorbar=dict(
                            title="flybys",
                            x=1.02,
                            len=0.28,
                            y=0.13,
                        ) if phase == 2 else None,
                        showscale=(phase == 2),
                    ),
                    name=f"phase {phase} points",
                    hovertemplate=(
                        f"phase {phase}<br>"
                        "argp: %{x:.1f}°<br>"
                        "M0: %{y:.1f}°<br>"
                        "flybys: %{marker.color}<extra></extra>"
                    ),
                ), row=3, col=1)

            # best solution marker
            fig.add_trace(go.Scatter(
                x=[np.degrees(self.best_argp)],
                y=[np.degrees(self.best_M0)],
                mode='markers',
                marker=dict(color='tomato', size=14, symbol='star'),
                name='best solution',
                hovertemplate=(
                    f"BEST<br>"
                    f"argp: {np.degrees(self.best_argp):.1f}°<br>"
                    f"M0: {np.degrees(self.best_M0):.1f}°<br>"
                    f"e: {self.best_e:.4f}<br>"
                    f"flybys: {self.best_count}<extra></extra>"
                ),
            ), row=3, col=1)

        fig.update_layout(
            title=dict(
                text=(
                    f"{self.central_body['name'].capitalize()} system — "
                    f"Flyby optimization: {self.target_body.capitalize()}<br>"
                    f"<sup>{self.start_date} → {self.end_date}   "
                    f"alt: {self.alt_min_km}–{self.alt_max_km} km   "
                    f"result: {self.best_count} flybys   "
                    f"a = {self.res_a:.0f} km (fixed)</sup>"
                ),
                font=dict(size=13, color='black'),
                x=0.5, xanchor='center',
            ),
            paper_bgcolor='white',
            plot_bgcolor='white',
            font=dict(color='black'),
            width=1100,
            height=900,
            legend=dict(
                x=1.02, y=1,
                bordercolor='lightgrey',
                borderwidth=1,
                font=dict(color='black', size=9),
            ),
        )

        fig.update_xaxes(tickfont=dict(color='black'), gridcolor='lightgrey', showgrid=True)
        fig.update_yaxes(tickfont=dict(color='black'), gridcolor='lightgrey', showgrid=True)

        fig.update_xaxes(title=dict(text="Time (days)",  font=dict(color='black')), row=1, col=1)
        fig.update_xaxes(title=dict(text="Evaluation",   font=dict(color='black')), row=2, col=1)
        fig.update_xaxes(title=dict(text="argp (deg)",   font=dict(color='black')), row=3, col=1)
        fig.update_yaxes(title=dict(text="Altitude (km)",font=dict(color='black')), row=1, col=1)
        fig.update_yaxes(title=dict(text="Flyby count",  font=dict(color='black')), row=2, col=1)
        fig.update_yaxes(title=dict(text="M0 (deg)",     font=dict(color='black')), row=3, col=1)

        if save_html:
            fig.write_html(save_html)
            print(f"Saved to {save_html}")

        fig.show()
        return fig