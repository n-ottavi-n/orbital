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
                 a_bounds_pct=2.0,
                 steps_opt=3000,
                 steps_final=100000,
                 popsize=15,
                 maxiter=100,
                 seed=42):
        """
        Parameters
        ----------
        target_body    : str — SPICE body name e.g. "EUROPA"
        central_body   : planetary_data dict e.g. planetary_data.jupiter
        start_date     : str — propagation start
        end_date       : str — propagation end
        resonant_orbit : list [a, e, i, M0, argp, raan] in km and radians
                         from get_resonant_orbit()
        perturbations  : list of perturbation dicts
        alt_min_km     : minimum flyby altitude (km)
        alt_max_km     : maximum flyby altitude (km)
        a_bounds_pct   : search range around resonant a in percent
        steps_opt      : propagation steps during optimization (coarse)
        steps_final    : propagation steps for final solution (fine)
        popsize        : differential evolution population size
        maxiter        : differential evolution max iterations
        seed           : random seed for reproducibility
        """

        self.target_body    = target_body
        self.central_body   = central_body
        self.start_date     = start_date
        self.end_date       = end_date
        self.resonant_orbit = resonant_orbit
        self.perturbations  = perturbations
        self.alt_min_km     = alt_min_km
        self.alt_max_km     = alt_max_km
        self.a_bounds_pct   = a_bounds_pct
        self.steps_opt      = steps_opt
        self.steps_final    = steps_final
        self.popsize        = popsize
        self.maxiter        = maxiter
        self.seed           = seed

        # body radius for altitude computation
        self.body_radius = pd.get_body_radius(target_body)

        # results populated by optimize()
        self.best_a         = None
        self.best_M0        = None
        self.best_count     = 0
        self.best_flybys    = []
        self.final_prop     = None
        self.history        = []   # (iteration, a, M0, count) per evaluation

        # progress tracking
        self._eval_count    = 0
        self._t_start       = None
        self._best_so_far   = 0

    # --------------------------------------------------
    def _propagate(self, a, M0, steps):
        """Run propagator with given a and M0, rest from resonant orbit."""
        orbit = self.resonant_orbit.copy()
        orbit[0] = a
        orbit[3] = M0

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
    def find_all_flybys(self, rs, times):
        """
        Find all local minima in distance to target body
        within altitude bounds.
        """
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
                self.alt_min_km <= dist[i] <= self.alt_max_km):

                flybys.append({
                    "idx":      i,
                    "time_et":  times[i],
                    "time_utc": spice.et2utc(float(times[i]), "C", 0),
                    "altitude": dist[i],
                })

        return flybys

    # --------------------------------------------------
    def _objective(self, u):
        """Objective function for differential evolution."""
        a, M0 = u

        self._eval_count += 1

        try:
            prop   = self._propagate(a, M0, self.steps_opt)
            flybys = self.find_all_flybys(prop.rs, prop.times)
            count  = len(flybys)
        except Exception as e:
            count  = 0

        # update best
        if count > self._best_so_far:
            self._best_so_far = count

        # record history
        self.history.append((self._eval_count, a, M0, count))

        # progress report every 10 evaluations
        if self._eval_count % 10 == 0:
            elapsed   = time.time() - self._t_start
            rate      = self._eval_count / elapsed
            print(f"\r  eval {self._eval_count:4d}   "
                  f"current: {count:2d} flybys   "
                  f"best: {self._best_so_far:2d} flybys   "
                  f"a: {a:.0f} km   "
                  f"M0: {np.degrees(M0):.1f}°   "
                  f"elapsed: {elapsed:.0f}s   "
                  f"rate: {rate:.1f} eval/s",
                  end="", flush=True)

        return -count   # minimize negative count

    # --------------------------------------------------
    def _callback(self, xk, convergence):
        """Called by differential evolution after each generation."""
        a, M0 = xk
        gen   = self._eval_count // (self.popsize * 2)
        print(f"\n  Generation {gen:3d}   "
              f"best so far: {self._best_so_far:2d} flybys   "
              f"convergence: {convergence:.4f}   "
              f"a: {a:.0f} km   M0: {np.degrees(M0):.1f}°")

    # --------------------------------------------------
    def optimize(self):
        """
        Run differential evolution to maximize flyby count.
        """
        a_res   = self.resonant_orbit[0]
        a_min   = a_res * (1 - self.a_bounds_pct / 100)
        a_max   = a_res * (1 + self.a_bounds_pct / 100)

        bounds  = [
            (a_min, a_max),   # semi-major axis
            (0, 2 * np.pi),   # mean anomaly at epoch
        ]

        print(f"\nFlyby Optimizer — {self.central_body['name'].capitalize()} "
              f"system → {self.target_body.capitalize()}")
        print(f"Start: {self.start_date}   End: {self.end_date}")
        print(f"Resonant a: {a_res:.0f} km   "
              f"Search: ±{self.a_bounds_pct}%   "
              f"[{a_min:.0f}, {a_max:.0f}] km")
        print(f"Altitude bounds: {self.alt_min_km} – {self.alt_max_km} km")
        print(f"Steps (opt): {self.steps_opt}   "
              f"Steps (final): {self.steps_final}")
        print(f"Population: {self.popsize}   Max iter: {self.maxiter}")
        print(f"{'─'*70}")

        self._t_start    = time.time()
        self._eval_count = 0
        self._best_so_far = 0
        self.history     = []

        result = differential_evolution(
            self._objective,
            bounds=bounds,
            maxiter=self.maxiter,
            popsize=self.popsize,
            tol=0.01,
            seed=self.seed,
            callback=self._callback,
            polish=True,   # local refinement at the end
            workers=1,
        )

        print(f"\n{'─'*70}")
        print(f"Optimization complete in "
              f"{time.time() - self._t_start:.0f}s   "
              f"{self._eval_count} evaluations")

        # --------------------------------------------------
        # final high-resolution propagation
        # --------------------------------------------------
        self.best_a  = result.x[0]
        self.best_M0 = result.x[1]

        print(f"\nFinal propagation at full resolution "
              f"({self.steps_final} steps)...")

        self.final_prop  = self._propagate(
            self.best_a, self.best_M0, self.steps_final
        )
        self.best_flybys = self.find_all_flybys(
            self.final_prop.rs, self.final_prop.times
        )
        self.best_count  = len(self.best_flybys)

        print(f"Final result: {self.best_count} flybys")
        print(f"  a  = {self.best_a:.2f} km  "
              f"(resonant: {a_res:.2f} km, "
              f"Δa = {self.best_a - a_res:+.2f} km)")
        print(f"  M0 = {np.degrees(self.best_M0):.2f}°  "
              f"(resonant: {np.degrees(self.resonant_orbit[3]):.2f}°)")
        print(f"\nFlyby summary:")
        for k, fb in enumerate(self.best_flybys):
            print(f"  [{k+1:2d}]  {fb['time_utc']}   "
                  f"alt = {fb['altitude']:.0f} km")

        return result

    # --------------------------------------------------
    def plot_results(self, save_html=None):
        """
        Plot distance to target body over time with flyby markers,
        and optimization history.
        """
        if self.final_prop is None:
            print("Run optimize() first.")
            return

        prop   = self.final_prop
        times  = prop.times
        ts_days = (np.array(times) - times[0]) / 86400

        # distance to target
        r_body = np.array([
            spice.spkpos(self.target_body, float(t),
                         'J2000', 'NONE',
                         self.central_body['name'].upper())[0]
            for t in times
        ])
        dist = np.linalg.norm(prop.rs - r_body, axis=1) - self.body_radius

        # --------------------------------------------------
        # figure: two rows
        # top: distance over time
        # bottom: optimization history
        # --------------------------------------------------
        fig = make_subplots(
            rows=2, cols=1,
            subplot_titles=(
                f"Distance to {self.target_body.capitalize()} — "
                f"{self.best_count} flybys   "
                f"a={self.best_a:.0f} km   "
                f"M0={np.degrees(self.best_M0):.1f}°",
                "Optimization history"
            ),
            vertical_spacing=0.12,
            row_heights=[0.65, 0.35],
        )

        # distance trace
        fig.add_trace(go.Scatter(
            x=ts_days, y=dist,
            mode='lines',
            line=dict(color='royalblue', width=1),
            name='altitude',
            hovertemplate="day: %{x:.1f}<br>alt: %{y:.0f} km<extra></extra>",
        ), row=1, col=1)

        # altitude bounds
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

        # flyby markers
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

        # optimization history
        if self.history:
            h_evals  = [h[0] for h in self.history]
            h_counts = [h[3] for h in self.history]

            # running best
            running_best = []
            best = 0
            for c in h_counts:
                best = max(best, c)
                running_best.append(best)

            fig.add_trace(go.Scatter(
                x=h_evals, y=h_counts,
                mode='markers',
                marker=dict(color='steelblue', size=3, opacity=0.4),
                name='each evaluation',
                hovertemplate="eval: %{x}<br>flybys: %{y}<extra></extra>",
            ), row=2, col=1)

            fig.add_trace(go.Scatter(
                x=h_evals, y=running_best,
                mode='lines',
                line=dict(color='tomato', width=2),
                name='best so far',
                hovertemplate="eval: %{x}<br>best: %{y}<extra></extra>",
            ), row=2, col=1)

        # --------------------------------------------------
        # layout
        # --------------------------------------------------
        fig.update_layout(
            title=dict(
                text=(
                    f"{self.central_body['name'].capitalize()} system — "
                    f"Flyby optimization: {self.target_body.capitalize()}<br>"
                    f"<sup>{self.start_date} → {self.end_date}   "
                    f"alt bounds: {self.alt_min_km}–{self.alt_max_km} km   "
                    f"result: {self.best_count} flybys</sup>"
                ),
                font=dict(size=13, color='black'),
                x=0.5, xanchor='center',
            ),
            paper_bgcolor='white',
            plot_bgcolor='white',
            font=dict(color='black'),
            width=1100,
            height=800,
            legend=dict(
                x=1.02, y=1,
                bordercolor='lightgrey',
                borderwidth=1,
                font=dict(color='black', size=9),
            ),
        )

        fig.update_xaxes(
            title=dict(text="Time (days)", font=dict(color='black')),
            tickfont=dict(color='black'),
            gridcolor='lightgrey',
            showgrid=True,
        )
        fig.update_yaxes(
            tickfont=dict(color='black'),
            gridcolor='lightgrey',
            showgrid=True,
        )
        fig.update_yaxes(
            title=dict(text="Altitude (km)", font=dict(color='black')),
            row=1, col=1,
        )
        fig.update_yaxes(
            title=dict(text="Flyby count", font=dict(color='black')),
            row=2, col=1,
        )

        if save_html:
            fig.write_html(save_html)
            print(f"Saved to {save_html}")

        fig.show()
        return fig