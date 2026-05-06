import numpy as np
import spiceypy as spice
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from datetime import datetime
from trajectory_utils import lambert_branches, hohmann_phase_angle


def et_to_datetime(et):
    utc_str = spice.et2utc(et, "ISOC", 0)
    return datetime.strptime(utc_str, "%Y-%m-%dT%H:%M:%S")

def et_to_date_str(et):
    return et_to_datetime(et).strftime("%d %b %Y")


class porkchop_interface:

    def __init__(self, origin_body, target_body, launch_date, arrival_date,
                 max_c3, span=90, grid_res=5, mu=1.32712440018e11):

        base_launch  = spice.str2et(launch_date)
        base_arrival = spice.str2et(arrival_date)

        launch_span  = span * 86400
        arrival_span = span * 86400
        step         = grid_res * 86400

        self.mu          = mu
        self.origin_body = origin_body
        self.target_body = target_body
        self.max_c3      = max_c3

        self.launches = np.arange(base_launch - launch_span,
                                  base_launch + launch_span, step)
        self.arrivals = np.arange(base_arrival - arrival_span,
                                  base_arrival + arrival_span, step)

        self.max_tof = (self.arrivals[-1] - self.launches[0]) / 86400

    def compute(self):
        """
        Compute C3, v∞ arrival, TOF and minimum inclination grids.
        Separated from plot() so grids can be reused without recomputing.
        """
        from lambert_interface import lambert_interface
        from trajectory_utils import compute_arrival_inclination

        launches = self.launches
        arrivals = self.arrivals

        n_launch  = len(launches)
        n_arrival = len(arrivals)

        C3      = np.full((n_launch, n_arrival), np.nan)
        VINF    = np.full((n_launch, n_arrival), np.nan)
        TOF     = np.full((n_launch, n_arrival), np.nan)
        INC_MIN = np.full((n_launch, n_arrival), np.nan)

        tested = passed = rejected = 0

        for i, t0 in enumerate(launches):
            for j, t1 in enumerate(arrivals):

                tof_days = (t1 - t0) / 86400
                if tof_days <= 0 or tof_days > self.max_tof:
                    continue

                # phase angle filter
                r1, _ = spice.spkpos(self.origin_body, t0, 'ECLIPJ2000', 'NONE', 'SUN')
                r2, _ = spice.spkpos(self.target_body, t0, 'ECLIPJ2000', 'NONE', 'SUN')

                a1 = np.linalg.norm(r1)
                a2 = np.linalg.norm(r2)

                dot       = np.clip(np.dot(r1, r2) / (a1 * a2), -1, 1)
                phase     = np.degrees(np.arccos(dot))
                cross     = np.cross(r1, r2)
                if cross[2] < 0:
                    phase = 360 - phase

                phi_nom = hohmann_phase_angle(a1, a2, self.mu)
                diff    = abs((phase - phi_nom + 180) % 360 - 180)

                tested += 1
                if diff > 60:
                    rejected += 1
                    continue
                passed += 1

                try:
                    start_str = spice.et2utc(t0, "C", 0)
                    arr_str   = spice.et2utc(t1, "C", 0)

                    sim = lambert_interface(
                        self.origin_body, self.target_body,
                        start_str, arr_str, arr_str,
                        steps=500, perturbations=[]
                    )
                    sim.solve()

                    state_dep, _ = spice.spkezr(self.origin_body, t0,
                                                'ECLIPJ2000', 'NONE', 'SUN')
                    state_arr, _ = spice.spkezr(self.target_body, t1,
                                                'ECLIPJ2000', 'NONE', 'SUN')

                    v_planet_dep = state_dep[3:]
                    v_planet_arr = state_arr[3:]

                    sols = lambert_branches(sim.start_r, sim.end_r,
                                           sim.tof, mu=self.mu)

                    best_c3          = np.inf
                    best_vinf        = np.inf
                    best_vinf_vec    = None

                    for sol in sols:
                        v0 = sol[:3]
                        v1 = sol[3:]

                        v_inf_dep = v0 - v_planet_dep
                        v_inf_arr_vec = v1 - v_planet_arr
                        v_inf_arr = np.linalg.norm(v_inf_arr_vec)
                        c3        = np.dot(v_inf_dep, v_inf_dep)

                        best_c3 = min(best_c3, c3)
                        if v_inf_arr < best_vinf:
                            best_vinf     = v_inf_arr
                            best_vinf_vec = v_inf_arr_vec

                    C3[i, j]   = best_c3
                    VINF[i, j] = best_vinf
                    TOF[i, j]  = tof_days

                    if best_vinf_vec is not None:
                        INC_MIN[i, j] = compute_arrival_inclination(
                            best_vinf_vec, self.target_body, t1
                        )

                except Exception as e:
                    continue

        print(f"tested={tested}, passed={passed}, rejected={rejected}")

        self.C3      = C3
        self.VINF    = VINF
        self.TOF     = TOF
        self.INC_MIN = INC_MIN

    def plot(self, save_html=None):
        """
        Render interactive Plotly porkchop plot.
        Calls compute() if grids not already available.
        """
        if not hasattr(self, 'C3'):
            self.compute()

        launches = self.launches
        arrivals = self.arrivals

        # --------------------------------------------------
        # mask high C3
        # --------------------------------------------------
        C3_plot = np.copy(self.C3)
        C3_plot[C3_plot > self.max_c3] = np.nan
        C3_masked = np.where(np.isfinite(C3_plot), C3_plot, np.nan)

        # --------------------------------------------------
        # date label arrays for axes
        # --------------------------------------------------
        launch_dates  = [et_to_date_str(t) for t in launches]
        arrival_dates = [et_to_date_str(t) for t in arrivals]

        # --------------------------------------------------
        # dynamic contour levels
        # --------------------------------------------------
        c3_min  = np.nanmin(C3_masked)
        c3_max  = np.nanmax(C3_masked)
        c3_step = max(round((c3_max - c3_min) / 12), 1)
        c3_levels = int((c3_max - c3_min) / c3_step)

        vinf_min  = np.nanmin(self.VINF)
        vinf_max  = np.nanmax(self.VINF)
        vinf_step = max(round((vinf_max - vinf_min) / 10, 1), 0.5)

        tof_min   = np.nanmin(self.TOF[np.isfinite(C3_plot)])
        tof_step  = max(round((self.max_tof - tof_min) / 10 / 10) * 10, 10)

        inc_min_val = np.nanmin(self.INC_MIN)
        inc_max_val = np.nanmax(self.INC_MIN)
        inc_step    = max(round((inc_max_val - inc_min_val) / 8), 1)

        # --------------------------------------------------
        # hover text — show all quantities at each grid point
        # --------------------------------------------------
        hover = np.full_like(C3_masked, None, dtype=object)
        for i in range(len(launches)):
            for j in range(len(arrivals)):
                if np.isfinite(C3_masked[i, j]):
                    hover[i, j] = (
                        f"Launch:  {launch_dates[i]}<br>"
                        f"Arrival: {arrival_dates[j]}<br>"
                        f"TOF:     {self.TOF[i,j]:.0f} days<br>"
                        f"C3:      {self.C3[i,j]:.1f} km²/s²<br>"
                        f"v∞ arr:  {self.VINF[i,j]:.2f} km/s<br>"
                        f"inc min: {self.INC_MIN[i,j]:.1f}°"
                    )

        # --------------------------------------------------
        # build figure
        # --------------------------------------------------
        fig = go.Figure()

        # C3 heatmap
        fig.add_trace(go.Heatmap(
            z=C3_masked.T,
            x=launch_dates,
            y=arrival_dates,
            colorscale='Viridis_r',
            colorbar=dict(title="C3 (km²/s²)", x=1.02),
            hoverinfo='text',
            text=hover.T,
            name='C3',
            zmin=c3_min,
            zmax=c3_max,
        ))

        # C3 contours
        fig.add_trace(go.Contour(
            z=C3_masked.T,
            x=launch_dates,
            y=arrival_dates,
            colorscale=[[0, 'yellow'], [1, 'yellow']],
            showscale=False,
            contours=dict(
                start=c3_min,
                end=c3_max,
                size=c3_step,
                showlabels=True,
                labelfont=dict(size=10, color='white'),
                coloring='lines'
            ),
            line=dict(width=0.8),
            hoverinfo='skip',
            name='C3 contours',
        ))

        # v∞ arrival contours
        fig.add_trace(go.Contour(
            z=self.VINF.T,
            x=launch_dates,
            y=arrival_dates,
            colorscale=[[0, 'grey'], [1, 'grey']],
            showscale=False,
            contours=dict(
                start=vinf_min,
                end=vinf_max,
                size=vinf_step,
                showlabels=True,
                labelfont=dict(size=10, color='lightgrey'),
                coloring='lines'
            ),
            line=dict(width=1.0, dash='dot'),
            hoverinfo='skip',
            name='v∞ arrival',
        ))

        # TOF contours
        fig.add_trace(go.Contour(
            z=self.TOF.T,
            x=launch_dates,
            y=arrival_dates,
            colorscale=[[0, 'red'], [1, 'red']],
            showscale=False,
            contours=dict(
                start=tof_min,
                end=self.max_tof,
                size=tof_step,
                showlabels=True,
                labelfont=dict(size=10, color='red'),
                coloring='lines'
            ),
            line=dict(width=1.0, dash='dash'),
            hoverinfo='skip',
            name='TOF',
        ))

        # inclination contours
        fig.add_trace(go.Contour(
            z=self.INC_MIN.T,
            x=launch_dates,
            y=arrival_dates,
            colorscale=[[0, 'white'], [1, 'white']],
            showscale=False,
            contours=dict(
                start=inc_min_val,
                end=inc_max_val,
                size=inc_step,
                showlabels=True,
                labelfont=dict(size=10, color='white', shadow='1px 1px 2px black, -1px -1px 2px black'),
                coloring='lines'
            ),
            line=dict(width=0.8),
            hoverinfo='skip',
            name='min inclination',
        ))

        # --------------------------------------------------
        # layout
        # --------------------------------------------------
        start_year = et_to_datetime(launches[0]).year
        end_year   = et_to_datetime(arrivals[-1]).year

        fig.update_layout(
            title=dict(
                text=f"{self.origin_body} → {self.target_body} porkchop plot<br>"
                     f"<sup>Launch window: {launch_dates[0]} – {launch_dates[-1]}   "
                     f"Arrival window: {arrival_dates[0]} – {arrival_dates[-1]}</sup>",
                font=dict(size=16, color='black'),
                x=0.5,
                xanchor='center'
            ),
            xaxis=dict(
                title="Launch Date",
                tickangle=30,
                tickmode='array',
                tickvals=launch_dates[::max(1, len(launch_dates)//10)],   # ~10 ticks
                ticktext=launch_dates[::max(1, len(launch_dates)//10)],
                tickfont=dict(size=12, color='black'),
            ),
            yaxis=dict(
                title="Arrival Date",
                tickmode='array',
                tickvals=arrival_dates[::max(1, len(arrival_dates)//10)],  # ~10 ticks
                ticktext=arrival_dates[::max(1, len(arrival_dates)//10)],
                tickfont=dict(size=12, color='black'),
            ),
            width=1500,
            height=1050,
            paper_bgcolor='white',
            plot_bgcolor="#ffffff",
            font=dict(color='black'),
            legend=dict(
                x=1.08, y=1,
                bgcolor='rgba(0,0,0,0.5)',
                bordercolor='grey'
            )
        )

        if save_html:
            fig.write_html(save_html)
            print(f"Saved to {save_html}")

        fig.show()