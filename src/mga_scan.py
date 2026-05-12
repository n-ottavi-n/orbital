import numpy as np
import spiceypy as spice
import plotly.graph_objects as go
import planetary_data as pd
from lambert_interface import lambert_interface
from trajectory_utils import lambert_branches, hohmann_phase_angle
from datetime import datetime
import csv
import time
import pandas as pd_df

AU = 1.495978707e8  # km per AU


def et_to_utc(et):
    return spice.et2utc(float(et), "C", 0)

def et_to_datetime(et):
    utc_str = spice.et2utc(float(et), "ISOC", 0)
    return datetime.strptime(utc_str, "%Y-%m-%dT%H:%M:%S")


class MGAscan:

    def __init__(self,
                 sequence,
                 launch_window,
                 flyby_windows,
                 propagation_days,
                 vinf_bounds,
                 periapsis_radii,
                 min_flyby_alt_km,
                 grid_res_days=5,
                 mu_sun=1.32712440018e11):
        """
        Parameters
        ----------
        sequence         : list of body name strings
                           e.g. ["EARTH", "VENUS", "JUPITER"]
        launch_window    : (start_date_str, end_date_str)
                           e.g. ("DEC 01 2030", "MAR 01 2031")
        flyby_windows    : list of (start_date_str, end_date_str)
                           one per intermediate body
                           e.g. [("MAR 01 2031", "JUN 01 2031")]
        propagation_days : how long to propagate after last flyby
                           to check for target SOI crossing
        vinf_bounds      : dict of (min, max) v∞ per body in km/s
                           e.g. {"EARTH": (0, 6), "VENUS": (4, 8),
                                 "JUPITER": (0, 10)}
        periapsis_radii  : list of flyby periapsis distances in body radii
                           e.g. [1.05, 2, 3, 5, 10, 20, 30]
        min_flyby_alt_km : dict of minimum flyby altitude per body in km
                           e.g. {"VENUS": 300}
        grid_res_days    : grid resolution in days
        mu_sun           : gravitational parameter of Sun (km³/s²)
        """

        self.sequence          = sequence
        self.origin            = sequence[0]
        self.target            = sequence[-1]
        self.flyby_bodies      = sequence[1:-1]
        self.propagation_days  = propagation_days
        self.vinf_bounds       = vinf_bounds
        self.periapsis_radii   = periapsis_radii
        self.min_flyby_alt_km  = min_flyby_alt_km
        self.grid_res_days     = grid_res_days
        self.mu_sun            = mu_sun

        # --------------------------------------------------
        # convert window strings to ET arrays
        # --------------------------------------------------
        step = grid_res_days * 86400

        t0, t1          = spice.str2et(launch_window[0]), spice.str2et(launch_window[1])
        self.launches   = np.arange(t0, t1, step)

        self.flyby_grids = []
        for window in flyby_windows:
            t0, t1 = spice.str2et(window[0]), spice.str2et(window[1])
            self.flyby_grids.append(np.arange(t0, t1, step))

        # --------------------------------------------------
        # body data
        # --------------------------------------------------
        self.body_data = {}
        for body in sequence:
            data           = pd.get_body(body)
            radius         = pd.get_body_radius(body)
            state, _       = spice.spkezr(body, 0.0, "ECLIPJ2000", "NONE", "SUN")
            r_mag          = np.linalg.norm(state[:3])
            v_mag          = np.linalg.norm(state[3:])
            a_km           = 1.0 / (2.0/r_mag - v_mag**2/mu_sun)
            mu_body        = data['mu']

            # SOI radius
            r_soi = a_km * (mu_body / mu_sun) ** (2.0/5.0)

            self.body_data[body] = {
                'mu':     mu_body,
                'radius': radius,
                'a_km':   a_km,
                'r_soi':  r_soi,
            }

        # results populated by compute()
        self.candidates = []

    # --------------------------------------------------
    def _phase_angle_filter(self, origin, dest, t0):
        """Quick geometric filter — skip if phase angle too far from Hohmann."""
        r1, _ = spice.spkpos(origin, t0, "ECLIPJ2000", "NONE", "SUN")
        r2, _ = spice.spkpos(dest,   t0, "ECLIPJ2000", "NONE", "SUN")

        a1 = np.linalg.norm(r1)
        a2 = np.linalg.norm(r2)

        dot   = np.clip(np.dot(r1, r2) / (a1 * a2), -1, 1)
        phase = np.degrees(np.arccos(dot))
        cross = np.cross(r1, r2)
        if cross[2] < 0:
            phase = 360 - phase

        phi_nom = hohmann_phase_angle(a1, a2, self.mu_sun)
        diff    = abs((phase - phi_nom + 180) % 360 - 180)
        return diff <= 90   # wider tolerance than direct porkchop

    # --------------------------------------------------
    def _keplerian_propagate(self, r0, v0, dt):
        """
        Propagate state (r0, v0) forward by dt seconds
        using universal variable method.
        Returns (r, v) at t0 + dt.
        Pure Keplerian — no perturbations.
        Uses the existing lambert_universal infrastructure.
        """
        from scipy.optimize import brentq

        mu  = self.mu_sun
        r0  = np.array(r0)
        v0  = np.array(v0)

        r0m = np.linalg.norm(r0)
        v0m = np.linalg.norm(v0)

        # specific energy and semi-major axis
        energy = v0m**2 / 2 - mu / r0m
        if abs(energy) < 1e-10:
            # parabolic — rare, skip
            return r0, v0
        a = -mu / (2 * energy)

        # Lagrange coefficients via universal variable
        alpha = 1.0 / a
        sigma0 = np.dot(r0, v0) / np.sqrt(mu)

        def kepler_eq(chi):
            psi   = chi**2 * alpha
            c2    = self._stumpff_c2(psi)
            c3    = self._stumpff_c3(psi)
            r_chi = (r0m * (1 - psi * c2) +
                     sigma0 * chi * (1 - psi * c3) +
                     chi**2 * c2)
            return (r0m * chi * (1 - psi * c3) +
                    sigma0 * chi**2 * c2 +
                    chi**3 * c3 - np.sqrt(mu) * dt)

        # initial guess for chi
        chi0 = np.sqrt(mu) * abs(alpha) * dt

        try:
            chi = brentq(kepler_eq, -10 * abs(chi0), 10 * abs(chi0),
                         xtol=1e-8, maxiter=100)
        except Exception:
            return r0, v0

        psi = chi**2 * alpha
        c2  = self._stumpff_c2(psi)
        c3  = self._stumpff_c3(psi)

        r_chi = (r0m * (1 - psi * c2) +
                 sigma0 * chi * (1 - psi * c3) +
                 chi**2 * c2)

        f    = 1 - chi**2 * c2 / r0m
        g    = dt - chi**3 * c3 / np.sqrt(mu)
        gdot = 1 - chi**2 * c2 / r_chi
        fdot = np.sqrt(mu) * chi * (psi * c3 - 1) / (r_chi * r0m)

        r = f * r0 + g * v0
        v = fdot * r0 + gdot * v0

        return r, v

    # --------------------------------------------------
    def _stumpff_c2(self, psi):
        if psi > 1e-6:
            return (1 - np.cos(np.sqrt(psi))) / psi
        elif psi < -1e-6:
            return (1 - np.cosh(np.sqrt(-psi))) / psi
        return 0.5

    def _stumpff_c3(self, psi):
        if psi > 1e-6:
            sp = np.sqrt(psi)
            return (sp - np.sin(sp)) / psi**1.5
        elif psi < -1e-6:
            sp = np.sqrt(-psi)
            return (np.sinh(sp) - sp) / (-psi)**1.5
        return 1.0/6.0

    # --------------------------------------------------
    def _rotate_vinf(self, v_inf_in, turn_angle):
        """
        Rotate v_inf_in by turn_angle in the ecliptic plane.
        Fixed clock angle — coplanar assumption.
        Returns rotated v∞ vector.
        """
        # rotation axis = ecliptic north pole
        k = np.array([0.0, 0.0, 1.0])

        # Rodrigues rotation formula
        v     = v_inf_in
        v_rot = (v * np.cos(turn_angle) +
                 np.cross(k, v) * np.sin(turn_angle) +
                 k * np.dot(k, v) * (1 - np.cos(turn_angle)))
        return v_rot

    # --------------------------------------------------
    def _check_target_soi(self, r_sc, v_sc, t_flyby):
        """
        Analytically check if post-flyby orbit reaches target SOI.
        """
        target     = self.target
        r_soi      = self.body_data[target]['r_soi']
        a_target   = self.body_data[target]['a_km']
        mu         = self.mu_sun

        r = np.array(r_sc)
        v = np.array(v_sc)

        rmag = np.linalg.norm(r)
        vmag = np.linalg.norm(v)

        # orbital energy and semi-major axis
        energy = vmag**2 / 2 - mu / rmag

        a_km  = -mu / (2 * energy)
        h  = np.cross(r, v)
        hmag = np.linalg.norm(h)
        e_vec = np.cross(v, h) / mu - r / rmag
        e = np.linalg.norm(e_vec)

        p  = hmag**2 / mu
        rp = p / (1 + e)
        ra = p / (1 - e) if e < 1 else np.inf   # inf for hyperbola

        # does apoapsis reach target's orbit?
        if ra < a_target - r_soi:
            return None   # orbit doesn't reach target

        # estimate TOF to target distance using Kepler's equation
        # find true anomaly where r = a_target (approximately)
        p = a_km * (1 - e**2)

        # true anomaly at target distance
        cos_nu = np.clip((p / a_target - 1) / e, -1, 1)
        nu_target = np.arccos(cos_nu)

        # current true anomaly
        e_vec = np.cross(v, h) / mu - r / rmag
        cos_nu0 = np.clip(np.dot(e_vec / np.linalg.norm(e_vec), r / rmag), -1, 1)
        nu0 = np.arccos(cos_nu0)
        if np.dot(r, v) < 0:
            nu0 = 2 * np.pi - nu0

        # TOF calculation — handle both elliptic and hyperbolic
        if e < 1:
            # elliptic — eccentric anomaly
            def nu_to_E(nu):
                return 2 * np.arctan2(np.sqrt(1-e) * np.sin(nu/2),
                                    np.sqrt(1+e) * np.cos(nu/2))
            
            E0       = nu_to_E(nu0)
            E_target = nu_to_E(nu_target)
            n        = np.sqrt(mu / abs(a_km)**3)
            M0       = E0 - e * np.sin(E0)
            M_target = E_target - e * np.sin(E_target)
            if M_target < M0:
                M_target += 2 * np.pi
            tof = (M_target - M0) / n

        else:
            # hyperbolic — hyperbolic anomaly
            def nu_to_F(nu):
                return 2 * np.arctanh(np.sqrt((e-1)/(e+1)) * np.tan(nu/2))
            
            F0       = nu_to_F(nu0)
            F_target = nu_to_F(nu_target)
            n        = np.sqrt(mu / abs(a_km)**3)
            tof      = (e * np.sinh(F_target) - F_target -
                        e * np.sinh(F0) + F0) / n

        # check TOF is within propagation window
        if tof > self.propagation_days * 86400:
            return None

        t_arrival = t_flyby + tof

        # get target state at estimated arrival
        state_target, _ = spice.spkezr(target, t_arrival,
                                        "ECLIPJ2000", "NONE", "SUN")
        r_target = state_target[:3]
        v_target = state_target[3:]

        # propagate spacecraft to arrival time
        r_arr, v_arr = self._keplerian_propagate(r, v, tof)

        # check actual closest approach distance
        range_km = np.linalg.norm(r_arr - r_target)

        if range_km > r_soi * 1:   # generous check — within SOI
            return None

        # compute v∞ at target
        v_inf_vec   = v_arr - v_target
        vinf_target = np.linalg.norm(v_inf_vec)

        # check v∞ bounds
        vmin, vmax = self.vinf_bounds.get(target, (0, np.inf))
        if not (vmin <= vinf_target <= vmax):
            return None

        return {
            "t_closest":        t_arrival,
            "min_range_km":     range_km,
            "soi_fraction":     range_km / r_soi,
            "vinf_target_km_s": vinf_target,
            "vinf_target_vec":  v_inf_vec,
        }

    # --------------------------------------------------
    def compute(self, output_file=None):
        """
        Run the MGA scan.
        Populates self.candidates with all valid solutions.
        Optionally writes results to CSV file.
        """

        self.candidates = []

        # CSV setup
        if output_file:
            csvfile = open(output_file, 'w', newline='', encoding='utf-8')
            writer  = self._init_csv(csvfile)
        else:
            csvfile = None
            writer  = None

        total    = len(self.launches) * len(self.flyby_grids[0])
        tested   = 0
        found    = 0

        print(f"MGA scan: {' → '.join(self.sequence)}")
        print(f"Grid: {len(self.launches)} launch × "
              f"{len(self.flyby_grids[0])} flyby = {total} combinations")
        print(f"Periapsis radii: {self.periapsis_radii}")

        t_start  = time.time()
        total    = len(self.launches) * len(self.flyby_grids[0])
        tested   = 0
        found    = 0

        for i, t_launch in enumerate(self.launches):

            # phase angle filter on Leg 1
            if not self._phase_angle_filter(self.origin,
                                            self.flyby_bodies[0],
                                            t_launch):
                continue

            for j, t_flyby in enumerate(self.flyby_grids[0]):

                # must arrive after departure
                tof_leg1 = (t_flyby - t_launch) / 86400
                if tof_leg1 <= 0:
                    continue

                tested += 1

                if tested % 50 == 0:
                    elapsed   = time.time() - t_start
                    rate      = tested / elapsed           # combinations per second
                    remaining = (total - tested) / rate    # seconds
                    
                    print(f"\r  {tested}/{total} ({100*tested/total:.1f}%)   "
                        f"found: {found}   "
                        f"elapsed: {elapsed:.0f}s   "
                        f"remaining: {remaining:.0f}s   ",
                        end="", flush=True)

                # ------------------------------------------
                # Leg 1: Lambert origin → flyby body
                # ------------------------------------------
                try:
                    start_str  = et_to_utc(t_launch)
                    flyby_str  = et_to_utc(t_flyby)

                    sim = lambert_interface(
                        self.origin, self.flyby_bodies[0],
                        start_str, flyby_str, flyby_str,
                        steps=200, perturbations=[]
                    )
                    sim.solve()

                    # departure v∞
                    state_dep, _ = spice.spkezr(self.origin, t_launch,
                                                "ECLIPJ2000", "NONE", "SUN")
                    v_dep        = state_dep[3:]
                    v_inf_dep    = sim.v0 - v_dep
                    vinf_dep     = np.linalg.norm(v_inf_dep)
                    c3           = vinf_dep**2

                    # check departure v∞ bounds
                    vmin, vmax = self.vinf_bounds.get(self.origin, (0, np.inf))
                    if not (vmin <= vinf_dep <= vmax):
                        continue

                    # arrival v∞ at flyby body
                    state_flyby, _ = spice.spkezr(self.flyby_bodies[0],
                                                   t_flyby,
                                                   "ECLIPJ2000", "NONE", "SUN")
                    r_flyby_body   = state_flyby[:3]
                    v_flyby_body   = state_flyby[3:]
                    v_arr_leg1     = sim.v[3:]   # Lambert arrival velocity
                    v_inf_in       = v_arr_leg1 - v_flyby_body
                    vinf_flyby     = np.linalg.norm(v_inf_in)

                    # check flyby v∞ bounds
                    flyby_body = self.flyby_bodies[0]
                    vmin, vmax = self.vinf_bounds.get(flyby_body, (0, np.inf))
                    if not (vmin <= vinf_flyby <= vmax):
                        continue

                except Exception:
                    continue

                # ------------------------------------------
                # turning angle sweep at flyby body
                # ------------------------------------------
                body_radius = self.body_data[flyby_body]['radius']
                mu_flyby    = self.body_data[flyby_body]['mu']
                min_alt     = self.min_flyby_alt_km.get(flyby_body, 300)
                r_min       = body_radius + min_alt

                # spacecraft heliocentric state at flyby
                r_sc_flyby = sim.end_r   # position at flyby (from Lambert)

                for r_mult in self.periapsis_radii:

                    r_peri = max(r_mult * body_radius, r_min)

                    # maximum turn angle at this periapsis
                    sin_half = mu_flyby / (mu_flyby + r_peri * vinf_flyby**2)
                    if sin_half > 1:
                        continue
                    turn_max = 2 * np.arcsin(sin_half)

                    # both prograde and retrograde turns
                    for sign in [+1, -1]:
                        turn_angle = sign * turn_max


                        # rotate v∞_in to get v∞_out
                        v_inf_out = self._rotate_vinf(v_inf_in, turn_angle)

                        # heliocentric departure velocity after flyby
                        v_sc_out  = v_flyby_body + v_inf_out

                        # ------------------------------------------
                        # propagate Leg 2 and check target SOI
                        # ------------------------------------------
                        result = self._check_target_soi(
                            r_sc_flyby, v_sc_out, t_flyby
                        )

                        if result is None:
                            continue

                        found += 1

                        candidate = {
                            # sequence
                            "sequence":            " → ".join(self.sequence),

                            # departure
                            "t_launch_et":         t_launch,
                            "t_launch_utc":        et_to_utc(t_launch),
                            "vinf_dep_km_s":       vinf_dep,
                            "c3_km2_s2":           c3,
                            "vinf_dep_vec_x":      v_inf_dep[0],
                            "vinf_dep_vec_y":      v_inf_dep[1],
                            "vinf_dep_vec_z":      v_inf_dep[2],

                            # flyby
                            "flyby_body":          flyby_body,
                            "t_flyby_et":          t_flyby,
                            "t_flyby_utc":         et_to_utc(t_flyby),
                            "tof_leg1_days":       tof_leg1,
                            "vinf_flyby_km_s":     vinf_flyby,
                            "flyby_periapsis_km":  r_peri,
                            "flyby_alt_km":        r_peri - body_radius,
                            "turn_angle_deg":      np.degrees(turn_angle),
                            "vinf_in_vec_x":       v_inf_in[0],
                            "vinf_in_vec_y":       v_inf_in[1],
                            "vinf_in_vec_z":       v_inf_in[2],
                            "vinf_out_vec_x":      v_inf_out[0],
                            "vinf_out_vec_y":      v_inf_out[1],
                            "vinf_out_vec_z":      v_inf_out[2],

                            # arrival
                            "target_body":         self.target,
                            "t_arrival_et":        result["t_closest"],
                            "t_arrival_utc":       et_to_utc(result["t_closest"]),
                            "tof_leg2_days":       (result["t_closest"] - t_flyby) / 86400,
                            "tof_total_days":      (result["t_closest"] - t_launch) / 86400,
                            "vinf_target_km_s":    result["vinf_target_km_s"],
                            "min_range_km":        result["min_range_km"],
                            "soi_fraction":        result["soi_fraction"],
                            "vinf_target_vec_x":   result["vinf_target_vec"][0],
                            "vinf_target_vec_y":   result["vinf_target_vec"][1],
                            "vinf_target_vec_z":   result["vinf_target_vec"][2],
                        }

                        self.candidates.append(candidate)

                        if writer:
                            writer.writerow(candidate)

                if tested % 100 == 0:
                    print(f"  {tested}/{total} tested, {found} candidates found")

        print(f"\nScan complete: {tested} tested, {found} candidates found")

        if csvfile:
            csvfile.close()
            print(f"Results saved to {output_file}")

    # --------------------------------------------------
    def _init_csv(self, csvfile):
        """Initialize CSV writer with header."""
        fieldnames = [
            "sequence",
            "t_launch_et", "t_launch_utc",
            "vinf_dep_km_s", "c3_km2_s2",
            "vinf_dep_vec_x", "vinf_dep_vec_y", "vinf_dep_vec_z",
            "flyby_body",
            "t_flyby_et", "t_flyby_utc", "tof_leg1_days",
            "vinf_flyby_km_s", "flyby_periapsis_km", "flyby_alt_km",
            "turn_angle_deg",
            "vinf_in_vec_x", "vinf_in_vec_y", "vinf_in_vec_z",
            "vinf_out_vec_x", "vinf_out_vec_y", "vinf_out_vec_z",
            "target_body",
            "t_arrival_et", "t_arrival_utc",
            "tof_leg2_days", "tof_total_days",
            "vinf_target_km_s", "min_range_km", "soi_fraction",
            "vinf_target_vec_x", "vinf_target_vec_y", "vinf_target_vec_z",
        ]

        # write metadata header as comments
        csvfile.write(f"# MGA SCAN\n")
        csvfile.write(f"# Sequence: {' → '.join(self.sequence)}\n")
        csvfile.write(f"# v∞ bounds: {self.vinf_bounds}\n")
        csvfile.write(f"# Periapsis radii: {self.periapsis_radii}\n")
        csvfile.write(f"# Min flyby alt: {self.min_flyby_alt_km}\n")
        csvfile.write(f"# Grid resolution: {self.grid_res_days} days\n")
        csvfile.write(f"# Propagation: {self.propagation_days} days\n#\n")

        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        return writer

    # --------------------------------------------------
    def plot(self, save_html=None):
        """
        Plot MGA scan as scatter plot.
        x: launch date
        y: target arrival date
        color: v∞ at target
        size: inversely proportional to C3 (larger = cheaper departure)
        """
        if not self.candidates:
            print("No candidates — run compute() first.")
            return

        
        df  = pd_df.DataFrame(self.candidates)
        fig = self._build_scatter_figure(df)

        if save_html:
            fig.write_html(save_html)
            print(f"Saved to {save_html}")

        fig.show()
        return fig
    
    def launch_dashboard(self):
        """
        Launch an interactive dashboard:
        - left panel: MGA scatter plot
        - right panel: heliocentric trajectory of clicked candidate
        """
        from dash import Dash, dcc, html, Input, Output, no_update
        import plotly.graph_objects as go
        import pandas as pd_df

        app  = Dash(__name__)
        df   = pd_df.DataFrame(self.candidates)

        # build the scatter figure (reuse existing plot logic)
        scatter_fig = self._build_scatter_figure(df)

        app.layout = html.Div([

            html.H2(
                f"MGA Scan — {' → '.join(self.sequence)}",
                style={'textAlign': 'center', 'fontFamily': 'Arial',
                    'color': 'black'}
            ),

            html.Div([

                # left: scatter plot
                html.Div([
                    dcc.Graph(
                        id='scan-scatter',
                        figure=scatter_fig,
                        style={'height': '750px'},
                        clear_on_unhover=True,
                    )
                ], style={'width': '48%', 'display': 'inline-block',
                        'verticalAlign': 'top'}),

                # right: candidate trajectory
                html.Div([
                    dcc.Graph(
                        id='candidate-plot',
                        figure=go.Figure(),   # empty until clicked
                        style={'height': '750px'},
                    )
                ], style={'width': '50%', 'display': 'inline-block',
                        'verticalAlign': 'top', 'marginLeft': '1%'}),

            ]),

            # info bar below
            html.Div(
                id='candidate-info',
                style={'textAlign': 'center', 'fontFamily': 'Arial',
                    'fontSize': '13px', 'color': '#333',
                    'marginTop': '10px', 'padding': '8px',
                    'backgroundColor': '#f5f5f5',
                    'borderRadius': '4px'}
            ),

        ], style={'backgroundColor': 'white', 'padding': '20px'})

        # --------------------------------------------------
        # callback: click on scatter → update trajectory plot
        # --------------------------------------------------
        @app.callback(
            Output('candidate-plot', 'figure'),
            Output('candidate-info', 'children'),
            Input('scan-scatter', 'clickData'),
        )
        def update_trajectory(clickData):
            if clickData is None:
                return go.Figure(), "Click a candidate to see its trajectory."

            # recover candidate index from customdata
            idx = clickData['points'][0].get('customdata')
            if idx is None:
                return no_update, no_update

            c   = self.candidates[idx]
            fig = self.plot_candidate(c, _return_only=True)

            info = (
                f"Candidate {idx}   |   "
                f"Launch: {c['t_launch_utc']}   |   "
                f"Flyby: {c['t_flyby_utc']}   |   "
                f"Arrival: {c['t_arrival_utc']}   |   "
                f"C3 = {c['c3_km2_s2']:.1f} km²/s²   |   "
                f"v∞ target = {c['vinf_target_km_s']:.2f} km/s   |   "
                f"TOF = {c['tof_total_days']:.0f} days"
            )

            return fig, info

        print("Dashboard running at http://127.0.0.1:8050")
        app.run(debug=False)
    
    def plot_candidate(self, candidate, save_html=None, n_points=500, _return_only=False):
        """
        Plot heliocentric trajectory for a single MGA candidate.

        Parameters
        ----------
        candidate : int or dict
            Index into self.candidates, or a dict/pandas Series
            with the candidate data (e.g. a row from the CSV).
        save_html : str, optional
            Path to save interactive HTML.
        n_points  : int
            Number of points per trajectory leg.
        """
        import pandas as pd_df

        if isinstance(candidate, int):
            c = self.candidates[candidate]
        elif hasattr(candidate, 'to_dict'):
            c = candidate.to_dict()   # pandas Series
        else:
            c = candidate             # already a dict

        mu  = self.mu_sun
        AU  = 1.495978707e8

        # --------------------------------------------------
        # reconstruct departure state
        # --------------------------------------------------
        t_launch = float(c["t_launch_et"])
        t_flyby  = float(c["t_flyby_et"])
        t_arrival= float(c["t_arrival_et"])

        # Earth state at launch
        state_origin, _ = spice.spkezr(
            self.origin, t_launch, "ECLIPJ2000", "NONE", "SUN"
        )
        r_launch = state_origin[:3]
        v_origin = state_origin[3:]

        # departure heliocentric velocity = v_planet + v∞_dep
        v_inf_dep = np.array([
            c["vinf_dep_vec_x"],
            c["vinf_dep_vec_y"],
            c["vinf_dep_vec_z"]
        ])
        v_launch = v_origin + v_inf_dep

        # --------------------------------------------------
        # reconstruct post-flyby state
        # --------------------------------------------------
        state_flyby_body, _ = spice.spkezr(
            self.flyby_bodies[0], t_flyby, "ECLIPJ2000", "NONE", "SUN"
        )
        r_flyby      = state_flyby_body[:3]
        v_flyby_body = state_flyby_body[3:]

        v_inf_out = np.array([
            c["vinf_out_vec_x"],
            c["vinf_out_vec_y"],
            c["vinf_out_vec_z"]
        ])
        v_postflyby = v_flyby_body + v_inf_out

        # --------------------------------------------------
        # propagate leg 1: launch → flyby
        # --------------------------------------------------
        tof_leg1  = t_flyby - t_launch
        dt1       = tof_leg1 / n_points
        leg1      = np.zeros((n_points + 1, 3))
        leg1[0]   = r_launch
        r, v      = r_launch.copy(), v_launch.copy()

        for k in range(n_points):
            r, v = self._keplerian_propagate(r, v, dt1)
            leg1[k + 1] = r

        # --------------------------------------------------
        # propagate leg 2: flyby → arrival
        # --------------------------------------------------
        tof_leg2  = t_arrival - t_flyby
        dt2       = tof_leg2 / n_points
        leg2      = np.zeros((n_points + 1, 3))
        leg2[0]   = r_flyby
        r, v      = r_flyby.copy(), v_postflyby.copy()

        for k in range(n_points):
            r, v = self._keplerian_propagate(r, v, dt2)
            leg2[k + 1] = r

        # --------------------------------------------------
        # planet positions over full mission duration
        # --------------------------------------------------
        times_full = np.linspace(t_launch, t_arrival,
                                n_points * 2)

        def get_body_track(body_name):
            return np.array([
                spice.spkpos(body_name, float(t),
                            "ECLIPJ2000", "NONE", "SUN")[0]
                for t in times_full
            ]) / AU

        origin_track  = get_body_track(self.origin)
        flyby_track   = get_body_track(self.flyby_bodies[0])
        target_track  = get_body_track(self.target)

        # key positions
        r_dep_planet = np.array(
            spice.spkpos(self.origin, t_launch,
                        "ECLIPJ2000", "NONE", "SUN")[0]
        ) / AU

        r_flyby_planet = np.array(
            spice.spkpos(self.flyby_bodies[0], t_flyby,
                        "ECLIPJ2000", "NONE", "SUN")[0]
        ) / AU

        r_arr_planet = np.array(
            spice.spkpos(self.target, t_arrival,
                        "ECLIPJ2000", "NONE", "SUN")[0]
        ) / AU

        # convert trajectories to AU
        leg1_au = leg1 / AU
        leg2_au = leg2 / AU

        # --------------------------------------------------
        # figure
        # --------------------------------------------------
        fig = go.Figure()

        # Sun
        fig.add_trace(go.Scatter3d(
            x=[0], y=[0], z=[0],
            mode='markers',
            marker=dict(color='gold', size=8),
            name='Sun',
            hoverinfo='name',
        ))

        # planet tracks
        for track, name, color in [
            (origin_track,  self.origin.capitalize(),          'steelblue'),
            (flyby_track,   self.flyby_bodies[0].capitalize(), 'tomato'),
            (target_track,  self.target.capitalize(),          'sandybrown'),
        ]:
            fig.add_trace(go.Scatter3d(
                x=track[:, 0], y=track[:, 1], z=track[:, 2],
                mode='lines',
                line=dict(color=color, width=1),
                opacity=0.4,
                name=f"{name} orbit",
                hoverinfo='skip',
            ))

        # leg 1
        fig.add_trace(go.Scatter3d(
            x=leg1_au[:, 0], y=leg1_au[:, 1], z=leg1_au[:, 2],
            mode='lines',
            line=dict(color='royalblue', width=3),
            name=f"Leg 1: {self.origin.capitalize()} → "
                f"{self.flyby_bodies[0].capitalize()}",
            hoverinfo='skip',
        ))

        # leg 2
        fig.add_trace(go.Scatter3d(
            x=leg2_au[:, 0], y=leg2_au[:, 1], z=leg2_au[:, 2],
            mode='lines',
            line=dict(color='darkorange', width=3),
            name=f"Leg 2: {self.flyby_bodies[0].capitalize()} → "
                f"{self.target.capitalize()}",
            hoverinfo='skip',
        ))

        # departure point
        fig.add_trace(go.Scatter3d(
            x=[r_dep_planet[0]], y=[r_dep_planet[1]], z=[r_dep_planet[2]],
            mode='markers',
            marker=dict(color='steelblue', size=6, symbol='circle'),
            name=f"Departure  {c['t_launch_utc']}",
            hovertemplate=(
                f"<b>Departure</b><br>"
                f"{c['t_launch_utc']}<br>"
                f"C3 = {c['c3_km2_s2']:.1f} km²/s²<br>"
                f"v∞ = {c['vinf_dep_km_s']:.2f} km/s"
                "<extra></extra>"
            ),
        ))

        # flyby point
        fig.add_trace(go.Scatter3d(
            x=[r_flyby_planet[0]], y=[r_flyby_planet[1]], z=[r_flyby_planet[2]],
            mode='markers',
            marker=dict(color='tomato', size=6, symbol='diamond'),
            name=f"Flyby  {c['t_flyby_utc']}",
            hovertemplate=(
                f"<b>{self.flyby_bodies[0].capitalize()} flyby</b><br>"
                f"{c['t_flyby_utc']}<br>"
                f"v∞ = {c['vinf_flyby_km_s']:.2f} km/s<br>"
                f"alt = {c['flyby_alt_km']:.0f} km<br>"
                f"turn = {c['turn_angle_deg']:.1f}°"
                "<extra></extra>"
            ),
        ))

        # arrival point
        fig.add_trace(go.Scatter3d(
            x=[r_arr_planet[0]], y=[r_arr_planet[1]], z=[r_arr_planet[2]],
            mode='markers',
            marker=dict(color='sandybrown', size=6, symbol='square'),
            name=f"Arrival  {c['t_arrival_utc']}",
            hovertemplate=(
                f"<b>{self.target.capitalize()} arrival</b><br>"
                f"{c['t_arrival_utc']}<br>"
                f"v∞ = {c['vinf_target_km_s']:.2f} km/s<br>"
                f"TOF = {c['tof_total_days']:.0f} days"
                "<extra></extra>"
            ),
        ))

        # --------------------------------------------------
        # layout
        # --------------------------------------------------
        seq_str = " → ".join(self.sequence)

        fig.update_layout(
            title=dict(
                text=(
                    f"{seq_str} trajectory<br>"
                    f"<sup>Launch: {c['t_launch_utc']}   "
                    f"Flyby: {c['t_flyby_utc']}   "
                    f"Arrival: {c['t_arrival_utc']}<br>"
                    f"C3 = {c['c3_km2_s2']:.1f} km²/s²   "
                    f"v∞ flyby = {c['vinf_flyby_km_s']:.2f} km/s   "
                    f"v∞ target = {c['vinf_target_km_s']:.2f} km/s   "
                    f"TOF = {c['tof_total_days']:.0f} days</sup>"
                ),
                font=dict(size=13, color='black'),
                x=0.5,
                xanchor='center',
            ),
            scene=dict(
                xaxis=dict(title='x (AU)', backgroundcolor='white',
                        gridcolor='lightgrey', showbackground=True),
                yaxis=dict(title='y (AU)', backgroundcolor='white',
                        gridcolor='lightgrey', showbackground=True),
                zaxis=dict(title='z (AU)', backgroundcolor='white',
                        gridcolor='lightgrey', showbackground=True),
                bgcolor='white',
                aspectmode='data',
            ),
            paper_bgcolor='white',
            font=dict(color='black'),
            width=1000,
            height=800,
            legend=dict(
                x=0.01, y=0.99,
                bordercolor='lightgrey',
                borderwidth=1,
                font=dict(color='black', size=9),
            ),
        )

        if save_html:
            fig.write_html(save_html)
            print(f"Saved to {save_html}")
        if not _return_only:
            fig.show()
        return fig
    
    def _build_scatter_figure(self, df):
        """
        Build the MGA scatter figure from a candidates dataframe.
        Called by both plot() and launch_dashboard().
        """
        import pandas as pd_df

        launch_dates  = [et_to_datetime(t).strftime("%d %b %Y")
                        for t in df["t_launch_et"]]
        arrival_dates = [et_to_datetime(t).strftime("%d %b %Y")
                        for t in df["t_arrival_et"]]

        c3_vals  = df["c3_km2_s2"].values
        c3_norm  = (c3_vals - c3_vals.min()) / (c3_vals.max() - c3_vals.min() + 1e-10)
        sizes    = 6 + (1 - c3_norm) * 14

        hover = [
            f"<b>Index:</b>    {idx}<br>"
            f"<b>Launch:</b>   {r['t_launch_utc']}<br>"
            f"<b>Flyby:</b>    {r['t_flyby_utc']}<br>"
            f"<b>Arrival:</b>  {r['t_arrival_utc']}<br>"
            f"<b>TOF total:</b> {r['tof_total_days']:.0f} days<br>"
            f"<b>TOF leg 1:</b> {r['tof_leg1_days']:.0f} days<br>"
            f"<b>TOF leg 2:</b> {r['tof_leg2_days']:.0f} days<br>"
            f"<b>C3:</b>        {r['c3_km2_s2']:.1f} km²/s²<br>"
            f"<b>v∞ dep:</b>    {r['vinf_dep_km_s']:.2f} km/s<br>"
            f"<b>v∞ flyby:</b>  {r['vinf_flyby_km_s']:.2f} km/s<br>"
            f"<b>v∞ target:</b> {r['vinf_target_km_s']:.2f} km/s<br>"
            f"<b>Turn:</b>      {r['turn_angle_deg']:.1f}°<br>"
            f"<b>Flyby alt:</b> {r['flyby_alt_km']:.0f} km<br>"
            f"<b>SOI frac:</b>  {r['soi_fraction']:.3f}"
            for idx, (_, r) in enumerate(df.iterrows())
        ]

        n_cand     = len(df)
        c3_min_val = df["c3_km2_s2"].min()
        c3_max_val = df["c3_km2_s2"].max()
        v_min_val  = df["vinf_target_km_s"].min()
        v_max_val  = df["vinf_target_km_s"].max()
        seq_str    = " → ".join(self.sequence)

        fig = go.Figure()

        fig.add_trace(go.Scatter(
            x=launch_dates,
            y=arrival_dates,
            mode='markers',
            marker=dict(
                color=df["vinf_target_km_s"],
                colorscale='Viridis_r',
                size=sizes,
                colorbar=dict(
                    title=dict(
                        text="v∞ target (km/s)",
                        font=dict(color='black'),
                    ),
                    tickfont=dict(color='black'),
                ),
                showscale=True,
                line=dict(width=0.5, color='white'),
            ),
            customdata=list(range(len(df))),   # candidate index for click callback
            text=hover,
            hoverinfo='text',
            name='candidates',
        ))

        fig.update_layout(
            title=dict(
                text=(
                    f"MGA Scan — {seq_str}<br>"
                    f"<sup>{n_cand} candidates   "
                    f"C3: {c3_min_val:.1f}–{c3_max_val:.1f} km²/s²   "
                    f"v∞ target: {v_min_val:.2f}–{v_max_val:.2f} km/s   "
                    f"color: v∞ target   size: 1/C3</sup>"
                ),
                font=dict(size=14, color='black'),
                x=0.5,
                xanchor='center',
            ),
            xaxis=dict(
                title=dict(text="Launch Date", font=dict(color='black')),
                tickangle=30,
                tickfont=dict(color='black', size=9),
                gridcolor='lightgrey',
                showgrid=True,
                type='category',
                categoryorder='array',
                categoryarray=sorted(set(launch_dates),
                                    key=lambda d: datetime.strptime(d, "%d %b %Y")),
            ),
            yaxis=dict(
                title=dict(text="Target Arrival Date", font=dict(color='black')),
                tickfont=dict(color='black', size=9),
                gridcolor='lightgrey',
                showgrid=True,
                type='category',
                categoryorder='array',
                categoryarray=sorted(set(arrival_dates),
                                    key=lambda d: datetime.strptime(d, "%d %b %Y")),
            ),
            paper_bgcolor='white',
            plot_bgcolor='white',
            font=dict(color='black'),
            width=1000,
            height=700,
            legend=dict(
                x=1.08, y=1,
                bordercolor='lightgrey',
                borderwidth=1,
                font=dict(color='black'),
            ),
        )

        return fig


