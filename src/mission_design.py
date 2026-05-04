import numpy as np
import spiceypy as spice
import matplotlib.pyplot as plt
from lambert_interface import lambert_interface
from optimize_periapsis import optimize_periapsis
from arrival_orbit_elements import arrival_orbit_elements
from departure_orbit_elements import departure_orbit_elements
from plot_approach import plot_approach
from plot_departure import plot_departure
import planetary_data as pd
from tools import load_solar_system_kernels


class MissionDesign:

    def __init__(self,
                 origin,
                 destination,
                 launch_date,
                 arrival_date,
                 periapsis_alt_km=300.0,
                 inc_target_deg=None,
                 parking_alt_km=400.0,
                 arrival_window_days=30,
                 steps_coarse=1000,
                 steps_fine=5000,
                 max_iter_coarse=200,
                 max_iter_fine=150,
                 verbose=True):

        self.origin            = origin
        self.destination       = destination
        self.launch_date       = launch_date
        self.arrival_date      = arrival_date
        self.periapsis_alt_km  = periapsis_alt_km
        self.inc_target_deg    = inc_target_deg
        self.parking_alt_km    = parking_alt_km
        self.arrival_window_days = arrival_window_days
        self.steps_coarse      = steps_coarse
        self.steps_fine        = steps_fine
        self.max_iter_coarse   = max_iter_coarse
        self.max_iter_fine     = max_iter_fine
        self.verbose           = verbose

        # results — populated by solve()
        self.result_coarse     = None
        self.result            = None
        self.sim               = None
        self.departure         = None
        self.arrival           = None
        self.summary           = None

        # body data
        self.mu_origin         = pd.get_body(origin)['mu']
        self.mu_dest           = pd.get_body(destination)['mu']
        self.radius_origin     = pd.get_body(origin)['radius']
        self.radius_dest       = pd.get_body(destination)['radius']

        load_solar_system_kernels()

        # end date = arrival + window for propagation
        arrival_et  = spice.str2et(arrival_date)
        end_et      = arrival_et + (arrival_window_days + 10) * 86400
        self.end_date = spice.et2utc(end_et, "C", 0)

    # --------------------------------------------------
    def _log(self, msg):
        if self.verbose:
            print(f"[MissionDesign] {msg}")

    # --------------------------------------------------
    def solve(self):

        self._log(f"Starting mission design: {self.origin} → {self.destination}")
        self._log(f"Launch: {self.launch_date}   Arrival: {self.arrival_date}")
        self._log(f"Target periapsis: {self.periapsis_alt_km} km" +
                  (f"   Inclination: {self.inc_target_deg}°" if self.inc_target_deg else ""))

        # --------------------------------------------------
        # phase 1: coarse optimization
        # --------------------------------------------------
        self._log("Phase 1: coarse optimization...")

        self.result_coarse = optimize_periapsis(
            lambert_interface,
            self.origin, self.destination,
            self.launch_date, self.arrival_date, self.end_date,
            self.steps_coarse,
            perturbations=[],
            periapsis_altitude_km=self.periapsis_alt_km,
            mu_body=self.mu_dest,
            inc_target_deg=self.inc_target_deg,
            arrival_window_days=self.arrival_window_days,
            max_iter=self.max_iter_coarse
        )

        self._log(f"Coarse solution: "
                  f"periapsis alt = {self.result_coarse['periapsis_altitude_km']:.1f} km   "
                  f"dv = {self.result_coarse['dv_magnitude_ms']:.1f} m/s")

        # --------------------------------------------------
        # phase 2: fine optimization
        # --------------------------------------------------
        self._log("Phase 2: fine optimization...")

        self.result = optimize_periapsis(
            lambert_interface,
            self.origin, self.destination,
            self.launch_date, self.arrival_date, self.end_date,
            self.steps_fine,
            perturbations=[],
            periapsis_altitude_km=self.periapsis_alt_km,
            mu_body=self.mu_dest,
            inc_target_deg=self.inc_target_deg,
            arrival_window_days=self.arrival_window_days,
            max_iter=self.max_iter_fine,
            x0=self.result_coarse["u_rtn"],
            dt0=self.result_coarse["dt_arrival_days"]
        )

        self._log(f"Fine solution: "
                  f"periapsis alt = {self.result['periapsis_altitude_km']:.1f} km   "
                  f"dv = {self.result['dv_magnitude_ms']:.1f} m/s   "
                  f"{'✓' if self.result['success'] else '✗ iteration limit'}")

        # --------------------------------------------------
        # phase 3: compute elements
        # --------------------------------------------------
        self._log("Computing departure and arrival elements...")

        self.sim = self.result["sim_object"]

        self.departure = departure_orbit_elements(
            self.sim, self.origin, self.mu_origin,
            parking_alt_km=self.parking_alt_km
        )

        self.arrival = arrival_orbit_elements(
            self.sim, self.destination, self.mu_dest
        )

        # --------------------------------------------------
        # phase 4: assemble summary
        # --------------------------------------------------
        self._log("Assembling mission summary...")
        self._assemble_summary()
        self._log("Done.")

    # --------------------------------------------------
    def _assemble_summary(self):

        tof_total = (self.arrival["epoch_et"] - self.departure["injection_et"]) / 86400
        tof_helio = (self.arrival["epoch_et"] - self.departure["departure_et"]) / 86400

        self.summary = {
            # departure
            "origin":                  self.origin,
            "destination":             self.destination,
            "injection_utc":           self.departure["injection_utc"],
            "injection_dv_km_s":       self.departure["delta_v_km_s"],
            "parking_alt_km":          self.departure["parking_alt_km"],
            "parking_inc_deg":         self.departure["parking_inc_deg"],
            "c3_km2_s2":               self.departure["c3_km2_s2"],
            "vinf_dep_km_s":           self.departure["vinf_km_s"],
            "soi_exit_utc":            self.departure["departure_utc"],
            "t_soi_hours":             self.departure["t_soi_hours"],
            "asymptote_dec_deg":       self.departure["asymptote_dec_deg"],

            # transfer
            "tof_total_days":          tof_total,       # injection → arrival periapsis
            "tof_heliocentric_days":   tof_helio,       # SOI exit → arrival periapsis

            # arrival
            "arrival_utc":             spice.et2utc(self.arrival["epoch_et"], "C", 0),
            "arrival_shift_days":      self.result["dt_arrival_days"],
            "vinf_arr_km_s":           self.arrival["vinf_km_s"],
            "capture_dv_km_s":         (np.sqrt(self.arrival["vinf_km_s"]**2 +
                                        2 * self.mu_dest / self.arrival["periapsis_km"]) -
                                        np.sqrt(2 * self.mu_dest / self.arrival["periapsis_km"])),
            "periapsis_alt_km":        self.arrival["periapsis_altitude_km"],
            "inclination_deg":         self.arrival["inclination_deg"],
            "periapsis_lat_deg":       self.arrival["periapsis_latitude_deg"],
            "periapsis_lon_deg":       self.arrival["periapsis_longitude_deg"],
            "eccentricity":            self.arrival["eccentricity"],
            "raan_deg":                self.arrival["raan_deg"],

            # budget
            "total_dv_km_s":           (self.departure["delta_v_km_s"] +
                                        (np.sqrt(self.arrival["vinf_km_s"]**2 +
                                        2 * self.mu_dest / self.arrival["periapsis_km"]) -
                                        np.sqrt(2 * self.mu_dest / self.arrival["periapsis_km"]))),
        }

    # --------------------------------------------------
    def report(self):
        if self.summary is None:
            print("No solution yet — call solve() first.")
            return

        s = self.summary
        print("\n" + "="*60)
        print(f"  MISSION DESIGN REPORT")
        print(f"  {s['origin']} → {s['destination']}")
        print("="*60)

        print(f"\n  DEPARTURE")
        print(f"    Injection burn:    {s['injection_utc']}")
        print(f"    Parking orbit:     {s['parking_alt_km']:.0f} km  "
              f"i = {s['parking_inc_deg']:.2f}°")
        print(f"    Departure ΔV:      {s['injection_dv_km_s']:.3f} km/s")
        print(f"    C3:                {s['c3_km2_s2']:.3f} km²/s²")
        print(f"    v∞ departure:      {s['vinf_dep_km_s']:.3f} km/s")
        print(f"    Asymptote dec:     {s['asymptote_dec_deg']:.2f}°")
        print(f"    SOI exit:          {s['soi_exit_utc']}")
        print(f"    SOI transit:       {s['t_soi_hours']:.1f} h")

        print(f"\n  TRANSFER")
        print(f"    Total TOF:         {s['tof_total_days']:.1f} days  "
              f"(injection → arrival)")
        print(f"    Heliocentric TOF:  {s['tof_heliocentric_days']:.1f} days  "
              f"(SOI exit → arrival)")
        if s['arrival_shift_days'] != 0:
            print(f"    Arrival shift:     {s['arrival_shift_days']:+.2f} days")

        print(f"\n  ARRIVAL")
        print(f"    Arrival epoch:     {s['arrival_utc']}")
        print(f"    v∞ arrival:        {s['vinf_arr_km_s']:.3f} km/s")
        print(f"    Capture ΔV:        {s['capture_dv_km_s']:.3f} km/s  (parabolic)")
        print(f"    Periapsis alt:     {s['periapsis_alt_km']:.1f} km")
        print(f"    Inclination:       {s['inclination_deg']:.2f}°" +
              (f"  (target {self.inc_target_deg}°)" if self.inc_target_deg else ""))
        print(f"    Periapsis lat/lon: {s['periapsis_lat_deg']:.2f}°  "
              f"{s['periapsis_lon_deg']:.2f}°")
        print(f"    Eccentricity:      {s['eccentricity']:.4f}")
        print(f"    RAAN:              {s['raan_deg']:.2f}°")

        print(f"\n  MISSION ΔV BUDGET")
        print(f"    Departure burn:    {s['injection_dv_km_s']:.3f} km/s")
        print(f"    Capture burn:      {s['capture_dv_km_s']:.3f} km/s")
        print(f"    ─────────────────────────────")
        print(f"    Total:             {s['total_dv_km_s']:.3f} km/s")
        print("="*60 + "\n")

    # --------------------------------------------------
    def plot(self, save_report=False, report_path="mission_report.png"):

        if self.departure is None or self.arrival is None:
            print("No solution yet — call solve() first.")
            return

        fig1, ax1 = plot_departure(
            self.departure,
            body_name=self.origin.capitalize(),
            body_radius_km=self.radius_origin,
            mu_body=self.mu_origin
        )

        fig2, ax2 = plot_approach(
            self.arrival,
            body_name=self.destination.capitalize(),
            body_radius_km=self.radius_dest,
            mu_body=self.mu_dest
        )

        if save_report:
            self._save_report(fig1, fig2, report_path)
            summary_path = report_path.replace(".png", "_summary.txt")
            self._save_summary(summary_path)

        plt.show()

    # --------------------------------------------------
    def _save_report(self, fig1, fig2, path):
        """
        Composite departure and arrival plots into a single image.
        """
        import matplotlib.backends.backend_agg as agg
        import PIL.Image as Image
        import io

        def fig_to_image(fig):
            buf = io.BytesIO()
            fig.savefig(buf, format='png', dpi=150, bbox_inches='tight')
            buf.seek(0)
            return Image.open(buf).copy()

        img1 = fig_to_image(fig1)
        img2 = fig_to_image(fig2)

        # side by side
        total_width  = img1.width + img2.width
        total_height = max(img1.height, img2.height)

        combined = Image.new("RGB", (total_width, total_height), "white")
        combined.paste(img1, (0, 0))
        combined.paste(img2, (img1.width, 0))
        combined.save(path)

        self._log(f"Report saved to {path}")

    # -----------------------------------------------------
    def _save_summary(self, path):
        """
        Save mission summary to a formatted text file.
        """
        if self.summary is None:
            return

        with open(path, 'w', encoding='utf-8') as f:
            s = self.summary
            f.write("=" * 60 + "\n")
            f.write(f"  MISSION DESIGN REPORT\n")
            f.write(f"  {s['origin']} → {s['destination']}\n")
            f.write("=" * 60 + "\n")

            f.write(f"\n  DEPARTURE\n")
            f.write(f"    Injection burn:    {s['injection_utc']}\n")
            f.write(f"    Parking orbit:     {s['parking_alt_km']:.0f} km  "
                    f"i = {s['parking_inc_deg']:.2f}°\n")
            f.write(f"    Departure ΔV:      {s['injection_dv_km_s']:.3f} km/s\n")
            f.write(f"    C3:                {s['c3_km2_s2']:.3f} km²/s²\n")
            f.write(f"    v∞ departure:      {s['vinf_dep_km_s']:.3f} km/s\n")
            f.write(f"    Asymptote dec:     {s['asymptote_dec_deg']:.2f}°\n")
            f.write(f"    SOI exit:          {s['soi_exit_utc']}\n")
            f.write(f"    SOI transit:       {s['t_soi_hours']:.1f} h\n")

            f.write(f"\n  TRANSFER\n")
            f.write(f"    Total TOF:         {s['tof_total_days']:.1f} days\n")
            f.write(f"    Heliocentric TOF:  {s['tof_heliocentric_days']:.1f} days\n")
            if s['arrival_shift_days'] != 0:
                f.write(f"    Arrival shift:     {s['arrival_shift_days']:+.2f} days\n")

            f.write(f"\n  ARRIVAL\n")
            f.write(f"    Arrival epoch:     {s['arrival_utc']}\n")
            f.write(f"    v∞ arrival:        {s['vinf_arr_km_s']:.3f} km/s\n")
            f.write(f"    Capture ΔV:        {s['capture_dv_km_s']:.3f} km/s  (parabolic)\n")
            f.write(f"    Periapsis alt:     {s['periapsis_alt_km']:.1f} km\n")
            f.write(f"    Inclination:       {s['inclination_deg']:.2f}°" +
                    (f"  (target {self.inc_target_deg}°)\n" if self.inc_target_deg else "\n"))
            f.write(f"    Periapsis lat/lon: {s['periapsis_lat_deg']:.2f}°  "
                    f"{s['periapsis_lon_deg']:.2f}°\n")
            f.write(f"    Eccentricity:      {s['eccentricity']:.4f}\n")
            f.write(f"    RAAN:              {s['raan_deg']:.2f}°\n")

            f.write(f"\n  MISSION ΔV BUDGET\n")
            f.write(f"    Departure burn:    {s['injection_dv_km_s']:.3f} km/s\n")
            f.write(f"    Capture burn:      {s['capture_dv_km_s']:.3f} km/s\n")
            f.write(f"    {'─'*33}\n")
            f.write(f"    Total:             {s['total_dv_km_s']:.3f} km/s\n")
            f.write("=" * 60 + "\n")

        self._log(f"Summary saved to {path}")