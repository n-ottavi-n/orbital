from mga_scan import MGAscan

from tools import load_solar_system_kernels

import pandas as pd

load_solar_system_kernels()

departure="VENUS"
flyby="EARTH"
target="JUPITER"


scan = MGAscan(
    sequence          = [departure, flyby, target],
    launch_window     = ("JAN 01 2028", "DEC 30 2030"),
    flyby_windows     = [("FEB 01 2028", "DEC 28 2031")],
    propagation_days  = 3000,
    vinf_bounds       = {
        departure:   (0, 15),
        flyby: (3, 20),
        target: (0, 100)
    },
    periapsis_radii   = [1.05, 1.1, 1.2, 1.3, 1.4, 1.5, 2, 3,4, 5, 6, 7, 8, 9, 10, 12, 15, 20, 30],
    min_flyby_alt_km  = {flyby: 200},
    grid_res_days     =14,
)

scan.compute(output_file="reports/eve_scan.csv")
scan.launch_dashboard()