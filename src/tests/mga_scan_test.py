from mga_scan import MGAscan

from tools import load_solar_system_kernels

import pandas as pd

load_solar_system_kernels()



scan = MGAscan(
    sequence          = ["VENUS", "EARTH", "JUPITER"],
    launch_window     = ("JAN 01 2026", "DEC 30 2030"),
    flyby_windows     = [("JAN 01 2026", "DEC 30 2033")],
    propagation_days  = 3000,
    vinf_bounds       = {
        "VENUS":   (3, 15),
        "EARTH": (3, 20),
        "JUPITER": (0, 100)
    },
    periapsis_radii   = [1.05, 1.1, 1.2, 1.3, 1.4, 1.5, 2, 3,4, 5, 6, 7, 8, 9, 10, 12, 15, 20, 30],
    min_flyby_alt_km  = {"EARTH": 200},
    grid_res_days     =7,
)

scan.compute(output_file="reports/vej_scan.csv")
scan.launch_dashboard()