from mga_scan import MGAscan

from tools import load_solar_system_kernels

load_solar_system_kernels()



pc = MGAscan(
    sequence          = ["EARTH", "JUPITER", "SATURN"],
    launch_window     = ("JAN 01 1977", "DEC 30 1977"),
    flyby_windows     = [("JAN 01 1979", "DEC 30 1979")],
    propagation_days  = 3000,
    vinf_bounds       = {
        "EARTH":   (3, 13),
        "JUPITER": (3, 13),
        "SATURN": (0, 100)
    },
    periapsis_radii   = [1.05, 1.1, 1.2, 1.3, 1.4, 1.5, 2, 3,4, 5, 6, 7, 8, 9, 10, 12, 15, 20, 30],
    min_flyby_alt_km  = {"JUPITER": 5000},
    grid_res_days     =7,
)

pc.compute(output_file="reports/ejs_scan.csv")
pc.plot(save_html="reports/ejs_scan.html")