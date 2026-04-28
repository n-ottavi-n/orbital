from optimize_periapsis import optimize_periapsis
import spiceypy as spice
from lambert_interface import lambert_interface
import planetary_data
from arrival_orbit_elements import arrival_orbit_elements
from plot_approach import plot_approach
import matplotlib.pyplot as plt

spice.kclear()

import os
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
SPICE_DIR = os.path.normpath(os.path.join(BASE_DIR, '..', 'spice_solar_system'))
os.chdir(os.path.join(BASE_DIR, '..'))  # set working dir to src/
spice.furnsh(os.path.join(SPICE_DIR, 'solar_system_kernel.txt'))

start_date='JAN 30, 2028, 23:00 UTC'
arrival_date='DEC 09, 2028, 00:00 UTC'
end_date='DEC 30, 2028, 00:00 UTC'

origin='EARTH'
dest='CERES'

des_inc = 25 # desired inclination at arrival in deg
des_pe = 300 # desired periapsis at arrival in km

mu_body=planetary_data.ceres['mu']

perturbations=[planetary_data.earth, planetary_data.moon,]


# phase 1: fast coarse optimization
result_coarse = optimize_periapsis(
    lambert_interface,
    origin, dest,
    start_date, arrival_date, end_date,
    1000,
    perturbations,
    periapsis_altitude_km=des_pe,
    mu_body=mu_body,
    inc_target_deg=des_inc,
    arrival_window_days=30,
    max_iter=200
)

# phase 2: refine at full resolution from coarse solution
result = optimize_periapsis(
    lambert_interface,
    origin, dest,
    start_date, arrival_date, end_date,
    5000,
    perturbations,
    periapsis_altitude_km=des_pe,
    mu_body=mu_body,
    inc_target_deg=des_inc,
    arrival_window_days=30,
    max_iter=150,
    x0=result_coarse["u_rtn"],
    dt0=result_coarse["dt_arrival_days"]
)

print(result)


solved_sim = result["sim_object"]
arrival = arrival_orbit_elements(solved_sim, dest, mu_body)
print(arrival)

fig, ax = plot_approach(arrival, body_name=dest, body_radius_km=470)
plt.show()