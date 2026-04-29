from optimize_periapsis import optimize_periapsis
from lambert_interface import lambert_interface
import planetary_data
from planetary_data import get_body, get_bodies
from arrival_orbit_elements import arrival_orbit_elements
from plot_approach import plot_approach
import matplotlib.pyplot as plt
from tools import load_solar_system_kernels


load_solar_system_kernels()

start_date='JAN 30, 2028, 23:00 UTC'
arrival_date='DEC 09, 2028, 00:00 UTC'
end_date='DEC 30, 2028, 00:00 UTC'

origin='EARTH'
dest='CERES'
perts=['EARTH', 'MOON', 'MARS']

perturbations=get_bodies(perts)

des_inc = None # desired inclination at arrival in deg
des_pe = 300 # desired periapsis at arrival in km


mu_body = get_body(dest)['mu']
radius = get_body(dest)['radius']

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