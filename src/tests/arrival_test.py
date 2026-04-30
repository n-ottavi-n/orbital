from lambert_interface import lambert_interface
from planetary_data import get_body, get_bodies
from arrival_orbit_elements import arrival_orbit_elements
from plot_approach import plot_approach
import matplotlib.pyplot as plt
from tools import load_solar_system_kernels
from optimize_periapsis import optimize_periapsis


load_solar_system_kernels()

start_date='DEC 01, 2028, 23:00 UTC'
arrival_date='JUL 01, 2029, 00:00 UTC'
end_date='JUL 30, 2029, 00:00 UTC'

origin='EARTH'
dest='MARS'
perts=[]

perturbations=get_bodies(perts)

des_inc = 45 # desired inclination at arrival in deg
des_pe = 300 # desired periapsis at arrival in km


mu_body_dest = get_body(dest)['mu']
radius = get_body(dest)['radius']

# phase 1: fast coarse optimization
result = optimize_periapsis(
    lambert_interface,
    origin, dest,
    start_date, arrival_date, end_date,
    500,
    perturbations,
    periapsis_altitude_km=des_pe,
    mu_body=mu_body_dest,
    inc_target_deg=des_inc,
    arrival_window_days=30,
    max_iter=200
)


solved_sim = result["sim_object"]
arrival = arrival_orbit_elements(solved_sim, dest, mu_body_dest)
print(arrival)
fig, ax = plot_approach(arrival, body_name=dest, body_radius_km=radius)
plt.show()