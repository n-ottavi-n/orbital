from optimize_periapsis import optimize_periapsis
import spiceypy as spice
from lambert_interface import lambert_interface
import planetary_data
from arrival_orbit_elements import arrival_orbit_elements
from plot_approach import plot_approach
import matplotlib.pyplot as plt

spice.kclear()
spice.furnsh('../spice_solar_system/solar_system_kernel.txt')

start_date='DEC 09, 2028, 23:00 UTC'
arrival_date='JUL 15, 2029, 00:00 UTC'
end_date='AUG 12, 2029, 00:00 UTC'

origin='EARTH'
dest='MARS'

mu_body=planetary_data.mars['mu']

perturbations=[planetary_data.earth, planetary_data.moon,]


# phase 1: fast coarse optimization
result_coarse = optimize_periapsis(
    lambert_interface,
    origin, dest,
    start_date, arrival_date, end_date,
    1000,
    perturbations,
    periapsis_altitude_km=300,
    mu_body=mu_body,
    inc_target_deg=42,
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
    periapsis_altitude_km=300,
    mu_body=mu_body,
    inc_target_deg=42,
    arrival_window_days=30,
    max_iter=150,
    x0=result_coarse["u_rtn"],
    dt0=result_coarse["dt_arrival_days"]
)

print(result)


solved_sim = result["sim_object"]
arrival = arrival_orbit_elements(solved_sim, "MARS", mu_body)
print(arrival)

fig, ax = plot_approach(arrival, body_name="Mars", body_radius_km=3389.5)
plt.show()