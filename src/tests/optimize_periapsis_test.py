from optimize_periapsis import optimize_periapsis
from lambert_interface import lambert_interface
from planetary_data import get_body, get_bodies
from departure_orbit_elements import departure_orbit_elements
from arrival_orbit_elements import arrival_orbit_elements
from plot_departure import plot_departure
from plot_approach import plot_approach
from tools import load_solar_system_kernels


load_solar_system_kernels()

start_date='JAN 30, 2031, 23:00 UTC'
arrival_date='JUL 12, 2031, 00:00 UTC'
end_date='JUL 30, 2031, 00:00 UTC'

origin='EARTH'
dest='MARS'
perts=['EARTH', 'MOON']

perturbations=get_bodies(perts)

des_inc = 45 # desired inclination at arrival in deg
des_pe = 300 # desired periapsis at arrival in km


mu_body_dest = get_body(dest)['mu']
mu_body_origin = get_body(origin)['mu']
radius = get_body(dest)['radius']

# phase 1: fast coarse optimization
result_coarse = optimize_periapsis(
    lambert_interface,
    origin, dest,
    start_date, arrival_date, end_date,
    1000,
    perturbations,
    periapsis_altitude_km=des_pe,
    mu_body=mu_body_dest,
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
    mu_body=mu_body_dest,
    inc_target_deg=des_inc,
    arrival_window_days=30,
    max_iter=150,
    x0=result_coarse["u_rtn"],
    dt0=result_coarse["dt_arrival_days"]
)


solved_sim = result["sim_object"]

departure = departure_orbit_elements(solved_sim, origin, mu_body_origin)
arrival = arrival_orbit_elements(solved_sim, dest, mu_body_dest)

print(departure)
print(arrival)

total_flight_time_days = (arrival["epoch_et"] - departure["injection_et"]) / 86400
print("TOF: ",total_flight_time_days, " days")