
from lambert_interface import lambert_interface
from planetary_data import get_body, get_bodies
from departure_orbit_elements import departure_orbit_elements
from plot_departure import plot_departure
import matplotlib.pyplot as plt
from tools import load_solar_system_kernels


load_solar_system_kernels()

start_date='DEC 01, 2028, 23:00 UTC'
arrival_date='JUL 01, 2029, 00:00 UTC'
end_date='JUL 30, 2029, 00:00 UTC'

origin='EARTH'
dest='MARS'
perts=['EARTH', 'MOON']

perturbations=get_bodies(perts)

des_inc = 25 # desired inclination at arrival in deg
des_pe = 300 # desired periapsis at arrival in km
parking_alt = 300 #parking orbit altitude in km


mu_body_dest = get_body(dest)['mu']
mu_body_origin = get_body(origin)['mu']
radius = get_body(origin)['radius']

sim = lambert_interface(
    origin, dest,
    start_date, arrival_date, end_date,
    500,          # minimal steps, just enough to propagate
    perturbations=[]   # no perturbations — faster
)
sim.solve()
sim.propagate()

departure = departure_orbit_elements(sim, origin, mu_body_origin, parking_alt_km=parking_alt)
print(departure)
fig, ax = plot_departure(departure, body_name=origin, body_radius_km=radius, mu_body=mu_body_origin, frame=None)
plt.show()
