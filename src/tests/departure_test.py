
from optimize_periapsis import optimize_periapsis
from lambert_interface import lambert_interface
import planetary_data
from planetary_data import get_body, get_bodies
from departure_orbit_elements import departure_orbit_elements
from arrival_orbit_elements import arrival_orbit_elements
from plot_approach import plot_approach
import matplotlib.pyplot as plt
from tools import load_solar_system_kernels
import spiceypy as spice


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


mu_body_dest = get_body(dest)['mu']
mu_body_origin = get_body(origin)['mu']
radius = get_body(dest)['radius']

sim = lambert_interface(
    origin, dest,
    start_date, arrival_date, end_date,
    1000,          # minimal steps, just enough to propagate
    perturbations=[]   # no perturbations — faster
)
sim.solve()
sim.propagate()

departure = departure_orbit_elements(sim, origin, mu_body_origin)
print(departure)
