from optimize_periapsis import optimize_periapsis
import spiceypy as spice
from lambert_interface import lambert_interface
import planetary_data
from arrival_orbit_elements import arrival_orbit_elements

spice.kclear()
spice.furnsh('../spice_solar_system/solar_system_kernel.txt')

start_date='DEC 09, 2028, 23:00 UTC'
arrival_date='JUL 15, 2029, 00:00 UTC'
end_date='AUG 12, 2029, 00:00 UTC'

origin='EARTH'
dest='MARS'

perturbations=[planetary_data.earth, planetary_data.moon,]


result = optimize_periapsis(
    lambert_interface,
    origin,
    dest,
    start_date,
    arrival_date,
    end_date,
    5000,
    perturbations,
    periapsis_altitude_km=300,
    max_iter=80
)

print(result)

solved_sim = result["sim_object"]
arrival = arrival_orbit_elements(solved_sim, "MARS", planetary_data.mars['mu'])
print(arrival)