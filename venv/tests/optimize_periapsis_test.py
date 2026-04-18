from optimize_periapsis import optimize_periapsis
import spiceypy as spice
from lambert_interface import lambert_interface
import planetary_data

spice.kclear()
spice.furnsh('../spice_solar_system/solar_system_kernel.txt')

start_date='Oct 26, 2029, 04:00 UTC'
arrival_date='Apr 04, 2030, 00:00 UTC'
end_date='Apr 12, 2030, 00:00 UTC'

origin='EARTH'
dest='VENUS'

perturbations=[planetary_data.earth, planetary_data.moon, planetary_data.venus,]


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