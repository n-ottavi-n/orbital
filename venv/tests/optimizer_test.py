from optimizer import optimize_c3
from lambert_interface import lambert_interface
import spiceypy as spice

spice.kclear()
spice.furnsh('../spice_solar_system/solar_system_kernel.txt')

best = optimize_c3(
    lambert_interface,
    origin='EARTH',
    dest='VENUS',
    base_start='Mar 28, 2020, 12:00 UTC',
    base_arrival='Sep 20, 2020, 00:00 UTC',
    days_range=1,
    step_hours=1
)

print(best)