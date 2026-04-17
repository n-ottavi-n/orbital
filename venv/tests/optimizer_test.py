from optimizer import optimize_c3
from lambert_interface import lambert_interface
import spiceypy as spice

spice.kclear()
spice.furnsh('../spice_solar_system/solar_system_kernel.txt')

best = optimize_c3(
    lambert_interface,
    origin='EARTH',
    dest='JUPITER BARYCENTER',
    base_start='Dec 25, 2028, 12:00 UTC',
    base_arrival='Jul 01, 2031, 00:00 UTC',
    days_range=1,
    step_hours=1
)

print(best)