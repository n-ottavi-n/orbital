from porkchop_interface import porkchop_interface
import spiceypy as spice

spice.kclear()
spice.furnsh('../spice_solar_system/solar_system_kernel.txt')

origin_body = "EARTH" # must be SPICE bodies
target_body = "MARS"

launch_date = "2028 OCT 31 00:00:00"
arrival_date = "2029 MAR 12 00:00:00"

max_tof = 300 # in days solutions longer than this will not be tested
max_c3 = 50 # c3 values above this will be masked

int = porkchop_interface(origin_body, target_body, launch_date, arrival_date, max_tof, max_c3,span=300, grid_res=7)
int.plot()