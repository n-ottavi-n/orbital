from porkchop_interface import porkchop_interface
from tools import load_solar_system_kernels

load_solar_system_kernels()

origin_body = "EARTH" # must be SPICE bodies
target_body = "JUPITER"

launch_date = "2030 JAN 17 00:00:00"
arrival_date = "2032 APR 15 00:00:00"

max_c3 = 400 # c3 values above this will be masked

int = porkchop_interface(origin_body, target_body, launch_date, arrival_date, max_c3,span=180, grid_res=7)
int.plot()