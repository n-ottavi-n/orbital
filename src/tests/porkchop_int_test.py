from porkchop_interface import porkchop_interface
from tools import load_solar_system_kernels

load_solar_system_kernels()

origin_body = "EARTH" # must be SPICE bodies
target_body = "MARS"

launch_date = "2028 DEC 30 00:00:00"
arrival_date = "2029 JUN 12 00:00:00"

max_tof = 350 # in days solutions longer than this will not be tested
max_c3 = 50 # c3 values above this will be masked

int = porkchop_interface(origin_body, target_body, launch_date, arrival_date, max_tof, max_c3,span=60, grid_res=1)
int.plot()