from porkchop_interface import porkchop_interface
from tools import load_solar_system_kernels

load_solar_system_kernels()

origin_body = "EARTH" # must be SPICE bodies
target_body = "CERES"

launch_date = "2027 DEC 30 00:00:00"
arrival_date = "2028 OCT 12 00:00:00"

max_tof = 600 # in days solutions longer than this will not be tested
max_c3 = 200 # c3 values above this will be masked

int = porkchop_interface(origin_body, target_body, launch_date, arrival_date, max_tof, max_c3,span=180, grid_res=7)
int.plot()