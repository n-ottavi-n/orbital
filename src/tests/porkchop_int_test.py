from porkchop_interface import porkchop_interface
import spiceypy as spice

spice.kclear()

import os
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
SPICE_DIR = os.path.normpath(os.path.join(BASE_DIR, '..', 'spice_solar_system'))
os.chdir(os.path.join(BASE_DIR, '..'))  # set working dir to src/
spice.furnsh(os.path.join(SPICE_DIR, 'solar_system_kernel.txt'))


origin_body = "EARTH" # must be SPICE bodies
target_body = "CERES"

launch_date = "2027 DEC 30 00:00:00"
arrival_date = "2028 OCT 12 00:00:00"

max_tof = 600 # in days solutions longer than this will not be tested
max_c3 = 200 # c3 values above this will be masked

int = porkchop_interface(origin_body, target_body, launch_date, arrival_date, max_tof, max_c3,span=180, grid_res=7)
int.plot()