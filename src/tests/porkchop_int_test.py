from porkchop_interface_plotly import porkchop_interface
from tools import load_solar_system_kernels

load_solar_system_kernels()

origin_body = "EARTH" # must be SPICE bodies
target_body = "JUPITER"

launch_date = "2030 JAN 17 00:00:00"
arrival_date = "2032 APR 15 00:00:00"
'''
max_c3 = 400 # c3 values above this will be masked

int = porkchop_interface(origin_body, target_body, launch_date, arrival_date, max_c3,span=180, grid_res=7)
int.plot()'''

pc = porkchop_interface("EARTH", "MARS",
                        "DEC 09, 2028", "JUL 15, 2029",
                        max_c3=60, span=90, grid_res=5)
pc.plot(save_html="earth_mars_2028.html")