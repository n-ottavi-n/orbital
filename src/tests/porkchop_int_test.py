from porkchop_interface_plotly import porkchop_interface
from tools import load_solar_system_kernels

load_solar_system_kernels()

origin_body = "EARTH" # must be SPICE bodies
target_body = "VENUS"

launch_date = "2028 JAN 09 00:00:00"
arrival_date = "2028 APR 27 00:00:00"
'''
max_c3 = 400 # c3 values above this will be masked

int = porkchop_interface(origin_body, target_body, launch_date, arrival_date, max_c3,span=180, grid_res=7)
int.plot()'''

pc = porkchop_interface(origin_body, target_body,
                        launch_date, arrival_date,
                        max_c3=100, span=360, grid_res=7)
pc.plot(save_html="reports/porkchops/earth_venus_2029.html")

