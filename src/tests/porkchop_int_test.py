from porkchop_interface_plotly import porkchop_interface
from tools import load_solar_system_kernels

load_solar_system_kernels()

origin_body = "VENUS" # must be SPICE bodies
target_body = "MERCURY"

launch_date = "2031 NOV 21 00:00:00"
arrival_date = "2031 DEC 15 00:00:00"
'''
max_c3 = 400 # c3 values above this will be masked

int = porkchop_interface(origin_body, target_body, launch_date, arrival_date, max_c3,span=180, grid_res=7)
int.plot()'''

pc = porkchop_interface(origin_body, target_body,
                        launch_date, arrival_date,
                        max_c3=38, span=90, grid_res=1)
pc.plot(save_html="reports/porkchops/venus_mercury_2031.html")