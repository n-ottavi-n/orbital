import numpy as np
import spiceypy as spice
import matplotlib.pyplot as plt
from trajectory_utils import lambert_branches
from lambert_interface import lambert_interface
from datetime import datetime

spice.kclear()
spice.furnsh('../spice_solar_system/solar_system_kernel.txt')

def et_to_datetime(et):
    utc_str = spice.et2utc(et, "ISOC", 0)  # e.g. 2028-12-12T00:00:00
    return datetime.strptime(utc_str, "%Y-%m-%dT%H:%M:%S")

origin_body = "EARTH" # must be SPICE bodies
target_body = "VENUS"

base_launch = spice.str2et("2020 MAR 12 00:00:00")
base_arrival = spice.str2et("2020 SEP 12 00:00:00")

max_c3 = 20 # c3 values above this will be masked

span = 60 # will look at launch and arrival dates ±span days around the base_launch and base_arrival
launch_span = span * 86400
arrival_span = span * 86400

grid_res = 24 #grid resolution in hours
step = grid_res * 3600


launches = np.arange(base_launch - launch_span,
                      base_launch + launch_span,
                      step)

arrivals = np.arange(base_arrival - arrival_span,
                      base_arrival + arrival_span,
                      step)

C3 = np.zeros((len(launches), len(arrivals)))
C3[:] = np.nan

mu = 1.32712440018e11  # Sun GM

for i, t0 in enumerate(launches):
    for j, t1 in enumerate(arrivals):

        if t1 <= t0:
            continue

        try:
            start_str = spice.et2utc(t0, "C", 0)
            arr_str   = spice.et2utc(t1, "C", 0)

            sim = lambert_interface(
                origin_body, target_body,
                start_str,
                arr_str,
                arr_str,
                steps=500,
                perturbations=[]
            )

            sim.solve()

            state, _ = spice.spkezr("EARTH", t0,
                                    'ECLIPJ2000', 'NONE', 'SUN')
            v_planet = state[3:]

            sols = lambert_branches(sim.start_r,
                                    sim.end_r,
                                    sim.tof,
                                    mu=mu)

            best_c3 = np.inf

            for sol in sols:
                v0 = sol[:3]
                v_inf = v0 - v_planet
                c3 = np.dot(v_inf, v_inf)
                best_c3 = min(best_c3, c3)

            C3[i, j] = best_c3


        except Exception as e:

            print("FAIL:", i, j, e)

            continue

# mask high C3
C3_plot = np.copy(C3)
C3_plot[C3_plot > max_c3] = np.nan
C3_masked = np.ma.masked_invalid(C3_plot) # mask very large values for visibility

plt.figure(figsize=(10, 6))

# plot heatmap
cf = plt.contourf(launches/86400,
                  arrivals/86400,
                  C3_masked.T,
             levels=50)

plt.colorbar(cf, label="C3 (km^2/s^2)")

# plot C3 contours
levels = np.arange(8, 16, 0.5) #resolution of c3 values

cs = plt.contour(launches/86400,
                  arrivals/86400,
                  C3_masked.T,
                 levels=levels,
                 linewidths=0.8)

plt.clabel(cs, inline=True, fontsize=8, fmt="%.1f")

# plot TOF lines
TOF = np.zeros_like(C3)

for i, t0 in enumerate(launches):
    for j, t1 in enumerate(arrivals):
        TOF[i, j] = (t1 - t0) / 86400  # days

tof_levels = np.arange(100, 400, 15) #resolution of flight times

cs_tof = plt.contour(launches/86400,
                     arrivals/86400,
                     TOF.T,
                     levels=tof_levels,
                     linestyles='dashed',
                     colors='red')

plt.clabel(cs_tof, inline=True, fontsize=8, fmt="%d d", colors='red')

# define clean ticks for pretty graph
tick_idx = np.linspace(0, len(launches)-1, 6, dtype=int)

launch_ticks = launches[tick_idx] / 86400
launch_labels = [
    et_to_datetime(launches[i]).strftime("%d %b")
    for i in tick_idx
]

tick_idx_y = np.linspace(0, len(arrivals)-1, 6, dtype=int)

arrival_ticks = arrivals[tick_idx_y] / 86400
arrival_labels = [
    et_to_datetime(arrivals[i]).strftime("%d %b")
    for i in tick_idx_y
]

plt.xticks(launch_ticks, launch_labels, rotation=30)
plt.yticks(arrival_ticks, arrival_labels)

#legend and labels
start_year = et_to_datetime(launches[0]).year
end_year   = et_to_datetime(arrivals[-1]).year

plt.title(f"{origin_body} → {target_body} transfer\nLaunch: {start_year} | Arrival: {end_year}")

plt.xlabel("Launch Date")
plt.ylabel("Arrival Date")

plt.show()