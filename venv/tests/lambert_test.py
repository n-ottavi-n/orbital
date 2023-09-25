import math
import planetary_data
from lambert_tools import lambert_uv
import spiceypy as spice
import numpy as np
from plot_interplanetary_interface import plot_trajectory


spice.furnsh('../spice_solar_system/solar_system_kernel.txt')

start_date='Jul 30, 2020, 00:00 UTC'
end_date='Feb 18, 2021, 00:00 UTC'

steps=10000

etOne = spice.str2et(start_date)
etTwo = spice.str2et(end_date)
tof=etTwo-etOne
dt = (etTwo - etOne) / steps


# initilize times array
times = [x * (etTwo - etOne) / steps + etOne for x in range(steps)]  # [start_date : end_date] len=steps

obs = 'SOLAR SYSTEM BARYCENTER'
props=[]

perturbations=[]
bodies=['earth','mars']

states_bodies=[]
#get states of other bodies to plot e.g. perturbating bodies
for body in bodies:
    positions, lightTimes = spice.spkpos(body, times, 'ECLIPJ2000', 'NONE', obs)
    # append states to propagator states
    props.append(positions)


start_r=props[0][0]
end_r=props[1][-1]

v=lambert_uv(start_r,end_r,tof, 0.8, 4*math.pi**2,-4*math.pi**2)

v0=v[:3] #departure velocity vector
x0=start_r#departure position vector
xv=np.append(x0,v0)
xv=xv.tolist()

states=[['spacecraft'],[xv]]
print("deltaV departure: {} km/s".format(np.linalg.norm(v0)-29.78))

plot_trajectory(states,bodies,start_date,end_date,perturbations,steps=steps,show=True, animate=False)