import math
import planetary_data
from lambert_tools import lambert_uv
import spiceypy as spice
import numpy as np
from plot_interplanetary_interface import interplanetary_interface


spice.furnsh('../spice_solar_system/solar_system_kernel.txt')

start_date='Jul 30, 2020, 00:00 UTC'
arrival_date='Feb 18, 2021, 00:00 UTC'
end_date='Jul 30, 2022, 00:00 UTC'

steps=10000

etOne = spice.str2et(start_date)
etTwo = spice.str2et(end_date)
etArr = spice.str2et(arrival_date)
tof=etArr-etOne
dt = (etTwo - etOne) / steps


# initilize times array
times = [x * (etTwo - etOne) / steps + etOne for x in range(steps)]  # [start_date : end_date] len=steps

obs = 'SOLAR SYSTEM BARYCENTER'
props=[]

perturbations=[]
origin='earth'
dest='mars'
bodies=[origin,dest]

states_bodies=[]
#get states of other bodies to plot e.g. perturbating bodies
for body in bodies:
    positions, lightTimes = spice.spkpos(body, times, 'ECLIPJ2000', 'NONE', obs)

    # append states to propagator states
    props.append(positions)

position_tgt, lightTimes = spice.spkpos(dest, etArr, 'ECLIPJ2000', 'NONE', obs)


start_r=props[0][0] #position vector of origin at start date
#end_r=props[1][-1]
end_r=position_tgt # position vector of target at arrival date


v=lambert_uv(start_r,end_r,tof, 0.8, 4*math.pi**2,-4*math.pi**2)

v0=v[:3] #departure velocity vector
x0=start_r#departure position vector
xv=np.append(x0,v0)
xv=xv.tolist()

states=[['spacecraft'],[xv]]
origin_v=(props[0][1]-props[0][0])/dt #obital velocity vector at origin
deltaV=v0-origin_v #deltaV vector at origin
print("deltaV departure: {} km/s\ntotal: {} km/s".format(deltaV, np.linalg.norm(deltaV)))

spice.kclear()

int=interplanetary_interface(states,bodies,start_date,end_date,perturbations,steps=steps)
int.plot_trajectory(show=True, animate=False)