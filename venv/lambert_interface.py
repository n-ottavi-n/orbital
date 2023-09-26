import math
import planetary_data
from lambert_tools import lambert_uv
import spiceypy as spice
import numpy as np
from plot_interplanetary_interface import interplanetary_interface

start_date='Jul 30, 2020, 00:00 UTC'
arrival_date='Feb 18, 2021, 00:00 UTC'
end_date='Jul 30, 2022, 00:00 UTC'

steps=10000


class lambert_interface:
    def __init__(self,origin, dest, start_date,arrival_date,end_date,steps):

        spice.furnsh('../spice_solar_system/solar_system_kernel.txt')
        self.start_date=start_date
        self.end_date=end_date
        etOne = spice.str2et(start_date)
        etTwo = spice.str2et(end_date)
        self.etArr = spice.str2et(arrival_date)
        self.tof=self.etArr-etOne
        self.dt = (etTwo - etOne) / steps
        self.steps=steps

        # initilize times array
        self.times = [x * (etTwo - etOne) / steps + etOne for x in range(steps)]  # [start_date : end_date] len=steps

        self.obs = 'SOLAR SYSTEM BARYCENTER'


        self.perturbations=[]
        self.origin=origin
        self.dest=dest
        self.bodies=[origin,dest]

    def solve(self):
        props = []
        #get states of other bodies to plot e.g. perturbating bodies
        for body in self.bodies:
            positions, lightTimes = spice.spkpos(body, self.times, 'ECLIPJ2000', 'NONE', self.obs)
            # append states to propagator states
            props.append(positions)
        self.props=props

        position_tgt, lightTimes = spice.spkpos(self.dest, self.etArr, 'ECLIPJ2000', 'NONE', self.obs)


        self.start_r=props[0][0] #position vector of origin at start date
        #end_r=props[1][-1]
        self.end_r=position_tgt # position vector of target at arrival date


        self.v=lambert_uv(self.start_r,self.end_r,self.tof, 0.8, 4*math.pi**2,-4*math.pi**2)

    def plot(self,show=False, animate=False):

        self.v0=self.v[:3] #departure velocity vector
        x0=self.start_r#departure position vector
        xv=np.append(x0,self.v0)
        xv=xv.tolist()
        states=[['spacecraft'],[xv]]
        int=interplanetary_interface(states,self.bodies,self.start_date,self.end_date,self.perturbations,steps=self.steps)
        int.plot_trajectory(show=True, animate=False)

    def deltaV(self):
        origin_v = (self.props[0][1] - self.props[0][0]) / self.dt  # obital velocity vector at origin
        deltaV = self.v0 - origin_v  # deltaV vector at origin
        print("deltaV departure: {} km/s\ntotal: {} km/s".format(deltaV, np.linalg.norm(deltaV)))