import planetary_data as pd
from lambert_tools import lambert_uv, lambert_universal
import spiceypy as spice
import numpy as np
from interplanetary_interface import interplanetary_interface


class lambert_interface:
    def __init__(self,origin, dest, start_date,arrival_date,end_date,steps, perturbations, mu=pd.sun['mu'], arrival_et_offset=0.0):

        self.mu = mu
        self.start_date=start_date
        self.end_date=end_date

        etOne = spice.str2et(start_date)
        etTwo = spice.str2et(end_date)

        # arrival_et_offset is in seconds — optimizer shifts this
        self.etArr = spice.str2et(arrival_date) + arrival_et_offset
        self.tof=self.etArr-etOne
        self.dt = (etTwo - etOne) / steps
        self.steps=steps

        self.start_et = etOne
        self.arrival_et = self.etArr  # consistent with offset
        self.end_et = etTwo




        # initilize times array
        self.times = [x * (etTwo - etOne) / steps + etOne for x in range(steps)]  # [start_date : end_date] len=steps

        #self.obs = 'SOLAR SYSTEM BARYCENTER'
        self.obs = 'SUN'

        self.perturbations=perturbations
        self.origin=origin
        self.dest=dest
        self.bodies=[origin,dest]


    def solve(self, prograde=True):
        props = []
        #get states of other bodies to plot e.g. perturbating bodies
        for body in self.bodies:
            positions, lightTimes = spice.spkpos(body, self.times, 'ECLIPJ2000', 'NONE', self.obs)
            # append states to propagator states
            props.append(positions)
        self.props=props

        #position_tgt, lightTimes = spice.spkpos(self.dest, self.etArr, 'ECLIPJ2000', 'NONE', self.obs)

        self.start_r, _ = spice.spkpos(self.origin, self.start_et, 'ECLIPJ2000', 'NONE', self.obs)
        self.end_r, _ = spice.spkpos(self.dest, self.arrival_et, 'ECLIPJ2000', 'NONE', self.obs)

        #self.v=lambert_uv(self.start_r,self.end_r,self.tof, 0.8, 4*math.pi**2,-4*math.pi**2)
        self.v=lambert_universal(self.start_r,self.end_r,self.tof, self.mu, prograde=prograde)
        self.v0 = self.v[:3]

    def propagate(self,show=False, animate=False):

        #self.v0=self.v[:3] #departure velocity vector
        x0=self.start_r#departure position vector
        xv=np.append(x0,self.v0)
        xv=xv.tolist()
        states=[['spacecraft'],[xv]]
        int=interplanetary_interface(states,self.bodies,self.start_date,self.end_date,self.perturbations,steps=self.steps)
        int.plot_trajectory(show=False, animate=animate)
        self.int_props=int.props
        self.int_vels = int.vels
        self.int_states = int.states

    def plot(self,show=False, animate=False):

        self.v0=self.v[:3] #departure velocity vector
        x0=self.start_r#departure position vector
        xv=np.append(x0,self.v0)
        xv=xv.tolist()
        states=[['spacecraft'],[xv]]
        int=interplanetary_interface(states,self.bodies,self.start_date,self.end_date,self.perturbations,steps=self.steps)
        int.plot_trajectory(show=show, animate=animate)
        self.int_props=int.props


    def deltaV(self):

        state, _ = spice.spkezr(self.origin, self.times[0], 'ECLIPJ2000', 'NONE', self.obs)
        print("states spkezr: {} ".format(state))
        origin_v = state[3:]  # obital velocity vector at origin
        #origin_v = (self.props[0][1] - self.props[0][0]) / self.dt  # obital velocity vector at origin
        v_inf_vectors = self.v0 - origin_v
        v_inf_norm = np.linalg.norm(v_inf_vectors)  # v infinity km/s
        print("v_inf departure: {} km/s".format(v_inf_norm))
        c3 = v_inf_norm**2 # injection C3 km^2/s^2
        print("injection C3: {} km²/s²".format(c3))
        '''
        deltaV = self.v0 - origin_v  # deltaV vector at origin
        print("deltaV departure: {} km/s\ntotal: {} km/s".format(deltaV, np.linalg.norm(deltaV)))
'''
        positions_sc=self.int_props[0] # positions of spacecraft
        positions_tgt=self.int_props[-1] # positions of target body
        positions_rel=positions_sc-positions_tgt #vectors from sc to tgt
        ranges=np.linalg.norm(positions_rel, axis=1) #distances from sc to tgt
        closest_dist=min(ranges) #closest approach
        print("closest approach: {} km".format(closest_dist))

        return closest_dist