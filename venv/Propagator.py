import numpy as np
import matplotlib.pyplot as plt
import planetary_data
from scipy.integrate import ode
import tools as t
import planetary_data as pd
import matplotlib.animation as animation
import spiceypy as spice

class Propagator:
    '''
    inputs:
        if cartesian: [x0,y0,z0,vx0,vy0,vz0]
        if keplerian: a,e,i,m,aop,raan
    '''
    def __init__(self, state0, n_steps, start_date='Sep 16, 2023, 00:00 UTC', end_date='Sep 17, 2023, 00:00 UTC',spacecraft_data={}, coes=False, deg=False, cb=pd.earth, perts=[], integrator='dopri5'):
        '''
        @param state0: if cartesian: [x0,y0,z0,vx0,vy0,vz0]  if keplerian: [a,e,i,m,aop,raan]
        @param n_steps: number of timesteps for integration
        @param start_date: format:  '%b %d, %Y, %H:%M %Z'
        @param end_date: format:  '%b %d, %Y, %H:%M %Z'
        @param coes: bool: True if state0 is keplerian elements
        @param deg: bool: True if input keplerian elements in degrees
        @param cb: central body
        @param perts: perturbations to account for / ['j2','srp', [pd.mars,pd.venus,pd.moon,.....]]
        '''
        if coes:
            self.r0,self.v0=t.coes2rv(state0, deg=deg, mu=cb['mu'])
        else:
            self.r0=state0[:3]
            self.v0=state0[3:]

        #self.tspan = tspan
        self.n_steps=n_steps
        self.cb = cb
        self.perts=perts
        self.integrator=integrator

        #make a vector of times from start to end dates
        etOne = spice.str2et(start_date)
        etTwo = spice.str2et(end_date)
        self.times = [x * (etTwo - etOne) / n_steps + etOne for x in range(n_steps)]
        self.dt=self.times[1]-self.times[0]

        self.j2=False
        if "j2" in perts:
            self.j2=True

        self.srp=False
        if "srp" in perts:
            self.vectors_from_sun, lightTimes =spice.spkpos(self.cb['name'], self.times, 'J2000', 'NONE', "SUN")
            self.srp=True

        self.pert_bodies = True
        self.bodies=[]
        #check if list of perturbing bodies is present then store it.
        for pert in perts:
            if type(pert) is list:
                self.pert_bodies=True
                self.bodies=pert # list of pd objects

        self.body_names=[]
        self.r_cb2nb_lst=[]
        if self.bodies:
            for body in self.bodies:
                self.body_names.append(body['name']) #for future debugging
                r_cb2nb, lightTimes = spice.spkpos(body['name'], self.times, 'J2000', 'NONE', cb['name'])
                self.r_cb2nb_lst.append(r_cb2nb)

        self.sc_data=spacecraft_data


    def propagate(self):
        print("propagating...")
        #self.n_steps=int(np.ceil((self.tspan/self.dt)))

        self.ts = np.zeros((self.n_steps,1))
        self.ys = np.zeros((self.n_steps,6))
        self.y0 = self.r0+self.v0
        self.ts[0]=0
        self.ys[0]=self.y0
        self.steps=1

        self.solver=ode(self.diffy_q)
        self.solver.set_integrator(self.integrator)
        self.solver.set_initial_value(self.y0, 0)

        while self.solver.successful() and self.steps < self.n_steps:
            self.solver.integrate((self.solver.t + self.dt))
            self.ts[self.steps]=self.solver.t
            self.ys[self.steps]=self.solver.y
            self.steps += 1

        self.rs=self.ys[:,:3]
        self.vs = self.ys[:,3:]

    def diffy_q(self, t, y):
        rx,ry,rz,vx,vy,vz=y
        r=np.array([rx,ry,rz])
        v=np.array([vx,vy,vz])

        norm_r=np.linalg.norm(r)

        ax,ay,az=-r*self.cb['mu']/norm_r**3

        a=[ax,ay,az]

        if self.j2:
            j2=self.cb['j2']

            z2 = r[2] ** 2
            r2 = norm_r ** 2
            tx = r[0] / norm_r * (5 * z2 / r2 - 1)
            ty = r[1] / norm_r * (5 * z2 / r2 - 1)
            tz = r[2] / norm_r * (5 * z2 / r2 - 3)
            a += 1.5 * self.cb['j2'] * self.cb['mu'] * self.cb['radius'] ** 2 / norm_r ** 4 * np.array([tx, ty, tz])

        if self.srp:
            d=self.vectors_from_sun[self.steps, :]+r #vector from sun to spacecraft
            g1=1e8
            beta=((1+self.sc_data['reflectance'])*g1)/(self.sc_data['mass']/self.sc_data['area'])
            g=beta/((np.linalg.norm(d)**2)*1000)
            #print("before: ",a)
            a += g*d
            #print("after: ",a)
            #source: https://ntrs.nasa.gov/api/citations/20080012725/downloads/20080012725.pdf


        if self.pert_bodies:
            for i in range(len(self.bodies)):
                mu=self.bodies[i]['mu']
                r_sat2body = self.r_cb2nb_lst[i][self.steps, :] - r
                a += mu * ((r_sat2body / (np.linalg.norm(r_sat2body) ** 3)) - (self.r_cb2nb_lst[i][self.steps, :] / (np.linalg.norm(self.r_cb2nb_lst[i][self.steps, :]) ** 3)))

        return [vx, vy, vz, a[0], a[1], a[2]]

    def calculate_coes(self):
        print("calculating...")
        coes=[]
        for time in range(len(self.rs)):
            coes.append(t.rv2coes(self.rs[time],self.vs[time], deg=True))
        self.coes=np.array(coes)

    def plot_coes(self, hours=False, days=False):
        print("plotting...")

        ts=self.ts
        xlabel = "seconds"

        if hours:
            ts=ts/3600
            xlabel = "hours"

        if days:
            ts=ts/(3600*24)
            xlabel = "days"

        fig,axs=plt.subplots(2,3)

        # plot semi major axis
        axs[0,0].plot(ts, self.coes[:,0])
        axs[0,0].set_title('semi-major axis vs time')
        axs[0, 0].grid(True)
        axs[0, 0].set_ylabel("a (km)")
        axs[0, 0].set_xlabel(xlabel)

        # plot eccentricity
        axs[0, 1].plot(ts, self.coes[:, 1])
        axs[0, 1].set_title('eccentricity vs time')
        axs[0, 1].grid(True)
        axs[0, 1].set_ylabel("e")
        axs[0, 1].set_xlabel(xlabel)

        # plot inclination
        axs[0, 2].plot(ts, self.coes[:, 2])
        axs[0, 2].set_title('inclination vs time')
        axs[0, 2].grid(True)
        axs[0, 2].set_ylabel("i (deg)")
        axs[0, 2].set_xlabel(xlabel)

        # plot true anomaly
        axs[1, 0].plot(ts, self.coes[:, 3])
        axs[1, 0].set_title('true anomaly vs time')
        axs[1, 0].grid(True)
        axs[1, 0].set_ylabel("nu (deg)")
        axs[1, 0].set_xlabel(xlabel)

        # plot argument of periapsis
        #axs[1, 1].plot(ts, self.coes[:, 4])
        z = np.polyfit(ts.flatten(), self.coes[:, 4] , 1)
        p = np.poly1d(z)
        axs[1, 1].plot(ts, p(ts) )
        axs[1, 1].set_title('apsidal rotation')
        axs[1, 1].grid(True)
        axs[1, 1].set_ylabel("aop (deg)")
        axs[1, 1].set_xlabel(xlabel)

        # plot longitude of ascending node
        axs[1, 2].plot(ts, self.coes[:, 5])
        axs[1, 2].set_title('nodal regression')
        axs[1, 2].grid(True)
        axs[1, 2].set_ylabel("lan (deg)")
        axs[1, 2].set_xlabel(xlabel)

        plt.show()