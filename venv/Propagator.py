import numpy as np
import matplotlib.pyplot as plt
import planetary_data
from scipy.integrate import ode
import tools as t
import planetary_data as pd
import matplotlib.animation as animation
import spiceypy as spice

def null_perts():
    return {
        'j2':False,
        'aero':False,
        'moon':False,
        'sun':False
    }


class Propagator:
    '''
    inputs:
        if cartesian: [x0,y0,z0,vx0,vy0,vz0]
        if keplerian: a,e,i,m,aop,raan
    '''
    def __init__(self, state0, n_steps, start_date='Sep 16, 2023, 00:00 UTC', end_date='Sep 17, 2023, 00:00 UTC', coes=False, deg=False, cb=pd.earth, perts=null_perts()):
        '''
        @param state0: if cartesian: [x0,y0,z0,vx0,vy0,vz0]  if keplerian: [a,e,i,m,aop,raan]
        @param n_steps: number of timesteps for integration
        @param start_date: format:  '%b %d, %Y, %H:%M %Z'
        @param end_date: format:  '%b %d, %Y, %H:%M %Z'
        @param coes: bool: True if state0 is keplerian elements
        @param deg: bool: True if input keplerian elements in degrees
        @param cb: central body
        @param perts: perturbations to account for
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

        #make a vector of times from start to end dates
        etOne = spice.str2et(start_date)
        etTwo = spice.str2et(end_date)
        self.times = [x * (etTwo - etOne) / n_steps + etOne for x in range(n_steps)]
        self.dt=self.times[1]-self.times[0]

        if self.perts['moon']:
            self.mu_moon = planetary_data.moon['mu']
            #get moon position vectors from start to end dates
            self.r_cb2nb, lightTimes = spice.spkpos('MOON', self.times, 'J2000', 'NONE', 'EARTH BARYCENTER')


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
        self.solver.set_integrator('dopri5')
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

        if self.perts['j2']:
            j2=self.cb['j2']

            #r7=norm_r**7
            #x2=rx**2
            #y2=ry**2
            #z2=rz**2

            #ax += j2 * (rx / r7) * (6 * z2 - 1.5 * (x2 + y2))
            #ay += j2 * (ry / r7) * (6 * z2 - 1.5 * (x2 + y2))
            #az += j2 * (rz / r7) * (3 * z2 - 4.5 * (x2 + y2))

            z2 = r[2] ** 2
            r2 = norm_r ** 2
            tx = r[0] / norm_r * (5 * z2 / r2 - 1)
            ty = r[1] / norm_r * (5 * z2 / r2 - 1)
            tz = r[2] / norm_r * (5 * z2 / r2 - 3)
            a += 1.5 * self.cb['j2'] * self.cb['mu'] * self.cb['radius'] ** 2 / norm_r ** 4 * np.array([tx, ty, tz])

        if self.perts['moon']:
            r_sat2body=self.r_cb2nb[self.steps, :]-r
            a+=self.mu_moon*((r_sat2body/(np.linalg.norm(r_sat2body)**3))-(self.r_cb2nb[self.steps,:]/(np.linalg.norm(self.r_cb2nb[self.steps,:])**3)))


        #return [vx,vy,vz,ax,ay,az]
        return [vx, vy, vz, a[0], a[1], a[2]]

    def calculate_coes(self):
        print("calculating...")
        coes=[]
        for time in range(len(self.rs)):
            coes.append(t.rv2coes(self.rs[time],self.vs[time], deg=True))
        self.coes=np.array(coes)

    def plot_3d_animate(self, show_plot=False, save_plot=False, au_units=False):
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, projection='3d')
        ax.set_aspect('equal')

        # plot central body

        radius=self.cb['radius']

        if au_units:
            radius=radius/1.5e8
            self.rs=self.rs/1.5e8

        _u = np.linspace(0, 2 * np.pi, 100)
        _v = np.linspace(0, np.pi, 100)
        _x = radius * np.outer(np.cos(_u), np.sin(_v))
        _y = radius * np.outer(np.sin(_u), np.sin(_v))
        _z = radius * np.outer(np.ones(np.size(_u)), np.cos(_v))
        ax.plot_surface(_x, _y, _z, color='linen', alpha=0.5)

        # plot circular curves over the surface
        theta = np.linspace(0, 2 * np.pi, 100)
        z = np.zeros(100)
        x = radius * np.sin(theta)
        y = radius * np.cos(theta)

        ax.plot(x, y, z, color='black', alpha=0.75)
        ax.plot(z, x, y, color='black', alpha=0.75)

        ## add axis lines
        zeros = np.zeros(1000)
        line = np.linspace(-10, 10, 1000)

        ax.plot(line, zeros, zeros, color='black', alpha=0.75)
        ax.plot(zeros, line, zeros, color='black', alpha=0.75)
        ax.plot(zeros, zeros, line, color='black', alpha=0.75)

        # plot x y z vectors
        l = radius * 2

        max_val = np.max(np.abs(self.rs))

        ax.set_xlim([-max_val, max_val])
        ax.set_ylim([-max_val, max_val])
        ax.set_zlim([-max_val, max_val])

        ax.set_xlabel('X (km)')
        ax.set_ylabel('Y (km)')
        ax.set_zlabel('Z (km)')

        if au_units:
            ax.set_xlabel('X (au)')
            ax.set_ylabel('Y (au)')
            ax.set_zlabel('Z (au)')


        def update(num, data, line2, r_current):
            print("updating...")
            #plot trajectory
            line2.set_data(data[:num,:2].transpose())
            line2.set_3d_properties(data[:num,2])
            #plot current position
            r_current.set_data(data[num,:2].transpose())
            r_current.set_3d_properties(data[num, 2])

        N = self.n_steps
        data = self.rs

        ax.plot([self.rs[0, 0]], [self.rs[0, 1]], [self.rs[0, 2]], 'wo', label='initial position', zorder=5)
        r_current, =ax.plot([self.rs[0, 0]], [self.rs[0, 1]], [self.rs[0, 2]], 'yo', label='current position', zorder=5)

        #.plot([self.rs[-1, 0]], [self.rs[-1, 1]], [self.rs[-1, 2]], 'yo', label='current position', zorder=6)

        line2, = ax.plot([self.rs[0, 0]], [self.rs[0, 1]], [self.rs[0, 2]], label='spacecraft',zorder=4)  # initialize

        plt.legend()

        ani = animation.FuncAnimation(fig, update, N, fargs=(data, line2, r_current), interval=0.2, blit=False)

        if show_plot:
            plt.show()
        if save_plot:
            plt.savefig(title+'.png')

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