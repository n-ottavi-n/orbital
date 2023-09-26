import numpy as np
import tools as t
from Propagator import Propagator
import planetary_data as pd
import spiceypy as spice


class interplanetary_interface:

    def __init__(self,data, bodies, start_date, end_date, perturbations,central_body=pd.sun, sc_data={}, steps=1000,):
        self.data=data,
        self.bodies=bodies,
        self.start_date=start_date
        self.end_date=end_date
        self.perturbations=perturbations,
        self.cb=central_body,
        self.sc_data=sc_data,
        self.steps=steps
        spice.furnsh('../spice_solar_system/solar_system_kernel.txt')

    def plot_trajectory(self, steps=1000, animate=False, show=False, save=False, save_file='matplot003.gif'):
        '''
        An interface for plotting solar system trajectories from state vectors, plots can be animated and saved
        @param data: [[labels], [[x,y,z,vx,vy,vz]]] initial conditions
        @param bodies: bodies to plot using NASA SPICE data
        @param start_date: start of propagation
        @param end_date: end of propagation
        @param perturabtions: list of perturbations to account for ['j2', 'srp', [pd.jupiter,pd.moon,....]] bodies must be plantery_data objects
        @param central_body: planetary_data object
        @param steps: number of time steps for propagation
        @param animate: if True makes an animated 3d plot
        @param show: if True shows the plot
        @param save: if True saves the file under the name save_file
        @param save_file: see above
        @return: None
        '''
        # get names and state vectors from list of spacecraft
        labels = self.data[0][0]
        states = self.data[0][1]
        #initialize propagators array
        props=[]

        etOne = spice.str2et(self.start_date)
        etTwo = spice.str2et(self.end_date)
        dt = (etTwo - etOne) / self.steps

        # initilize times array
        times = [x * (etTwo - etOne) / steps + etOne for x in range(self.steps)]  # [start_date : end_date] len=steps

        cb=self.cb[0] # I have to but whyyyyyyy ????
        #loop through propagators
        for state in states:
            prop = Propagator(state, self.steps, start_date=self.start_date, end_date=self.end_date, coes=False, deg=True,cb=cb, perts=self.perturbations, spacecraft_data=self.sc_data)
            times = prop.times
            prop.propagate()
            props.append(prop.rs)

        obs = 'SOLAR SYSTEM BARYCENTER'

        states_bodies=[]
        #get states of other bodies to plot e.g. perturbating bodies
        bodies=self.bodies[0]
        for body in bodies:
            positions, lightTimes = spice.spkpos(body, times, 'ECLIPJ2000', 'NONE', obs)
            # append states to propagator states
            props.append(positions)
            # append other bodies labels to satellite names
            labels.append(body)

        self.props=np.array(props)


        #set distance to AU if central_body=sun
        au=False
        if cb['name']=='sun':
            au=True

        #set plot title
        plot_title=''

        #call the required plot
        if animate:
            t.plot_n_orbits_animate(self.props, step_t=dt, labels=labels, cb=cb, show_plot=show, save=save, au_units=au, save_file=save_file)
        else:
            t.plot_n_orbits(self.props, step_t=dt, labels=labels, cb=cb, show_plot=show, save_plot=save, au_units=au, save_file=save_file)



