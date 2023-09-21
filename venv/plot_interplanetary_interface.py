import numpy as np
import tools as t
from Propagator import Propagator
import planetary_data as pd
import spiceypy as spice


def plot_trajectory(data, bodies, start_date, end_date, perturbations,central_body=pd.sun, sc_data={}, steps=1000, animate=False, show=False, save=False, save_file='matplot003.gif'):
    '''
    An interface for plotting solar system trajectories from state vectors, plots can be animated and saved
    @param data: [[labels], [[x,y,z,vx,vy,vz]]]
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
    # get list of coes and labels from satellite_names
    labels = data[0]
    states = data[1]
    #initialize propagators array
    props=[]

    spice.furnsh('../spice_solar_system/solar_system_kernel.txt')

    etOne = spice.str2et(start_date)
    etTwo = spice.str2et(end_date)
    dt = (etTwo - etOne) / steps

    # initilize times array
    times = [x * (etTwo - etOne) / steps + etOne for x in range(steps)]  # [start_date : end_date] len=steps

    #loop through propagators
    for state in states:
        prop = Propagator(state, steps, start_date=start_date, end_date=end_date, coes=False, deg=True,cb=central_body, perts=perturbations, spacecraft_data=sc_data)
        times = prop.times
        prop.propagate()
        props.append(prop.rs)

    obs = 'SOLAR SYSTEM BARYCENTER'


    states_bodies=[]
    #get states of other bodies to plot e.g. perturbating bodies
    for body in bodies:

        positions, lightTimes = spice.spkpos(body, times, 'ECLIPJ2000', 'NONE', obs)
        # append states to propagator states
        props.append(positions)
        # append other bodies labels to satellite names
        labels.append(body)

    props=np.array(props)


    spice.kclear()

    #set distance to AU if central_body=sun
    au=False
    if central_body['name']=='sun':
        au=True

    #set plot title
    plot_title=''

    #call the required plot
    if animate:
        t.plot_n_orbits_animate(props, step_t=dt, labels=labels, cb=central_body, show_plot=show, save=save, au_units=au, save_file=save_file)
    else:
        t.plot_n_orbits(props, step_t=dt, labels=labels, cb=central_body, show_plot=show, save_plot=save, au_units=au, save_file=save_file)


