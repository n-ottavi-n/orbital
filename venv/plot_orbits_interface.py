import numpy as np
import tools as t
from Propagator import Propagator
import planetary_data as pd
import spiceypy as spice


def plot_orbits(satellite_names, bodies, start_date, end_date, perturbations,central_body=pd.earth, steps=400, animate=False, show=False, save=False, save_file='matplot003.gif'):
    '''
    An interface for plotting 3d orbits of any number of satellites and other bodies if they have a name and their tle is in a file in the /venv/data folder,
    plots can be animated and saved
    @param satellite_names: list of sat names as they appear in /venv/data txt files
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
    labels, states = t.get_sats_from_file(satellite_names, start_date)

    #initialize propagators array
    props=[]

    '''
    #get perturbations
    perts = null_perts()
    for p in perturbations:
        perts[p]=True
        '''

    if central_body['name']=='earth':
        spice.furnsh('../spice_lunar/earth_moon_kernel.txt')
    else:
        spice.furnsh('../spice_solar_system/solar_system_kernel.txt')

    etOne = spice.str2et(start_date)
    etTwo = spice.str2et(end_date)
    dt = (etTwo - etOne) / steps

    # initilize times array
    times = [x * (etTwo - etOne) / steps + etOne for x in range(steps)]  # [start_date : end_date] len=steps

    #loop through propagators
    for state in states:
        prop = Propagator(state, steps, start_date=start_date, end_date=end_date, coes=True, deg=True, perts=perturbations)
        times = prop.times
        prop.propagate()
        props.append(prop.rs)



    states_bodies=[]
    #get states of other bodies to plot e.g. perturbating bodies
    for body in bodies:
        positions, lightTimes = spice.spkpos(body, times, 'J2000', 'NONE', central_body['name'])
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


