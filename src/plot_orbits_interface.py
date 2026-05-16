import numpy as np
import tools as t
from Propagator import Propagator
import planetary_data as pd
import spiceypy as spice


def plot_orbits(satellites, bodies, start_date, end_date, perturbations,central_body=pd.earth, sc_data={}, steps=400, from_file=True, animate=False, show=False, save=False, save_file='matplot003.gif', frame_step=10, plot_coes=False):
    '''
    An interface for plotting 3d orbits of any number of satellites and other bodies if they have a name and their tle is in a file in the /venv/data folder,
    plots can be animated and saved
    @param satellites: list of satellite names as they appear in /venv/data txt files or [[string satname, a, e, i, M0, argp, raan]] if from_file=False
    @param bodies: bodies to plot using NASA SPICE data
    @param start_date: start of propagation
    @param end_date: end of propagation
    @param perturabtions: list of perturbations to account for ['j2', 'srp', [pd.jupiter,pd.moon,....]] bodies must be plantery_data objects
    @param central_body: planetary_data object
    @param steps: number of time steps for propagation
    @param from_file: if True, satellites are loaded from files, else they must be given as coes 
    @param animate: if True makes an animated 3d plot
    @param show: if True shows the plot
    @param save: if True saves the file under the name save_file
    @param save_file: see above
    @return: None
    '''
    # get list of coes and labels from satellite_names
    if from_file:
        labels, states = t.get_sats_from_file(satellites, start_date)
    else:
        states = []
        labels = [sat[0] for sat in satellites]
        for sat in satellites:
            coes = sat[1:7]
            states.append(coes)
    props=[]

    if central_body['name']=='earth':
        t.load_earth_moon_kernels()
    else:
        t.load_solar_system_kernels()

    etOne = spice.str2et(start_date)
    etTwo = spice.str2et(end_date)
    dt = (etTwo - etOne) / steps

    # initilize times array
    times = [x * (etTwo - etOne) / steps + etOne for x in range(steps)]  # [start_date : end_date] len=steps

    print("central body: ", central_body['name'])
    #loop through propagators
    for state in states:
        prop = Propagator(state, steps, start_date=start_date, end_date=end_date, coes=True, deg=True, cb=central_body, perts=perturbations, spacecraft_data=sc_data)
        times = prop.times
        prop.propagate()
        props.append(prop.rs)
        if plot_coes:
            prop.calculate_coes()
            prop.plot_coes(days=True)

    #print(props)


    obs = central_body['name']

    states_bodies=[]
    #get states of other bodies to plot e.g. perturbating bodies
    for body in bodies:

        positions, lightTimes = spice.spkpos(body, times, 'J2000', 'NONE', obs)
        # append states to propagator states
        props.append(positions)
        # append other bodies labels to satellite names
        labels.append(body)

    props=np.array(props)

    spice.kclear()

    au = False
    #set distance to AU if central_body=sun
    if obs == 'sun':
        au = True

    #set plot title
    plot_title=''

    #call the required plot
    #return t.plot_n_orbits(props, step_t=dt, labels=labels, cb=central_body, show_plot=show, save_plot=save, au_units=au, save_file=save_file)
  
    if animate:
        t.plot_n_orbits_animate(props, frame_step=frame_step, step_t=dt, labels=labels, cb=central_body, show_plot=show, save=save, au_units=au, save_file=save_file)
    else:
        t.plot_n_orbits(props, step_t=dt, labels=labels, cb=central_body, show_plot=show, save_plot=save, au_units=au, save_file=save_file)


