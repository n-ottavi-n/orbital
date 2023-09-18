import tools as t
from Propagator import Propagator
from Propagator import null_perts
import planetary_data as pd
import spiceypy as spice

perturbations=['j2','moon']



def plot_orbits(satellite_names, bodies, start_date, end_date, perturabtions,central_body=pd.earth, steps=400, animate=False, show=False, save=False, save_file='matplot003.gif'):
    '''
    An interface for plotting 3d orbits of any number of satellites and other bodies if they have a name and their tle is in a file in the /venv/data folder
    plots can be animated and saved
    @param satellite_names: list of sat names as they appear in /venv/data txt files
    @param bodies: bodies to plot using NASA SPICE data
    @param start_date: start of propagation
    @param end_date: end of propagation
    @param perturabtions: list of perturbations to account for ['j2', 'moon', 'srp', 'jupiter',....] terms must be keys in planetary data dicts
    @param central_body: obvious
    @param steps: number of time steps for propagation
    @param animate: if True makes an animated 3d plot
    @param show: if True shows the plot
    @param save: if True saves the file under the name save_file
    @param save_file: see above
    @return: None
    '''
    # get list of tles and labels from satellite_names
    labels, states = t.get_sats_from_file(satellite_names)

    #initialize propagators array
    props=[]

    #initilize times array
    times=[] # [start_date : end_date] len=steps

    #get perturbations
    perts = null_perts()
    for p in perturbations:
        perts[p]=True


    #loop through propagators
    for state in states:
        prop = Propagator(state, steps, start_date=start_date, end_date=end_date, coes=True, deg=True, perts=perts)
        dt = prop.dt
        times = prop.times
        prop.propagate()
        props.append(prop.rs)

    #spice.furnsh
    spice.furnsh('spice_lunar/earth_moon_kernel.txt') #modify for more genricity

    states_bodies=[]
    #get states of other bodies to plot e.g. perturbating bodies
    for body in bodies:
        positions, lightTimes = spice.spkpos(body, times, 'J2000', 'NONE', central_body['name'])
        # append states to propagator states
        props.append(positions)
        # append other bodies labels to satellite names
        labels.append(body)

    spice.kclear()

    #set distance to AU if central_body=sun
    au=False
    if central_body['name']=='sun':
        au=True

    #set plot title
    plot_title=''

    #call the required plot
    if animate:
        t.plot_n_orbits_animate(props, step_t=dt, labels=labels, cb=central_body, show_plot=show, save=save, save_file=save_file)
    else:
        t.plot_n_orbits(props, labels=labels, cb=central_body, show_plot=show, save_plot=save, au_units=au) #update function
