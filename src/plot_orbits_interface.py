import numpy as np
import tools as t
from Propagator import Propagator
import planetary_data as pd
import spiceypy as spice


def plot_orbits(satellites, bodies, start_date, end_date, perturbations,central_body=pd.earth, sc_data={}, steps=400, from_file=True, plot_3d=False, animate=False, show=False, save=False, save_file='matplot003.gif', frame_step=10, plot_coes=False, plot_distance_body=None):
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
        if plot_distance_body:
            prop.plot_distance_to_body(plot_distance_body, days=True, title=f"{plot_distance_body} Probe — Distance to {plot_distance_body}",
            save_html=f"distance_{plot_distance_body.lower()}.html")
        if plot_coes:
            prop.calculate_coes()
            prop.plot_coes_plotly(days=True)

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

    # rotate into central body equatorial frame
    if central_body['name'] != 'earth' and central_body['name'] != 'sun':
        iau_frame = f"IAU_{central_body['name'].upper()}"
        epoch_mid = times[len(times)//2]   # use midpoint epoch for rotation
        try:
            tipm = spice.pxform("J2000", iau_frame, epoch_mid)
            rotated = []
            for traj in props:
                # traj is Nx3 — rotate each position vector
                rotated.append(np.array([tipm @ r for r in traj]))
            props = rotated
        except Exception as e:
            print(f"Warning: could not rotate to {iau_frame}: {e}")
            print("Plotting in J2000 frame")

    props=np.array(props)


    au = False
    #set distance to AU if central_body=sun
    if obs == 'sun':
        au = True

    #set plot title
    # build title from parameters
    sat_names  = ", ".join([s[0] for s in satellites]) if not from_file else ", ".join(satellites)
    body_names = ", ".join(bodies) if bodies else "none"
    pert_names = []
    for p in perturbations:
        if isinstance(p, list):
            pert_names += [b['name'].capitalize() for b in p if isinstance(b, dict)]
        elif isinstance(p, str) and p not in ['j2', 'srp']:
            pert_names.append(p.capitalize())
    pert_str = ", ".join(pert_names) if pert_names else "none"

    title = (
        f"{central_body['name'].capitalize()} system — "
        f"{sat_names}<br>"
        f"<sup>bodies: {body_names}   "
        f"perturbations: {pert_str}   "
        f"{start_date} → {end_date}</sup>"
    )


    #call the required plot
    #return t.plot_n_orbits(props, step_t=dt, labels=labels, cb=central_body, show_plot=show, save_plot=save, au_units=au, save_file=save_file)
    if plot_3d:
        if animate:
            t.plot_n_orbits_animate(props, frame_step=frame_step, step_t=dt, labels=labels, cb=central_body, show_plot=show, save=save, au_units=au, save_file=save_file)
        else:
            t.plot_n_orbits(props, step_t=dt, labels=labels, cb=central_body, show_plot=show, save_plot=save, au_units=au, save_file=save_file)
    else:
        t.plot_2d(props, labels, central_body, times=times, save_html=save_file if save else None, title=title)

    spice.kclear()


