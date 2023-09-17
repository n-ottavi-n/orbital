import tools as t
from tools import ecc_anom
from tools import coes2rv
import math
import planetary_data
import numpy as np
from Propagator import Propagator
from Propagator import null_perts
import matplotlib.animation as animation
import spiceypy as spice
from datetime import datetime

cb=planetary_data.earth

j2perts=null_perts()
j2perts['j2']=True
j2perts['moon']=True


spice.furnsh('spice_lunar/earth_moon_kernel.txt')


t0='Sep 17, 2023, 00:00 UTC'
tf='Sep 18, 2023, 00:00 UTC'
states, names = t.n_tle2coes("data/molniya.txt", 3, t0)


states0=np.array(states, dtype=float)
labels=np.array(names)

props=[]
times=[]

for state in states0:
    prop=Propagator(state,400,start_date=t0, end_date=tf, coes=True,deg=True, perts=j2perts)
    dt=prop.dt
    times=prop.times
    prop.propagate()
    props.append(prop.rs)
    #prop.calculate_coes()
    #prop.plot_coes(hours=True)
    #prop.plot_3d_animate(show_plot=True)





positions_m, lightTimes = spice.spkpos('MOON', times, 'J2000', 'NONE', 'EARTH')

positions_m=np.array(positions_m)
props.append(positions_m)

labels=np.append(labels,'moon')

t.plot_n_orbits_animate(props, step_t=dt, labels=labels,cb=cb, show_plot=True, save=False)
