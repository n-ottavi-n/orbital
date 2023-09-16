import spiceypy as spice
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import tools as t
import planetary_data


spice.furnsh('spice_lunar/earth_moon_kernel.txt')

step = 400
# we are going to get positions between these two dates
utc = ['Jan 26, 1994, 22:00 UTC', 'Feb 02, 1994']


# get et values one and two, we could vectorize str2et
etOne = spice.str2et(utc[0])
etTwo = spice.str2et(utc[1])

dt=(etTwo-etOne)/step
# get times
times = [x*(etTwo-etOne)/step + etOne for x in range(step)]

positions_m, lightTimes = spice.spkpos('MOON', times, 'J2000', 'NONE', 'EARTH BARYCENTER')
print(positions_m)
#perturbation
perts=['moon']

positions=np.array([positions_m])
labels=['moon']
# Clean up the kernels
spice.kclear()

#t.plot_n_orbits_animate(positions, dt, labels, cb=planetary_data.earth, show_plot=True,save=False, interval=0.01)