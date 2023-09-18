import spiceypy as spice
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import tools as t
import planetary_data


spice.furnsh('spice_solar_system/solar_system_kernel.txt')

step = 400
# we are going to get positions between these two dates
utc = ['Sep 6, 1977', 'Sep 10, 1989']

# get et values one and two, we could vectorize str2et
etOne = spice.str2et(utc[0])
etTwo = spice.str2et(utc[1])
#
dt=(etTwo-etOne)/step
# get times
times = [x*(etTwo-etOne)/step + etOne for x in range(step)]


# check the documentation on spkpos before continuing
#help(spice.spkpos)

#Run spkpos as a vectorized function
positions_e, lightTimes = spice.spkpos('EARTH BARYCENTER', times, 'J2000', 'NONE', 'SUN')
positions_v, lightTimes = spice.spkpos('VENUS BARYCENTER', times, 'J2000', 'NONE', 'SUN')
positions_j, lightTimes = spice.spkpos('JUPITER BARYCENTER', times, 'J2000', 'NONE', 'SUN')
positions_s, lightTimes = spice.spkpos('SATURN BARYCENTER', times, 'J2000', 'NONE', 'SUN')
positions_u, lightTimes = spice.spkpos('URANUS BARYCENTER', times, 'J2000', 'NONE', 'SUN')
positions_n, lightTimes = spice.spkpos('NEPTUNE BARYCENTER', times, 'J2000', 'NONE', 'SUN')
positions_vgr1, lightTimes = spice.spkpos('VOYAGER 1', times, 'J2000', 'NONE', 'SUN')
positions_vgr2, lightTimes = spice.spkpos('VOYAGER 2', times, 'J2000', 'NONE', 'SUN')

positions=np.array([positions_e, positions_v, positions_j, positions_s, positions_u, positions_n, positions_vgr1, positions_vgr2])
labels=['earth','venus','jupiter','saturn','uranus','neptune','vger1','vger2']
# Clean up the kernels
spice.kclear()

t.plot_n_orbits_animate(positions, dt, labels, cb=planetary_data.sun, show_plot=True,save=False,  au_units=True, interval=0.01)