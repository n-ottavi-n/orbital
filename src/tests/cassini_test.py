import spiceypy as spice
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import tools as t
import planetary_data


spice.furnsh('spice_cassini/cassini_kernel.txt')

step = 4000
# we are going to get positions between these two dates
utc = ['Jun 20, 2004', 'Dec 1, 2005']

# get et values one and two, we could vectorize str2et
etOne = spice.str2et(utc[0])
etTwo = spice.str2et(utc[1])
#print("ET One: {}, ET Two: {}".format(etOne, etTwo))
dt=(etTwo-etOne)/step

# get times
times = [x*(etTwo-etOne)/step + etOne for x in range(step)]

# check the documentation on spkpos before continuing
#help(spice.spkpos)

#Run spkpos as a vectorized function
positions_cass, lightTimes = spice.spkpos('Cassini', times, 'J2000', 'NONE', 'SATURN BARYCENTER')
positions_d, lightTimes = spice.spkpos('DIONE', times, 'J2000', 'NONE', 'SATURN BARYCENTER')
positions_e, lightTimes = spice.spkpos('ENCELADUS', times, 'J2000', 'NONE', 'SATURN BARYCENTER')
positions_t, lightTimes = spice.spkpos('TITAN', times, 'J2000', 'NONE', 'SATURN BARYCENTER')

positions=[positions_cass,positions_d,positions_e,positions_t]

# Clean up the kernels
spice.kclear()

labels=['cassini','dione','enceladus','titan']

t.plot_n_orbits_animate(positions, dt, labels, cb=planetary_data.Saturn, show_plot=False, save=True,  au_units=False, interval=0.1)
#t.plot_3d(positions, cb=planetary_data.Saturn, show_plot=True)