import numpy as np
import matplotlib.pyplot as plt
import planetary_data
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D
from Propagator import Propagator
import tools as t
from math import sqrt


cb=planetary_data.earth

##############ORBIT 1##############

r_mag=cb['radius']+400 #orbit radius at t0
v_mag=sqrt(cb['mu']/r_mag) # orbit velocity at t0

r00=[3856.19559189663, 6505.740175761596, 3098.372955020122]
v00=[-2.3843272769197865, 5.1330459987985995, -4.468810000455476]

state00=r00+v00
label0='iss'

##############ORBIT 2##############

r_mag=cb['radius']+3850 #orbit radius at t0
v_mag=sqrt(cb['mu']/r_mag) # orbit velocity at t0

r01=[r_mag,0,0]
v01=[0,v_mag,0]

state01=r01+v01
label1='moon'
####################################

#integration parameters
n_days=1.0
tspan=3600*24*n_days
dt=100

#propagation
prop0=Propagator(state00,tspan, dt, cb=cb)
prop1=Propagator(state01,tspan, dt, cb=cb)


prop0.propagate()
prop1.propagate()

t.plot_n_orbits([prop0.rs,prop1.rs], labels=[label0,label1],cb=cb, show_plot=True)
