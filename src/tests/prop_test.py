import numpy as np
import matplotlib.pyplot as plt
import planetary_data
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D
from Propagator import Propagator
import tools as t


cb=planetary_data.earth

r_mag=cb['radius']+3000000
v_mag=np.sqrt(cb['mu']/r_mag)

r0=[r_mag,r_mag*0.1,-r_mag*0.1]
v0=[0, v_mag/2, v_mag*0.3]
state0=r0+v0

tspan=3600*24*1.0

dt=1

#keplerian elements in degrees
rp=35757.0
ra=35814
a=cb['radius']+(rp+ra)/2
e=0.0006665
i=2.82
m=306.78 #mean anomaly at epoch
aop=214.65
raan=288.76

state00=[a,e,i,m,aop,raan]
label0='satellite'

state01=[cb['radius'],0,0,0,0,0]
label1="equator"

#propagation
prop0=Propagator(state00,tspan, dt, coes=True, deg=True, cb=cb)
prop1=Propagator(state01,tspan, dt, coes=True, deg=True, cb=cb)


prop0.propagate()
prop1.propagate()

t.plot_n_orbits([prop0.rs,prop1.rs], labels=[label0,label1],cb=cb, show_plot=True)
'''
prop=Propagator(state0,tspan, dt, coes=True, deg=True, cb=cb)
prop.propagate()
prop.plot_3d(show_plot=True)

'''