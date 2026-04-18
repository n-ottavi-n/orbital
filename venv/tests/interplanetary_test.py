import planetary_data
from interplanetary_interface import interplanetary_interface
import spiceypy as spice

spice.kclear()
spice.furnsh('../spice_solar_system/solar_system_kernel.txt')

perturbations=[planetary_data.earth,planetary_data.venus, planetary_data.jupiter]
bodies=['earth','venus']

#t0='Sep 24, 2023, 00:00 UTC'
#tf='Sep 24, 2024, 00:00 UTC'
t0='Oct 26, 2029, 04:00 UTC'
tf='Apr 04, 2030, 00:00 UTC'
cb=planetary_data.sun
states=[['venus 2029'],[[1.25272836e+08,  8.01426306e+07, -6.20899724e+03, -14.98060275,  22.67518145,  -0.48325323]]]
#states=[['mars 2020'],[[9.280670739824386E+07, -1.186408488879065E+08 , 8.107387848985940E+04, 2.685659452407585E+01 , 1.945466108357662E+01 , 1.330899689848824E+00]]]
#states=[['cassini'],[[-6.249878234110750E+07 ,-9.026835124549742E+07 , 2.336516112486362E+06, 3.736602844997936E+01 ,-1.582751157992247E+01 ,-2.364222029667019E+00]]]
#states=[[],[]]

sc_data={
    'mass':5000, #kg
    'reflectance':0,
    'area':12 #m²
}



int=interplanetary_interface(states,bodies,t0,tf, perturbations, central_body=cb, sc_data=sc_data, steps=10000)
int.plot_trajectory(animate=False,show=True, save=False)