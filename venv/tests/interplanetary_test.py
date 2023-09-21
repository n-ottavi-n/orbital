import planetary_data
from plot_interplanetary_interface import plot_trajectory

perturbations=[[planetary_data.earth,planetary_data.mars, planetary_data.jupiter]]
bodies=['earth','mars']

t0='Jul 31, 2020, 00:00 UTC'
tf='Feb 17, 2021, 00:00 UTC'
cb=planetary_data.sun

states=[['mars 2020'],[[9.280670739824386E+07, -1.186408488879065E+08 , 8.107387848985940E+04, 2.685659452407585E+01 , 1.945466108357662E+01 , 1.330899689848824E+00]]]
#states=[['cassini'],[[-6.249878234110750E+07 ,-9.026835124549742E+07 , 2.336516112486362E+06, 3.736602844997936E+01 ,-1.582751157992247E+01 ,-2.364222029667019E+00]]]

sc_data={
    'mass':5000, #kg
    'reflectance':0.2,
    'area':12 #m²
}




plot_trajectory(states,bodies,t0,tf, perturbations, central_body=cb, sc_data=sc_data, steps=100000, animate=False,show=True, save=False)