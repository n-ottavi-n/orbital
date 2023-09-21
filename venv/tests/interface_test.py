import planetary_data
from plot_orbits_interface import plot_orbits

#sats=['CHINASAT 18','YAOGAN-40 01C','BEIDOU 3 G4']
sats=['CHINASAT 18']
#perturbations=[[planetary_data.mars]]
perturbations=['j2', [planetary_data.moon]]
bodies=['moon']
#bodies=['earth', 'mars']
t0='Jul 31, 2023, 00:00 UTC'
tf='Aug 07, 2023, 00:00 UTC'
cb=planetary_data.earth

sc_data={
    'mass':5000, #kg
    'reflectance':0.2,
    'area':12 #m²
}



plot_orbits(sats,bodies,t0,tf, perturbations, central_body=cb, sc_data=sc_data, steps=4000, animate=True,show=True)