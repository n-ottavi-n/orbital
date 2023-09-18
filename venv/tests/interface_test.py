import planetary_data
from plot_orbits_interface import plot_orbits

#sats=['CHINASAT 18','YAOGAN-40 01C','BEIDOU 3 G4']
#perturbations=['j2','moon']
#bodies=['moon']
bodies=['venus','earth', 'mars']
t0='Sep 17, 2022, 00:00 UTC'
tf='Sep 18, 2023, 00:00 UTC'
sats=[]
cb=planetary_data.sun
perturbations=[]
plot_orbits(sats,bodies,t0,tf, perturbations, central_body=cb, animate=True,show=True)