import planetary_data
from plot_orbits_interface import plot_orbits
import spiceypy as spice
import tools as t
'''
#sats=['CHINASAT 18','YAOGAN-40 01C','BEIDOU 3 G4']
sats=['STARLINK-30378']
#perturbations=[[planetary_data.mars]]
perturbations=['j2', [planetary_data.moon]]

bodies=['earth']
t0='Jul 31, 2023, 00:00 UTC'
tf='Jul 31, 2024, 00:00 UTC'
cb=planetary_data.earth
'''

t.load_solar_system_kernels()

sv, _lighttimes=spice.spkezr("EUROPA", spice.str2et('Jan 01, 2024, 00:00 UTC'), 'J2000', 'NONE', "JUPITER")
elems=spice.oscelt(sv, spice.str2et('Jan 01, 2024, 00:00 UTC'), planetary_data.jupiter['mu'])
print("perifocal distance: ", elems[0])
print("eccentricity: ", elems[1])
print("inclination: ", elems[2]*spice.dpr())
print("longitude of ascending node: ", elems[3]*spice.dpr())
print("argument of periapsis: ", elems[4]*spice.dpr())
print("mean anomaly: ", elems[5]*spice.dpr())

sats=[['europa_probe', 1664000, 0.6, 25.426, 130, 225.071, 357.011]]

bodies=['europa', 'io', 'ganymede', 'callisto']
t0='Jan 01, 2024, 00:00 UTC'
tf='Jan 01, 2025, 00:00 UTC'
cb=planetary_data.jupiter
perturbations=[[planetary_data.io, planetary_data.europa, planetary_data.ganymede, planetary_data.callisto]] # perturbations can be given as strings of body names or as planetary_data objects in a list, for example: ['j2', [pd.mars, pd.venus], 'srp'] or ['j2', 'mars', 'venus', 'srp'] both are valid
sc_data={
    'mass':5000, #kg
    'reflectance':0.2,
    'area':12 #m²
}

plot_orbits(sats, bodies, t0, tf, perturbations, central_body=cb, sc_data=sc_data, steps=50000, from_file=False, animate=True, show=True, frame_step=50)