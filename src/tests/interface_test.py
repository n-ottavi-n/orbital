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

t0='Jan 03, 2024, 11:20 UTC'
tf='Jan 01, 2029, 00:00 UTC'

t0_et = spice.str2et(t0)
sv, _lighttimes=spice.spkezr("EUROPA", t0_et, 'J2000', 'NONE', "JUPITER")
elems=spice.oscelt(sv, t0_et, planetary_data.jupiter['mu'])

print("perifocal distance: ", elems[0])
print("eccentricity: ", elems[1])
print("inclination: ", elems[2]*spice.dpr())
print("longitude of ascending node: ", elems[3]*spice.dpr())
print("argument of periapsis: ", elems[4]*spice.dpr())
print("mean anomaly: ", elems[5]*spice.dpr())

res_orbit = t.get_resonant_orbit("EUROPA", [1, 5], t0, planetary_data.jupiter)
print("resonant orbit: ", res_orbit)

sats=[['europa_probe', res_orbit[0], res_orbit[1], res_orbit[2], 0.5, res_orbit[4], res_orbit[5]]]

bodies=['europa']

cb=planetary_data.jupiter
perturbations=[[planetary_data.io, planetary_data.europa, planetary_data.ganymede, planetary_data.callisto]] # perturbations can be given as strings of body names or as planetary_data objects in a list, for example: ['j2', [pd.mars, pd.venus], 'srp'] or ['j2', 'mars', 'venus', 'srp'] both are valid
sc_data={
    'mass':5000, #kg
    'reflectance':0.2,
    'area':12 #m²
}

plot_orbits(sats, bodies, t0, tf, perturbations, central_body=cb, sc_data=sc_data, steps=100000, from_file=False, animate=False, show=True, frame_step=50, plot_coes=False)