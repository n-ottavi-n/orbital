from lambert_interface import lambert_interface
import planetary_data
import spiceypy as spice

spice.kclear()
spice.furnsh('../spice_solar_system/solar_system_kernel.txt')

start_date='Dec 24, 2028, 21:00 UTC'
arrival_date='jul 02, 2031, 00:00 UTC'
end_date='Oct 12, 2031, 00:00 UTC'

steps=10000

perturbations=[planetary_data.earth, planetary_data.jupiter, planetary_data.moon]

origin='EARTH'
dest='JUPITER BARYCENTER'


int=lambert_interface(origin,dest,start_date,arrival_date,end_date,steps,perturbations)
int.solve()
int.plot(show=True, animate=False)
int.deltaV()