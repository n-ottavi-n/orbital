from lambert_interface import lambert_interface
import planetary_data
import spiceypy as spice

spice.kclear()
spice.furnsh('../spice_solar_system/solar_system_kernel.txt')

start_date='Oct 26, 2029, 04:00 UTC'
arrival_date='Apr 04, 2030, 00:00 UTC'
end_date='Apr 12, 2030, 00:00 UTC'

steps=10000

perturbations=[planetary_data.earth, planetary_data.venus, planetary_data.jupiter, planetary_data.moon]

origin='EARTH'
dest='VENUS'


int=lambert_interface(origin,dest,start_date,arrival_date,end_date,steps,perturbations)
int.solve()
int.plot(show=True, animate=False)
int.deltaV()