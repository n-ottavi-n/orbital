from lambert_interface import lambert_interface

start_date='Jul 30, 2020, 00:00 UTC'
arrival_date='Feb 18, 2021, 00:00 UTC'
end_date='Jul 30, 2021, 00:00 UTC'

steps=10000

int=lambert_interface('earth','mars',start_date,arrival_date,end_date,steps)
int.solve()
int.plot(show=True, animate=False)
int.deltaV()