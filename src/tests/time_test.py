import spiceypy as spice
import datetime


t0=datetime.datetime.now()

yr=t0.year%100
day_of_year = t0.timetuple().tm_yday  # returns 1 for January 1st

dt = datetime.timedelta(hours=t0.hour, minutes=t0.minute, seconds=t0.second)
secs_per_day = 24*60*60    # hours * mins * secs
frac_day=dt.total_seconds()/secs_per_day

res=yr*1000+day_of_year+frac_day
