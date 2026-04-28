import tools as t
from tools import ecc_anom
from tools import coes2rv
import math
import planetary_data
import numpy as np
from Propagator import Propagator
from Propagator import null_perts
import matplotlib.animation as animation
import spiceypy as spice
from datetime import datetime

cb=planetary_data.earth

spice.furnsh('spice_lunar/earth_moon_kernel.txt')


t0='Sep 17, 2023, 00:00 UTC'
tf='Oct 17, 2023, 00:00 UTC'
states, names = t.n_tle2coes("data/molniya.txt", 1, t0)

elems=[]
#propagate with j2 perts only
j2perts=null_perts()
j2perts['j2']=True

prop1=Propagator(states[0],400,start_date=t0, end_date=tf, coes=True,deg=True, perts=j2perts)
prop1.propagate()
prop1.calculate_coes()
coes1=prop1.coes
ts=prop1.ts

elems.append(coes1)

#propagate with j2 and moon perts
j2perts['moon']=True

prop2=Propagator(states[0],400,start_date=t0, end_date=tf, coes=True,deg=True, perts=j2perts)
prop2.propagate()
prop2.calculate_coes()
coes2=prop2.coes

elems.append(coes2)

labels=['j2', 'j2+moon']

t.plot_pert_coes(elems, ts, labels, days=True)