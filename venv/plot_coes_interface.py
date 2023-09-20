import numpy as np
import tools as t
from Propagator import Propagator
import planetary_data as pd
import spiceypy as spice


def plot_coes():
    sc_data = {
        'mass': 100,  # kg
        'reflectance': 0.2,
        'area': 10  # m²
    }


    cb = pd.earth

    spice.furnsh('../spice_lunar/earth_moon_kernel.txt')

    t0 = 'Sep 17, 2023, 00:00 UTC'
    tf = 'Dec 17, 2023, 00:00 UTC'
    states, names = t.n_tle2coes("../data/some_geo.txt", 1, t0)

    elems = []
    # propagate with j2 and moon
    perturbations = ['j2',[pd.moon]]

    prop1 = Propagator(states[0], 400,spacecraft_data=sc_data, start_date=t0, end_date=tf, coes=True,cb=cb, deg=True, perts=perturbations)
    prop1.propagate()
    prop1.calculate_coes()
    coes1 = prop1.coes
    ts = prop1.ts

    elems.append(coes1)

    # propagate srp perts
    perturbations = ['j2','srp',[pd.moon]]

    prop2 = Propagator(states[0], 400,spacecraft_data=sc_data, start_date=t0, end_date=tf, coes=True,cb=cb, deg=True, perts=perturbations)
    prop2.propagate()
    prop2.calculate_coes()
    coes2 = prop2.coes

    elems.append(coes2)

    labels = ['base','+srp']

    t.plot_pert_coes(elems, ts, labels, days=True)

