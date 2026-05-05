import spiceypy as spice
import numpy as np

#source values: https://fr.mathworks.com/help/aeroblks/zonalharmonicgravitymodel.html
bodies = {
    'sun':{
        'name':"sun",
        'mass':1.989e30,
        'mu':1.32712e11,
        'radius':695700.0
    },

    'mercury':{
        'name':"mercury",
        'mass':0,
        'mu':2.2032e4,
        'radius':2439.0,
        'j2':0.00006
    },

    'venus':{
        'name':"venus",
        'mass':0,
        'mu':3.257e5,
        'radius':6052.0,
        'j2':0.000027
    },

    'earth':{
        'name':"earth",
        'mass':5.972e24,
        'mu':398600.4415,
        'radius':6378.136, #semi-major axis of reference ellipsoid
        'j2':0.0010826269,
        'j3':-0.0000025323
    },

    'moon':{
        'name':"moon",
        'mass':0.07346e24,
        'mu':4902.799,
        'radius':1738.1,
        'j2':0.0002027,
        'j3':0
    },

    'mars':{
        'name':"mars",
        'mass':0,
        'mu':4.305e4,
        'radius':3397.2,
        'j2': 0.001964
    },

    'ceres':{
        'name':"ceres",
        'mass':0,
        'mu':62,
        'radius':469.7, #mean radius
        'j2':0
    },

    'jupiter':{
        'name':"jupiter barycenter",
        'mass':0,
        'mu':1.268e8,
        'radius':71492.0,
        'j2':0.01475
    },

    'saturn':{
        'name':"saturn",
        'mass':0,
        'mu':3.794e7,
        'radius':60268.0,
        'j2':0.01645
    },

    '54404085':{
        'name':"54404085",
        'mass':0,
        'mu':1.582e-13,
        'radius':6,
        'j2':0
    },
}


def get_body(name):
    """Case-insensitive body lookup by name or designation."""
    key = str(name).lower()
    if key in bodies:
        return bodies[key]
    # try exact match for numeric designations
    if str(name) in bodies:
        return bodies[str(name)]
    raise ValueError(f"Body '{name}' not found in planetary_data")

def get_bodies(names):
    """Return list of body dicts for a list of names."""
    return [get_body(n) for n in names]

def get_body_radius(body_name):
    """
    Get mean body radius in km.
    Tries SPICE kernel pool first, falls back to planetary_data dict.
    """
    try:
        radii = spice.bodvrd(body_name, "RADII", 3)[1]
        return float(np.mean(radii))
    except:
        return get_body(body_name)['radius']
    
# convenience references for backwards compatibility
mercury   = bodies['mercury']
venuus   = bodies['venus']
sun   = bodies['sun']
earth = bodies['earth']
mars  = bodies['mars']
moon  = bodies.get('moon')
jupiter   = bodies['jupiter']
saturn   = bodies['saturn']
ceres   = bodies['ceres']