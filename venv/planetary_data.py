#source values: https://fr.mathworks.com/help/aeroblks/zonalharmonicgravitymodel.html

sun={
    'name':"sun",
    'mass':1.989e30,
    'mu':1.32712e11,
    'radius':695700.0
}

mercury={
    'name':"mercury",
    'mass':0,
    'mu':	2.2032e4,
    'radius':2439.0,
    'j2':0.00006
}

venus={
    'name':"venus",
    'mass':0,
    'mu':3.257e5,
    'radius':6052.0,
    'j2':0.000027
}


earth={
    'name':"earth",
    'mass':5.972e24,
    'mu':398600.4415,
    'radius':6378.136, #semi-major axis of reference ellipsoid
    'j2':0.0010826269,
    'j3':-0.0000025323
}

moon={
    'name':"moon",
    'mass':0.07346e24,
    'mu':4902.799,
    'radius':1738.1,
    'j2':0.0002027,
    'j3':0
}

mars={
    'name':"mars",
    'mass':0,
    'mu':4.305e4,
    'radius':3397.2,
    'j2': 0.001964
}

jupiter={
    'name':"jupiter",
    'mass':0,
    'mu':1.268e8,
    'radius':71492.0,
    'j2':0.01475
}

saturn={
    'name':"saturn",
    'mass':0,
    'mu':3.794e7,
    'radius':60268.0,
    'j2':0.01645
}