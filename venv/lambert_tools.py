import math

import planetary_data as pd
import numpy as np


def lambert_uv(r0,r, delta_t, psi_0,psi_u,psi_l,tm=1, tol=1e-6, mu=pd.sun['mu']):
    print("solving lambert...")
    solved=False
    r0_norm=np.linalg.norm(r0)
    r_norm=np.linalg.norm(r)
    gamma=np.inner(r0,r)/(r_norm*r0_norm)
    beta=tm*math.sqrt(1-(gamma**2))
    A=tm*math.sqrt(r_norm*r0_norm*(1+gamma))
    if A==0:
        print("error A=0")
        return 0
    psi = psi_0
    for i in range(100):
        c2=C2(psi)
        c3=C3(psi)
        B=r0_norm+r_norm+((1/math.sqrt(c2))*(A*(psi*c3-1)))
        if A>0 and B<0:
            psi_l+=math.pi
            B*=-1
        chi=math.sqrt(B/c2)
        delta_tilde=(1/math.sqrt(mu))*(((chi**3)*c3)+(A*math.sqrt(B)))
        if abs(delta_t-delta_tilde)<tol:
            solved=True
            break
        if delta_tilde<=delta_t:
            psi_l=psi
        else:
            psi_u=psi
        psi=(psi_u+psi_l)/2

    if not solved:
        print("ERROR no convergence")
        return 0

    F=1-(B/r0_norm)
    G=A*math.sqrt(B/mu)
    G_dot=1-(B/r_norm)
    x,y,z=r
    x0,y0,z0=r0
    vx0 = (1 / G) * (x - (x0 * F))
    vy0 = (1 / G) * (y - (y0 * F))
    vz0 = (1 / G) * (z - (z0 * F))
    vx = (1 / G) * ((G_dot * x) - x0)
    vy = (1 / G) * ((G_dot * y) - y0)
    vz = (1 / G) * ((G_dot * z) - z0)
    return [vx0,vy0,vz0,vx,vy,vz]


def C2(psi):
    return (1-math.cos(math.sqrt(psi)))/psi

def C3(psi):
    return (math.sqrt(psi)-math.sin(math.sqrt(psi)))/math.sqrt(psi**3)

