import math

import planetary_data as pd
import numpy as np

def lambert_universal(r0, r, dt, mu=pd.sun['mu'], prograde=True, max_iter=100, tol=1e-6):

    sqrt_mu = math.sqrt(mu)
    r0 = np.array(r0)
    r = np.array(r)

    r0_norm = np.linalg.norm(r0)
    r_norm = np.linalg.norm(r)

    # --- Geometry ---
    dot = np.dot(r0, r) / (r0_norm * r_norm)
    dot = np.clip(dot, -1.0, 1.0)  # numerical safety
    dtheta = math.acos(dot)

    h = np.cross(r0, r)
    sign = np.sign(h[2]) if np.linalg.norm(h) > 0 else 1.0

    if prograde:
        if sign < 0:
            dtheta = 2 * math.pi - dtheta
    else:
        if sign > 0:
            dtheta = 2 * math.pi - dtheta

    # --- A parameter ---
    A = math.sin(dtheta) * math.sqrt(r0_norm * r_norm / (1 - math.cos(dtheta)))

    if abs(A) < 1e-8:
        raise RuntimeError("Trajectory undefined (A ≈ 0)")

    # --- Initial bounds for psi ---
    psi = 0.0
    psi_low = -4 * math.pi**2
    psi_up = 4 * math.pi**2

    # --- Iteration ---
    for _ in range(max_iter):
        c2 = max(C2(psi), 1e-8)
        c3 = C3(psi)

        y = r0_norm + r_norm + A * (psi * c3 - 1) / math.sqrt(c2)
        y = max(y, 1e-8)


        # Ensure physical solution
        if A > 0 and y < 0:
            psi_low += 0.1
            psi = (psi_up + psi_low) / 2
            continue

        chi = math.sqrt(y / c2)
        dt_new = (chi**3 * c3 + A * math.sqrt(y)) / sqrt_mu

        if abs(dt_new - dt) / dt <tol :
            break

        if dt_new <= dt:
            psi_low = psi
        else:
            psi_up = psi

        psi = (psi_up + psi_low) / 2
    else:
        raise RuntimeError("Lambert solver did not converge")

    # --- Lagrange coefficients ---
    f = 1 - y / r0_norm
    g = A * math.sqrt(y / mu)
    gdot = 1 - y / r_norm

    v0 = (r - f * r0) / g
    v = (gdot * r - r0) / g

    return np.concatenate((v0, v))


# --- Stumpff functions ---
def C2(psi):
    if psi > 1e-6:
        return (1 - math.cos(math.sqrt(psi))) / psi
    elif psi < -1e-6:
        return (1 - math.cosh(math.sqrt(-psi))) / psi
    else:
        return 0.5


def C3(psi):
    if psi > 1e-6:
        return (math.sqrt(psi) - math.sin(math.sqrt(psi))) / (psi**1.5)
    elif psi < -1e-6:
        return (math.sinh(math.sqrt(-psi)) - math.sqrt(-psi)) / ((-psi)**1.5)
    else:
        return 1 / 6


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

'''
def C2(psi):
    return (1-math.cos(math.sqrt(psi)))/psi

def C3(psi):
    return (math.sqrt(psi)-math.sin(math.sqrt(psi)))/math.sqrt(psi**3)

'''