import math

import planetary_data as pd
import numpy as np

def lambert_universal(r0, r, dt, mu=pd.sun['mu'], prograde=True, max_iter=100, tol=1e-6):
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
        dt_new = (chi**3 * c3 + A * math.sqrt(y)) / math.sqrt(mu)

        if abs(dt_new - dt) / dt <1e-8 :
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


def lambert_branches(r0, r, dt, mu, prograde=True):
    """
    Returns all candidate Lambert solutions for a given geometry.
    """

    solutions = []

    # short-way / prograde or retrograde
    try:
        sol = lambert_universal(r0, r, dt, mu=mu, prograde=prograde)
        solutions.append(sol)
    except:
        pass

    # long-way branch (flip transfer angle)
    try:
        sol = lambert_universal(r0, r, dt, mu=mu, prograde=not prograde)
        solutions.append(sol)
    except:
        pass

    return solutions

def compute_c3(v0, v_planet):
    v_inf = v0 - v_planet
    return np.dot(v_inf, v_inf)