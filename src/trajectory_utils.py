import math

import planetary_data as pd
import numpy as np
import spiceypy as spice

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

def hohmann_phase_angle(a1, a2, mu):
    t_H = np.pi * np.sqrt((a1 + a2)**3 / (8 * mu))
    n2 = np.sqrt(mu / a2**3)
    phi = np.pi - n2 * t_H
    return np.degrees(phi)

def compute_arrival_inclination(v_inf_arr_vec, target_body, epoch, mu_body=None):
    """
    Compute arrival inclination relative to target body equatorial plane
    from the hyperbolic excess velocity vector.
    The asymptote direction of the hyperbola determines the orbital plane,
    which sets the inclination independently of periapsis altitude.
    """
    # unit vector of incoming asymptote (reversed v_inf points INTO the SOI)
    v_inf_hat = -v_inf_arr_vec / np.linalg.norm(v_inf_arr_vec)

    # get body pole direction in ECLIPJ2000
    # the z-axis of IAU_BODY frame is the north pole
    tipm = spice.pxform(f"IAU_{target_body.upper()}", "ECLIPJ2000", epoch)
    pole = tipm @ np.array([0.0, 0.0, 1.0])
    pole_hat = pole / np.linalg.norm(pole)

    # inclination = angle between orbital plane normal and pole
    # orbital plane normal is perpendicular to v_inf (and to some periapsis vector)
    # but the asymptote alone constrains inclination to within a range
    # the angle between v_inf and the equatorial plane gives the declination of the asymptote
    # inc = 90° - |dec| for equatorial approach, but more precisely:
    dec = math.degrees(math.asin(np.clip(np.dot(v_inf_hat, pole_hat), -1, 1)))

    # inclination is bounded by: |dec| <= inc <= 180° - |dec|
    # minimum inclination achievable = |declination of asymptote|
    inc_min = abs(dec)

    return inc_min