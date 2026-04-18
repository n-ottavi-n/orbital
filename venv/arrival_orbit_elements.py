import numpy as np
import spiceypy as spice
import math


def arrival_orbit_elements(sim, body, mu_body):

    # --------------------------------------------------
    # spacecraft propagated states
    # --------------------------------------------------
    #sc_states = sim.int_props[0]     # Nx6 REQUIRED
    r_sc = sim.int_props[0]
    v_sc = sim.int_vels[0]

    # target positions from propagation
    r_body = sim.int_props[-1]

    # --------------------------------------------------
    # closest approach
    # --------------------------------------------------
    rel = r_sc - r_body
    ranges = np.linalg.norm(rel, axis=1)

    i = np.argmin(ranges)
    epoch = sim.times[i]

    # --------------------------------------------------
    # target velocity from SPICE
    # --------------------------------------------------
    state_body, _ = spice.spkezr(
        body,
        epoch,
        "ECLIPJ2000",
        "NONE",
        "SUN"
    )

    v_body = state_body[3:]

    # --------------------------------------------------
    # relative state in body frame
    # --------------------------------------------------
    r = r_sc[i] - r_body[i]
    v = v_sc[i] - v_body

    # --------------------------------------------------
    # body radius
    # --------------------------------------------------
    radii = spice.bodvrd(body, "RADII", 3)[1]
    body_radius = np.mean(radii)


    # --------------------------------------------------
    # orbital elements
    # --------------------------------------------------
    rmag = np.linalg.norm(r)
    vmag = np.linalg.norm(v)

    h = np.cross(r, v)
    hmag = np.linalg.norm(h)

    k = np.array([0.0, 0.0, 1.0])

    n = np.cross(k, h)
    nmag = np.linalg.norm(n)

    e_vec = np.cross(v, h)/mu_body - r/rmag
    e = np.linalg.norm(e_vec)

    energy = vmag**2 / 2 - mu_body / rmag

    if abs(e - 1.0) < 1e-10:
        a = np.inf
    else:
        a = -mu_body / (2 * energy)

    inc = math.degrees(math.acos(np.clip(h[2]/hmag, -1, 1)))

    if nmag < 1e-12:
        raan = 0.0
    else:
        raan = math.degrees(math.atan2(n[1], n[0])) % 360

    if nmag < 1e-12 or e < 1e-12:
        argp = 0.0
    else:
        argp = math.degrees(
            math.acos(np.clip(np.dot(n, e_vec)/(nmag*e), -1, 1))
        )
        if e_vec[2] < 0:
            argp = 360 - argp

    if e < 1e-12:
        nu = 0.0
    else:
        nu = math.degrees(
            math.acos(np.clip(np.dot(e_vec, r)/(e*rmag), -1, 1))
        )
        if np.dot(r, v) < 0:
            nu = 360 - nu

    if np.isfinite(a):
        rp = a * (1 - e)
        ra = a * (1 + e) if e < 1 else np.inf
    else:
        rp = hmag**2 / mu_body / (1 + e)
        ra = np.inf

    return {
        "epoch_et": epoch,
        "range_km": rmag,

        "semi_major_axis_km": a,
        "eccentricity": e,
        "inclination_deg": inc,
        "raan_deg": raan,
        "arg_periapsis_deg": argp,
        "true_anomaly_deg": nu,

        "periapsis_km": rp,
        "apoapsis_km": ra,

        "periapsis_altitude_km": abs(rp - body_radius),
        "apoapsis_altitude_km": ra - body_radius if np.isfinite(ra) else np.inf,

        "vinf_km_s": math.sqrt(2*energy) if energy > 0 else 0.0,

        "r_vec": r,
        "v_vec": v
    }


'''import numpy as np
import spiceypy as spice
import math


def arrival_orbit_elements(sim, body, mu_body):
    """
    Extract spacecraft state at closest approach and compute
    classical orbital elements in destination-centered frame.

    Parameters
    ----------
    sim : solved lambert_interface object
        Must contain:
            sim.int_props[0]   spacecraft positions
            sim.int_props[-1]  target positions
            sim.times          epochs

    body : str
        Target body name (mars, venus, jupiter...)

    mu_body : float
        GM of target body (km^3/s^2)

    Returns
    -------
    dict
    """

    # --------------------------------------------------
    # 1. Find closest approach index
    # --------------------------------------------------
    r_sc = sim.int_props[0]
    r_body = sim.int_props[-1]

    rel = r_sc - r_body
    ranges = np.linalg.norm(rel, axis=1)

    i = np.argmin(ranges)

    epoch = sim.times[i]

    # --------------------------------------------------
    # 2. Get velocities from SPICE
    # --------------------------------------------------
    state_body, _ = spice.spkezr(body, epoch, "ECLIPJ2000", "NONE", "SUN")

    v_body = state_body[3:]

    # spacecraft velocity from finite difference
    if i == 0:
        v_sc = (r_sc[i+1] - r_sc[i]) / (sim.times[i+1] - sim.times[i])
    elif i == len(r_sc) - 1:
        v_sc = (r_sc[i] - r_sc[i-1]) / (sim.times[i] - sim.times[i-1])
    else:
        dt = sim.times[i+1] - sim.times[i-1]
        v_sc = (r_sc[i+1] - r_sc[i-1]) / dt

    # --------------------------------------------------
    # 3. Relative state in body frame
    # --------------------------------------------------
    r = r_sc[i] - r_body[i]
    v = v_sc - v_body

    # --------------------------------------------------
    # 4. Classical orbital elements
    # --------------------------------------------------
    rmag = np.linalg.norm(r)
    vmag = np.linalg.norm(v)

    h = np.cross(r, v)
    hmag = np.linalg.norm(h)

    k = np.array([0.0, 0.0, 1.0])

    n = np.cross(k, h)
    nmag = np.linalg.norm(n)

    e_vec = (np.cross(v, h) / mu_body) - (r / rmag)
    e = np.linalg.norm(e_vec)

    energy = vmag**2 / 2 - mu_body / rmag

    if abs(e - 1.0) < 1e-8:
        a = np.inf
    else:
        a = -mu_body / (2 * energy)

    i_rad = math.acos(h[2] / hmag)
    inc = math.degrees(i_rad)

    # RAAN
    if nmag < 1e-12:
        raan = 0.0
    else:
        raan = math.degrees(math.atan2(n[1], n[0])) % 360

    # argument of periapsis
    if nmag < 1e-12 or e < 1e-12:
        argp = 0.0
    else:
        argp = math.degrees(
            math.acos(np.dot(n, e_vec) / (nmag * e))
        )
        if e_vec[2] < 0:
            argp = 360 - argp

    # true anomaly
    if e < 1e-12:
        nu = 0.0
    else:
        nu = math.degrees(
            math.acos(np.dot(e_vec, r) / (e * rmag))
        )
        if np.dot(r, v) < 0:
            nu = 360 - nu

    # periapsis / apoapsis
    rp = a * (1 - e) if np.isfinite(a) else hmag**2 / mu_body / (1 + e)
    ra = a * (1 + e) if (np.isfinite(a) and e < 1) else np.inf

    # mean body radius
    radii = spice.bodvrd(body, "RADII", 3)[1]
    body_radius = np.mean(radii)

    periapsis_altitude = rp - body_radius

    if np.isfinite(ra):
        apoapsis_altitude = ra - body_radius
    else:
        apoapsis_altitude = np.inf

    return {
        "epoch_et": epoch,
        "range_km": rmag,
        "semi_major_axis_km": a,
        "eccentricity": e,
        "inclination_deg": inc,
        "raan_deg": raan,
        "arg_periapsis_deg": argp,
        "true_anomaly_deg": nu,
        "periapsis_km": rp,
        "apoapsis_km": ra,
        "body_radius_km": body_radius,
        "periapsis_altitude_km": periapsis_altitude,
        "apoapsis_altitude_km": apoapsis_altitude,
        "vinf_km_s": math.sqrt(max(0, 2 * energy)) if energy > 0 else 0.0,
        "r_vec": r,
        "v_vec": v
    }'''