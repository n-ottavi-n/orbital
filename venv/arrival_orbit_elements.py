import numpy as np
import spiceypy as spice
import math

from scipy.optimize import minimize_scalar
from scipy.interpolate import CubicSpline


import numpy as np
import spiceypy as spice
import math


def arrival_orbit_elements(sim, body, mu_body):

    # --------------------------------------------------
    # spacecraft propagated states
    # --------------------------------------------------
    r_sc = sim.int_props[0]
    v_sc = sim.int_vels[0]

    # --------------------------------------------------
    # closest approach using propagated body position
    # --------------------------------------------------
    r_body_prop = sim.int_props[-1]                  # CHANGED: renamed from r_body
    rel = r_sc - r_body_prop
    ranges = np.linalg.norm(rel, axis=1)

    i = np.argmin(ranges)
    epoch = sim.times[i]

    # --------------------------------------------------
    # body state from SPICE (position AND velocity)      # CHANGED: was velocity-only
    # --------------------------------------------------
    state_body, _ = spice.spkezr(                       # CHANGED: now also gets position
        body,
        epoch,
        "ECLIPJ2000",
        "NONE",
        "SUN"
    )

    r_body_spice = state_body[:3]                        # CHANGED: new — SPICE position
    v_body = state_body[3:]

    # --------------------------------------------------
    # relative state in body frame
    # --------------------------------------------------
    r = r_sc[i] - r_body_spice                          # CHANGED: was r_sc[i] - r_body[i]
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

        "periapsis_altitude_km": rp - body_radius,
        "apoapsis_altitude_km": ra - body_radius if np.isfinite(ra) else np.inf,

        "vinf_km_s": math.sqrt(2*energy) if energy > 0 else 0.0,

        "r_vec": r,
        "v_vec": v
    }
