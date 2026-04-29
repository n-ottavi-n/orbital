import numpy as np
import spiceypy as spice
import math

from scipy.optimize import minimize_scalar
from scipy.interpolate import CubicSpline


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
    # rotate relative state into Mars equatorial frame
    # --------------------------------------------------
    # MARSIAU is Mars-centred, Mars-equator aligned, inertial
    # use the epoch of closest approach
    tipm = spice.pxform("ECLIPJ2000", f"IAU_{body.upper()}", epoch)

    r_eq = tipm @ r  # r in Mars equatorial frame
    v_eq = tipm @ v  # v in Mars equatorial frame

    # --------------------------------------------------
    # orbital elements  (replace r,v with r_eq, v_eq below)
    # --------------------------------------------------
    rmag = np.linalg.norm(r_eq)
    vmag = np.linalg.norm(v_eq)

    h = np.cross(r_eq, v_eq)
    hmag = np.linalg.norm(h)

    k = np.array([0.0, 0.0, 1.0])

    n = np.cross(k, h)
    nmag = np.linalg.norm(n)

    #e_vec = np.cross(v, h)/mu_body - r/rmag
    e_vec = np.cross(v_eq, h) / mu_body - r_eq / rmag
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
            math.acos(np.clip(np.dot(e_vec, r_eq)/(e*rmag), -1, 1))
        )
        if np.dot(r_eq, v_eq) < 0:
            nu = 360 - nu

    if np.isfinite(a):
        rp = a * (1 - e)
        ra = a * (1 + e) if e < 1 else np.inf
    else:
        rp = hmag**2 / mu_body / (1 + e)
        ra = np.inf

    # --------------------------------------------------
    # latitude and longitude at periapsis
    # --------------------------------------------------
    # periapsis position in body-centred inertial frame
    r_peri = rp * (e_vec / np.linalg.norm(e_vec))

    # rotation matrix from inertial (ECLIPJ2000) to body-fixed frame
    # IAU_MARS rotates with the planet — gives you geographic lat/lon


    # convert to lat/lon
    # IAU_MARS rotates with Mars — but r_eq is in the instantaneous IAU_MARS frame
    # so r_peri is already body-fixed
    r_peri_norm = np.linalg.norm(r_peri)
    lat = math.degrees(math.asin(r_peri[2] / r_peri_norm))
    lon = math.degrees(math.atan2(r_peri[1], r_peri[0])) % 360

    return {
        "epoch_et": epoch,
        "arrival_utc": spice.et2utc(epoch, "C", 0),
        "range_km": rmag, # minimum |r_sc - r_Mars| in simulation

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

        "r_vec": r_eq,
        "v_vec": v_eq,
        "periapsis_latitude_deg": lat,
        "periapsis_longitude_deg": lon,
    }
