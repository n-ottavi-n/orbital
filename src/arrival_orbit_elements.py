import numpy as np
import spiceypy as spice
import math
import planetary_data as pd


def arrival_orbit_elements(sim, body, mu_body, mu_central=1.32712440018e11):

    # --------------------------------------------------
    # spacecraft propagated states
    # --------------------------------------------------
    r_sc  = sim.int_props[0]
    v_sc  = sim.int_vels[0]
    times = np.array(sim.times)

    # --------------------------------------------------
    # body states from SPICE at each timestep
    # --------------------------------------------------
    r_body = np.array([
        spice.spkpos(body, float(t), "ECLIPJ2000", "NONE", "SUN")[0]
        for t in times
    ])
    v_body = np.array([
        spice.spkezr(body, float(t), "ECLIPJ2000", "NONE", "SUN")[0][3:]
        for t in times
    ])

    # --------------------------------------------------
    # SOI radius — computed dynamically
    # --------------------------------------------------
    state_body0, _ = spice.spkezr(body, float(times[0]), "ECLIPJ2000", "NONE", "SUN")
    r_body_sun     = np.linalg.norm(state_body0[:3])
    v_body_sun     = np.linalg.norm(state_body0[3:])
    a_body         = 1.0 / (2.0 / r_body_sun - v_body_sun**2 / mu_central)
    r_soi          = a_body * (mu_body / mu_central) ** (2.0 / 5.0)

    #print(f"SOI radius: {r_soi:.0f} km")

    # --------------------------------------------------
    # relative distance spacecraft — body
    # --------------------------------------------------
    rel    = r_sc - r_body
    ranges = np.linalg.norm(rel, axis=1)

    # --------------------------------------------------
    # find SOI entry point
    # spacecraft may or may not enter SOI depending on trajectory
    # --------------------------------------------------
    soi_entries = np.where(ranges < r_soi)[0]

    if len(soi_entries) > 0:
        # use first SOI entry — heliocentric propagation still valid here
        i     = soi_entries[0]
        mode  = "SOI entry"
    else:
        # spacecraft never enters SOI — use closest approach
        i     = np.argmin(ranges)
        mode  = "closest approach (outside SOI)"

    epoch = float(times[i])
    #print(f"Using {mode} at index {i}, range = {ranges[i]:.0f} km")

    # --------------------------------------------------
    # body state from SPICE at selected epoch
    # --------------------------------------------------
    state_body_ep, _ = spice.spkezr(body, epoch, "ECLIPJ2000", "NONE", "SUN")
    r_body_ep        = state_body_ep[:3]
    v_body_ep        = state_body_ep[3:]

    # --------------------------------------------------
    # relative state at selected epoch
    # --------------------------------------------------
    r = r_sc[i] - r_body_ep
    v = v_sc[i] - v_body_ep

    # --------------------------------------------------
    # body radius
    # --------------------------------------------------
    body_radius = pd.get_body_radius(body)

    # --------------------------------------------------
    # rotate into body equatorial frame
    # --------------------------------------------------
    try:
        tipm = spice.pxform("ECLIPJ2000", f"IAU_{body.upper()}", epoch)
        r_eq = tipm @ r
        v_eq = tipm @ v
    except:
        # no IAU frame available — use ecliptic frame
        r_eq = r
        v_eq = v

    # --------------------------------------------------
    # orbital elements in equatorial frame
    # --------------------------------------------------
    rmag = np.linalg.norm(r_eq)
    vmag = np.linalg.norm(v_eq)

    h     = np.cross(r_eq, v_eq)
    hmag  = np.linalg.norm(h)

    k     = np.array([0.0, 0.0, 1.0])
    n     = np.cross(k, h)
    nmag  = np.linalg.norm(n)

    e_vec = np.cross(v_eq, h) / mu_body - r_eq / rmag
    e     = np.linalg.norm(e_vec)

    energy = vmag**2 / 2 - mu_body / rmag

    if abs(e - 1.0) < 1e-10:
        a = np.inf
    else:
        a = -mu_body / (2 * energy)

    inc = math.degrees(math.acos(np.clip(h[2] / hmag, -1, 1)))

    if nmag < 1e-12:
        raan = 0.0
    else:
        raan = math.degrees(math.atan2(n[1], n[0])) % 360

    if nmag < 1e-12 or e < 1e-12:
        argp = 0.0
    else:
        argp = math.degrees(
            math.acos(np.clip(np.dot(n, e_vec) / (nmag * e), -1, 1))
        )
        if e_vec[2] < 0:
            argp = 360 - argp

    if e < 1e-12:
        nu = 0.0
    else:
        nu = math.degrees(
            math.acos(np.clip(np.dot(e_vec, r_eq) / (e * rmag), -1, 1))
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
    # v∞ and energy
    # --------------------------------------------------
    vinf   = math.sqrt(max(2 * energy, 0)) if energy > 0 else 0.0

    # --------------------------------------------------
    # latitude and longitude at periapsis
    # r_peri is in equatorial frame = already body-fixed
    # --------------------------------------------------
    if e > 1e-8:
        r_peri = rp * (e_vec / np.linalg.norm(e_vec))
    else:
        r_peri = r_eq / rmag * rp

    r_peri_norm = np.linalg.norm(r_peri)
    lat = math.degrees(math.asin(np.clip(r_peri[2] / r_peri_norm, -1, 1)))
    lon = math.degrees(math.atan2(r_peri[1], r_peri[0])) % 360

    # time from SOI entry to periapsis — hyperbolic TOF equation
    # at SOI entry: true anomaly = nu (already computed)
    # at periapsis: true anomaly = 0

    nu_soi   = math.radians(nu)   # true anomaly at SOI entry

    # hyperbolic anomaly at SOI entry
    F_soi    = 2.0 * math.atanh(
        math.sqrt((e - 1.0) / (e + 1.0)) * math.tan(nu_soi / 2.0)
    )

    # time from periapsis to SOI entry (then negate for SOI entry to periapsis)
    t_to_peri = math.sqrt((-a)**3 / mu_body) * (e * math.sinh(F_soi) - F_soi)

    # periapsis epoch
    periapsis_et  = epoch - t_to_peri   # subtract because spacecraft is inbound
    periapsis_utc = spice.et2utc(periapsis_et, "C", 0)

    return {
        "epoch_et":              epoch,
        "arrival_utc":           spice.et2utc(epoch, "C", 0),
        "periapsis_et":          periapsis_et,
        "periapsis_utc":         periapsis_utc,
        "t_soi_to_periapsis_hours": t_to_peri,
        "range_km":              rmag,
        "mode":                  mode,              # SOI entry or closest approach

        "semi_major_axis_km":    a,
        "eccentricity":          e,
        "inclination_deg":       inc,
        "raan_deg":              raan,
        "arg_periapsis_deg":     argp,
        "true_anomaly_deg":      nu,

        "periapsis_km":          rp,
        "apoapsis_km":           ra,

        "periapsis_altitude_km": rp - body_radius,
        "apoapsis_altitude_km":  ra - body_radius if np.isfinite(ra) else np.inf,

        "vinf_km_s":             vinf,

        "r_vec":                 r_eq,
        "v_vec":                 v_eq,

        "periapsis_latitude_deg":  lat,
        "periapsis_longitude_deg": lon,

        "soi_entry_idx": i,    # index into sim.times where spacecraft enters SOI
        "r_soi_km":              r_soi,             # for reference
    }