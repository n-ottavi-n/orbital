import numpy as np
import spiceypy as spice
import math


def departure_orbit_elements(sim, body, mu_body, parking_alt_km=400.0,
                             argp_deg=0.0, mu_central=1.32712440018e11):
    """
    Compute departure hyperbola and parking orbit elements from a propagated simulation.

    Parameters
    ----------
    sim            : lambert_interface object (after solve())
    body           : SPICE name of departure body (e.g. "EARTH")
    mu_body        : gravitational parameter of departure body (km³/s²)
    parking_alt_km : circular parking orbit altitude (km)
    argp_deg       : argument of periapsis of departure hyperbola (deg)
                     currently fixed — reserved for future optimization
    mu_central     : gravitational parameter of central body (km³/s²)
                     defaults to Sun - used for SOI computation
    """

    # --------------------------------------------------
    # spacecraft state at departure from Lambert solution
    # --------------------------------------------------
    r_sc_helio = sim.start_r        # departure position (ECLIPJ2000, heliocentric)
    v_sc_helio = sim.v0             # corrected departure velocity
    epoch      = sim.start_et       # SOI exit epoch (ET)

    # --------------------------------------------------
    # body state from SPICE at departure epoch
    # --------------------------------------------------
    state_body, _ = spice.spkezr(body, epoch, "ECLIPJ2000", "NONE", "SUN")
    v_body_spice  = state_body[3:]
    r_body_spice  = state_body[:3]

    # --------------------------------------------------
    # hyperbolic excess velocity vector (ECLIPJ2000)
    # --------------------------------------------------
    v_inf_vec = v_sc_helio - v_body_spice
    vinf      = np.linalg.norm(v_inf_vec)
    c3        = vinf**2

    # --------------------------------------------------
    # body radius and parking orbit
    # --------------------------------------------------
    radii       = spice.bodvrd(body, "RADII", 3)[1]
    body_radius = np.mean(radii)
    r_parking   = body_radius + parking_alt_km

    # --------------------------------------------------
    # departure burn delta-V
    # --------------------------------------------------
    v_circ  = math.sqrt(mu_body / r_parking)
    v_peri  = math.sqrt(vinf**2 + 2 * mu_body / r_parking)
    delta_v = v_peri - v_circ

    # --------------------------------------------------
    # departure hyperbola elements
    # --------------------------------------------------
    e_hyp = 1 + (r_parking * vinf**2) / mu_body
    a_hyp = -mu_body / c3                            # negative (hyperbola)
    p_hyp = r_parking * (e_hyp + 1)                 # semi-latus rectum

    # --------------------------------------------------
    # SOI radius — computed dynamically from vis-viva
    # r_SOI = a_body * (mu_body / mu_central)^(2/5)
    # --------------------------------------------------
    r_mag_body = np.linalg.norm(r_body_spice)
    v_mag_body = np.linalg.norm(v_body_spice)
    a_body     = 1.0 / (2.0 / r_mag_body - v_mag_body**2 / mu_central)
    r_soi      = a_body * (mu_body / mu_central) ** (2.0 / 5.0)

    # --------------------------------------------------
    # time from periapsis to SOI crossing
    # using hyperbolic time of flight equation
    # --------------------------------------------------
    cos_nu_soi = np.clip((p_hyp / r_soi - 1.0) / e_hyp, -1.0, 1.0)
    nu_soi     = math.acos(cos_nu_soi)

    # hyperbolic anomaly at SOI
    F_soi = 2.0 * math.atanh(
        math.sqrt((e_hyp - 1.0) / (e_hyp + 1.0)) * math.tan(nu_soi / 2.0)
    )

    # time from periapsis to SOI (seconds)
    t_soi = math.sqrt((-a_hyp)**3 / mu_body) * (e_hyp * math.sinh(F_soi) - F_soi)

    # --------------------------------------------------
    # injection epoch = SOI exit time - travel time
    # --------------------------------------------------
    injection_et  = epoch - t_soi
    injection_utc = spice.et2utc(injection_et, "C", 0)

    # --------------------------------------------------
    # rotate into body equatorial frame (IAU_BODY)
    # --------------------------------------------------
    iau_frame = f"IAU_{body.upper()}"
    tipm      = spice.pxform("ECLIPJ2000", iau_frame, epoch)
    v_inf_eq  = tipm @ v_inf_vec

    # --------------------------------------------------
    # asymptote right ascension and declination
    # --------------------------------------------------
    v_inf_hat = v_inf_eq / np.linalg.norm(v_inf_eq)
    dec       = math.degrees(math.asin(np.clip(v_inf_hat[2], -1, 1)))
    ra        = math.degrees(math.atan2(v_inf_hat[1], v_inf_hat[0])) % 360
    inc_min   = abs(dec)

    # --------------------------------------------------
    # parking orbit plane — minimum inclination solution
    # project north pole onto plane perpendicular to v_inf_hat
    # --------------------------------------------------
    north      = np.array([0.0, 0.0, 1.0])
    north_perp = north - np.dot(north, v_inf_hat) * v_inf_hat

    if np.linalg.norm(north_perp) < 1e-8:
        x_axis     = np.array([1.0, 0.0, 0.0])
        north_perp = x_axis - np.dot(x_axis, v_inf_hat) * v_inf_hat

    h_hat = north_perp / np.linalg.norm(north_perp)

    # periapsis direction in orbital plane
    e_hat = np.cross(h_hat, v_inf_hat)
    e_hat /= np.linalg.norm(e_hat)

    # --------------------------------------------------
    # apply argp rotation around h_hat
    # reserved for future optimization — currently fixed at argp_deg=0
    # --------------------------------------------------
    argp_rad = math.radians(argp_deg)
    q_hat    = np.cross(h_hat, e_hat)
    e_hat_r  = math.cos(argp_rad) * e_hat + math.sin(argp_rad) * q_hat

    # --------------------------------------------------
    # periapsis position and velocity vectors (equatorial frame)
    # --------------------------------------------------
    r_peri_vec = r_parking * e_hat_r
    v_peri_vec = v_peri * np.cross(h_hat, e_hat_r)

    # --------------------------------------------------
    # parking orbit elements
    # --------------------------------------------------
    k    = np.array([0.0, 0.0, 1.0])
    n    = np.cross(k, h_hat)
    nmag = np.linalg.norm(n)

    inc_park = math.degrees(math.acos(np.clip(h_hat[2], -1, 1)))

    if nmag < 1e-12:
        raan = 0.0
    else:
        raan = math.degrees(math.atan2(n[1], n[0])) % 360

    # --------------------------------------------------
    # periapsis latitude and longitude (body-fixed)
    # --------------------------------------------------
    r_peri_norm = np.linalg.norm(r_peri_vec)
    lat = math.degrees(math.asin(np.clip(r_peri_vec[2] / r_peri_norm, -1, 1)))
    lon = math.degrees(math.atan2(r_peri_vec[1], r_peri_vec[0])) % 360

    return {
        # timing
        "departure_et":            epoch,           # SOI exit epoch (ET)
        "departure_utc":           spice.et2utc(epoch, "C", 0),
        "injection_et":            injection_et,    # parking orbit burn epoch (ET)
        "injection_utc":           injection_utc,   # human-readable injection time
        "t_soi_hours":             t_soi / 3600,    # travel time periapsis → SOI (hours)

        # SOI
        "r_soi_km":                r_soi,           # dynamically computed SOI radius

        # hyperbolic excess
        "c3_km2_s2":               c3,
        "vinf_km_s":               vinf,
        "vinf_vec_eq":             v_inf_eq,

        # asymptote geometry
        "asymptote_ra_deg":        ra,
        "asymptote_dec_deg":       dec,
        "asymptote_inc_min_deg":   inc_min,

        # departure hyperbola
        "hyperbola_e":             e_hyp,
        "hyperbola_a_km":          a_hyp,
        "delta_v_km_s":            delta_v,
        "v_circ_km_s":             v_circ,
        "v_peri_km_s":             v_peri,

        # parking orbit
        "parking_alt_km":          parking_alt_km,
        "parking_radius_km":       r_parking,
        "parking_inc_deg":         inc_park,
        "parking_raan_deg":        raan,

        # periapsis geometry
        "periapsis_latitude_deg":  lat,
        "periapsis_longitude_deg": lon,
        "r_vec":                   r_peri_vec,      # for plot_departure
        "v_vec":                   v_peri_vec,      # for plot_departure

        # reserved for future optimization
        "argp_deg":                argp_deg,
    }