import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import spiceypy as spice


def plot_departure(elements: dict, body_name: str = "origin", body_radius_km: float = None,
                   mu_body: float = None, frame: str = "IAU"):

    r_vec          = np.array(elements["r_vec"])
    v_vec          = np.array(elements["v_vec"])
    e_hyp          = elements["hyperbola_e"]
    a_hyp          = elements["hyperbola_a_km"]
    rp             = elements["parking_radius_km"]
    vinf           = elements["vinf_km_s"]
    delta_v        = elements["delta_v_km_s"]
    inc            = elements["parking_inc_deg"]
    parking_alt    = elements["parking_alt_km"]
    v_inf_eq       = np.array(elements["vinf_vec_eq"])
    injection_utc  = elements["injection_utc"]
    t_soi_hours    = elements["t_soi_hours"]
    r_soi_km       = elements["r_soi_km"]

    if body_radius_km is None:
        body_radius_km = rp - parking_alt

    if mu_body is None:
        # recover from vis-viva: a = -mu/c3
        c3     = elements["c3_km2_s2"]
        mu_body = -a_hyp * c3

    # --------------------------------------------------
    # frame selection
    # --------------------------------------------------
    if frame == "J2000":
        # rotate from IAU_EARTH back to J2000 (ECLIPJ2000)
        epoch  = elements["departure_et"]
        tipm   = spice.pxform(f"IAU_{body_name.upper()}", "ECLIPJ2000", epoch)
        r_vec  = tipm @ r_vec
        v_vec  = tipm @ v_vec
        # v_inf is already in IAU_EARTH in elements — rotate it too
        v_inf_eq = tipm @ v_inf_eq
        frame_label = "ECLIPJ2000"
    else:
        frame_label = "IAU_EARTH"

    # --------------------------------------------------
    # perifocal frame from r_vec and v_vec
    # --------------------------------------------------
    h_vec = np.cross(r_vec, v_vec)
    hmag  = np.linalg.norm(h_vec)
    h_hat = h_vec / hmag
    rmag  = np.linalg.norm(r_vec)

    e_vec = np.cross(v_vec, h_vec) / mu_body - r_vec / rmag
    e_hat = e_vec / np.linalg.norm(e_vec) if e_hyp > 1e-8 else r_vec / rmag
    q_hat = np.cross(h_hat, e_hat)

    # --------------------------------------------------
    # hyperbola sweep — departure branch only (nu: 0 → nu_inf)
    # --------------------------------------------------
    nu_inf = np.arccos(-1.0 / e_hyp)
    nu_lim = nu_inf - np.radians(5)
    nu_arr = np.linspace(0, nu_lim, 400)   # departure: 0 (periapsis) outward

    p_orb = abs(a_hyp) * (e_hyp**2 - 1)   # semi-latus rectum for hyperbola

    def radius_at(nu):
        return p_orb / (1 + e_hyp * np.cos(nu))

    pts = np.array([
        radius_at(nu) * (np.cos(nu) * e_hat + np.sin(nu) * q_hat)
        for nu in nu_arr
    ])

    # --------------------------------------------------
    # asymptote ray (outgoing direction)
    # --------------------------------------------------
    nu_inf    = np.arccos(-1.0 / e_hyp)
    v_inf_hat = -(np.cos(nu_inf) * e_hat + np.sin(nu_inf) * q_hat)

    # --------------------------------------------------
    # parking orbit circle
    # --------------------------------------------------
    th         = np.linspace(0, 2 * np.pi, 300)
    park_pts   = np.array([
        rp * (np.cos(t) * e_hat + np.sin(t) * q_hat)
        for t in th
    ])

    # --------------------------------------------------
    # scale
    # --------------------------------------------------
    plot_radius = rp * 1.8

    # --------------------------------------------------
    # figure
    # --------------------------------------------------
    fig = plt.figure(figsize=(9, 8))
    ax  = fig.add_subplot(111, projection="3d")

    # --- body sphere ---
    u, v = np.mgrid[0:2*np.pi:60j, 0:np.pi:30j]
    R  = body_radius_km
    bx = R * np.cos(u) * np.sin(v)
    by = R * np.sin(u) * np.sin(v)
    bz = R * np.cos(v)
    ax.plot_surface(bx, by, bz, color="steelblue", alpha=0.5, linewidth=0, zorder=1)
    th_eq = np.linspace(0, 2 * np.pi, 200)
    ax.plot(R * np.cos(th_eq), R * np.sin(th_eq), np.zeros_like(th_eq),
            color="steelblue", lw=0.7, alpha=0.8)

    # --- parking orbit ---
    ax.plot(park_pts[:, 0], park_pts[:, 1], park_pts[:, 2],
            color="seagreen", lw=1.2, linestyle="--",
            label=f"parking orbit  {parking_alt:.0f} km  i={inc:.1f}°", zorder=2)

    # --- departure hyperbola ---
    ax.plot(pts[:, 0], pts[:, 1], pts[:, 2],
            color="royalblue", lw=1.8, label="departure trajectory", zorder=3)

    # --- periapsis marker ---
    peri_pt = rp * e_hat
    ax.scatter(*peri_pt, color="tomato", s=60, zorder=5,
               label=f"periapsis  alt={parking_alt:.0f} km")

    # --- velocity vector at periapsis ---
    v_peri_dir = q_hat * np.sign(np.dot(v_vec, q_hat))
    vscale     = plot_radius * 0.18 / np.linalg.norm(v_vec)
    ax.quiver(*peri_pt, *(v_peri_dir * np.linalg.norm(v_vec) * vscale),
              color="orange", linewidth=2, arrow_length_ratio=0.18,
              label=f"v periapsis  {elements['v_peri_km_s']:.2f} km/s")

    # center of hyperbola is offset from focus (Earth) along -e_hat by a*e
    center     = -abs(a_hyp) * e_hyp * e_hat   # hyperbola center in body frame

    # asymptote passes through center, not through Earth
    # both ray endpoints use the same v_inf_hat
    ray_start = pts[-1]
    ray_end   = ray_start + rp * 0.5 * v_inf_hat

    ax.plot([ray_start[0], ray_end[0]],
            [ray_start[1], ray_end[1]],
            [ray_start[2], ray_end[2]],
            "--", color="darkorange", lw=1.2, alpha=0.8,
            label=f"v∞ = {vinf:.2f} km/s")

    # --- body-frame axes ---
    ax_len = body_radius_km * 1.6
    for vec, col, lbl in zip(
        [np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1])],
        ["red", "green", "royalblue"],
        ["x", "y", "z"]
    ):
        ax.quiver(0, 0, 0, *(vec * ax_len), color=col, linewidth=0.8,
                  arrow_length_ratio=0.25, alpha=0.6)
        ax.text(*(vec * ax_len * 1.15), lbl, color=col, fontsize=8, alpha=0.8)

    # --- formatting ---
    lim = plot_radius
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_zlim(-lim, lim)
    ax.set_box_aspect([1, 1, 1])

    ax.set_xlabel("x (km)", labelpad=6, fontsize=9)
    ax.set_ylabel("y (km)", labelpad=6, fontsize=9)
    ax.set_zlabel("z (km)", labelpad=6, fontsize=9)

    ax.set_facecolor("white")
    fig.patch.set_facecolor("white")
    ax.tick_params(colors="black", labelsize=7)
    for pane in [ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane]:
        pane.fill = False
        pane.set_edgecolor("#cccccc")
    ax.grid(True, color="#dddddd", linewidth=0.4)

    # --- title ---
    title = (
        f"{body_name} departure  —  {frame_label}\n"
        f"C3 = {elements['c3_km2_s2']:.2f} km²/s²   "
        f"v∞ = {vinf:.3f} km/s   "
        f"ΔV = {delta_v:.3f} km/s\n"
        f"SOI transit: {t_soi_hours:.1f} h   "
        f"i = {elements['parking_inc_deg']:.2f}°\n"
        f"injection: {injection_utc}   "
    )
    ax.set_title(title, color="black", fontsize=10, pad=10)
    ax.legend(loc="upper left", fontsize=8,
              facecolor="white", edgecolor="#cccccc", labelcolor="black")

    plt.tight_layout()

    return fig, ax