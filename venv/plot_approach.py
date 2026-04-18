import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def plot_approach(elements: dict, body_name: str = "target", body_radius_km: float = None):
    """
    Plot an approach trajectory in the body-centred frame.

    Parameters
    ----------
    elements : dict
        Output dict from arrival_orbit_elements().
    body_name : str
        Label for the central body (e.g. "Mars").
    body_radius_km : float, optional
        Override body radius. Falls back to elements["periapsis_km"] - elements["periapsis_altitude_km"].
    """

    # ------------------------------------------------------------------
    # unpack elements
    # ------------------------------------------------------------------
    r_vec = np.array(elements["r_vec"])          # position at closest approach (body frame)
    v_vec = np.array(elements["v_vec"])          # velocity at closest approach (body frame)
    e     = elements["eccentricity"]
    a     = elements["semi_major_axis_km"]
    rp    = elements["periapsis_km"]
    inc   = np.radians(elements["inclination_deg"])
    raan  = np.radians(elements["raan_deg"])
    argp  = np.radians(elements["arg_periapsis_deg"])
    nu_ca = np.radians(elements["true_anomaly_deg"])   # true anomaly at closest approach

    if body_radius_km is None:
        body_radius_km = rp - elements["periapsis_altitude_km"]

    mu = elements.get("mu_km3_s2", None)
    # recompute mu from vis-viva if not stored: mu = -a * 2 * energy, energy from v/r
    if mu is None:
        rmag = np.linalg.norm(r_vec)
        vmag = np.linalg.norm(v_vec)
        energy = vmag**2 / 2 - np.linalg.norm(r_vec) * 0   # can't recover mu cleanly without it
        # fallback: derive from h and conic
        h_vec = np.cross(r_vec, v_vec)
        hmag  = np.linalg.norm(h_vec)
        p     = hmag**2 / (mu if mu else 1)
        # best we can do: use p = a(1-e^2) → mu = hmag^2 / p
        if abs(e - 1.0) > 1e-6:
            p_conic = abs(a) * abs(1 - e**2)
        else:
            p_conic = hmag**2   # parabolic: p = h²/mu → mu = h²/p but we don't have mu yet
        mu = hmag**2 / p_conic

    # ------------------------------------------------------------------
    # perifocal frame unit vectors (ê_p toward periapsis, ê_q 90° ahead)
    # ------------------------------------------------------------------
    h_vec = np.cross(r_vec, v_vec)
    hmag  = np.linalg.norm(h_vec)
    h_hat = h_vec / hmag

    # eccentricity vector points toward periapsis
    rmag  = np.linalg.norm(r_vec)
    e_vec = np.cross(v_vec, h_vec) / mu - r_vec / rmag
    e_hat = e_vec / np.linalg.norm(e_vec) if e > 1e-8 else r_vec / rmag

    q_hat = np.cross(h_hat, e_hat)   # completes right-handed perifocal frame

    # ------------------------------------------------------------------
    # true anomaly sweep
    # ------------------------------------------------------------------
    if e < 1.0:
        # ellipse: full orbit
        nu_arr = np.linspace(0, 2 * np.pi, 500)
    elif abs(e - 1.0) < 1e-6:
        # parabola: ±170°
        nu_arr = np.linspace(-np.radians(170), np.radians(170), 500)
    else:
        # hyperbola: clip at asymptote ±5° inside
        nu_inf = np.arccos(-1.0 / e)
        nu_lim = nu_inf - np.radians(5)
        nu_arr = np.linspace(-nu_lim, nu_lim, 500)

    # ------------------------------------------------------------------
    # orbit curve in body-centred Cartesian
    # ------------------------------------------------------------------
    p_orb = abs(a) * abs(1 - e**2) if abs(e - 1.0) > 1e-6 else hmag**2 / mu

    def radius_at(nu):
        return p_orb / (1 + e * np.cos(nu))

    pts = np.array([
        radius_at(nu) * (np.cos(nu) * e_hat + np.sin(nu) * q_hat)
        for nu in nu_arr
    ])

    # ------------------------------------------------------------------
    # asymptote direction for hyperbola
    # ------------------------------------------------------------------
    asymptote = None
    if e > 1.0:
        nu_inf = np.arccos(-1.0 / e)
        # incoming asymptote direction (unit vector)
        asymptote = np.cos(-nu_inf) * e_hat + np.sin(-nu_inf) * q_hat

    # ------------------------------------------------------------------
    # scale: show a few rp out
    # ------------------------------------------------------------------
    #plot_radius = max(rp * 4, body_radius_km * 5)
    plot_radius = rp * 1.8

    # closest-approach point on the orbit curve (true anomaly = nu_ca)
    ca_pt = radius_at(nu_ca) * (np.cos(nu_ca) * e_hat + np.sin(nu_ca) * q_hat)

    # ------------------------------------------------------------------
    # figure
    # ------------------------------------------------------------------
    fig = plt.figure(figsize=(9, 8))
    ax  = fig.add_subplot(111, projection="3d")

    # --- body sphere ---
    u, v = np.mgrid[0:2*np.pi:60j, 0:np.pi:30j]
    R = body_radius_km
    bx = R * np.cos(u) * np.sin(v)
    by = R * np.sin(u) * np.sin(v)
    bz = R * np.cos(v)
    ax.plot_surface(bx, by, bz, color="steelblue", alpha=0.35, linewidth=0, zorder=1)
    # equator ring
    th = np.linspace(0, 2*np.pi, 200)
    ax.plot(R*np.cos(th), R*np.sin(th), np.zeros_like(th),
            color="steelblue", lw=0.7, alpha=0.6)

    # --- orbit curve ---
    ax.plot(pts[:, 0], pts[:, 1], pts[:, 2],
            color="gold", lw=1.8, label="approach trajectory", zorder=3)

    # --- periapsis marker ---
    peri_pt = rp * e_hat
    ax.scatter(*peri_pt, color="tomato", s=60, zorder=5, label=f"periapsis  {rp:.0f} km")

    # --- closest-approach marker ---
    ax.scatter(*ca_pt, color="lime", s=80, marker="*", zorder=6,
               label=f"closest approach  {elements['range_km']:.0f} km")

    # --- velocity vector at closest approach (scaled for visibility) ---
    vscale = plot_radius * 0.18 / np.linalg.norm(v_vec)
    ax.quiver(*ca_pt, *(v_vec * vscale),
              color="white", linewidth=1.2, arrow_length_ratio=0.18,
              label=f"v  {np.linalg.norm(v_vec):.2f} km/s")

    # --- v∞ asymptote ray for hyperbolas ---
    if asymptote is not None:
        ray_start = -asymptote * plot_radius * 0.9
        ray_end   = -asymptote * body_radius_km * 1.5
        ax.plot([ray_start[0], ray_end[0]],
                [ray_start[1], ray_end[1]],
                [ray_start[2], ray_end[2]],
                "--", color="orange", lw=1.2, alpha=0.7,
                label=f"v∞ = {elements['vinf_km_s']:.2f} km/s")

    # --- body-frame axes (small) ---
    ax_len = body_radius_km * 1.6
    for vec, col, lbl in zip(
        [np.array([1,0,0]), np.array([0,1,0]), np.array([0,0,1])],
        ["#ff6666", "#66ff66", "#6699ff"],
        ["x", "y", "z"]
    ):
        ax.quiver(0, 0, 0, *(vec*ax_len), color=col, linewidth=0.8,
                  arrow_length_ratio=0.25, alpha=0.5)
        ax.text(*(vec*ax_len*1.15), lbl, color=col, fontsize=8, alpha=0.7)

    # --- formatting ---
    lim = plot_radius
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_zlim(-lim, lim)

    ax.set_box_aspect([1, 1, 1])

    ax.set_xlabel("x (km)", labelpad=6, fontsize=9)
    ax.set_ylabel("y (km)", labelpad=6, fontsize=9)
    ax.set_zlabel("z (km)", labelpad=6, fontsize=9)

    ax.set_facecolor("#0d0d1a")
    fig.patch.set_facecolor("#0d0d1a")
    ax.tick_params(colors="gray", labelsize=7)
    for pane in [ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane]:
        pane.fill = False
        pane.set_edgecolor("#333355")

    ax.grid(True, color="#222244", linewidth=0.4)

    e_label = ("hyperbolic" if e > 1 else "elliptic" if e < 1 else "parabolic")
    title = (
        f"{body_name} approach  —  {e_label}\n"
        f"e = {e:.4f}   rp = {rp:.0f} km   "
        f"alt = {elements['periapsis_altitude_km']:.0f} km   "
        f"i = {elements['inclination_deg']:.1f}°"
    )
    ax.set_title(title, color="white", fontsize=10, pad=10)
    ax.legend(loc="upper left", fontsize=8,
              facecolor="#1a1a2e", edgecolor="#333355", labelcolor="white")

    plt.tight_layout()
    return fig, ax