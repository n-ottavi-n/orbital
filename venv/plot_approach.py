import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def plot_approach(elements: dict, body_name: str = "target", body_radius_km: float = None, mu_body: float = None):

    r_vec = np.array(elements["r_vec"])
    v_vec = np.array(elements["v_vec"])
    e     = elements["eccentricity"]
    a     = elements["semi_major_axis_km"]
    rp    = elements["periapsis_km"]
    nu_ca = np.radians(elements["true_anomaly_deg"])
    vinf  = elements["vinf_km_s"]

    if body_radius_km is None:
        body_radius_km = rp - elements["periapsis_altitude_km"]

    # --------------------------------------------------
    # mu recovery
    # --------------------------------------------------
    h_vec = np.cross(r_vec, v_vec)
    hmag  = np.linalg.norm(h_vec)
    h_hat = h_vec / hmag

    rmag  = np.linalg.norm(r_vec)
    vmag  = np.linalg.norm(v_vec)

    if mu_body is None:
        if abs(e - 1.0) > 1e-6:
            p_conic = abs(a) * abs(1 - e**2)
        else:
            p_conic = hmag**2
        mu_body = hmag**2 / p_conic

    # --------------------------------------------------
    # capture delta-V
    # --------------------------------------------------
    v_periapsis = np.sqrt(vinf**2 + 2 * mu_body / rp)
    v_parabolic  = np.sqrt(2 * mu_body / rp)              # escape velocity at periapsis
    delta_v_capture = v_periapsis - v_parabolic

    # --------------------------------------------------
    # perifocal frame
    # --------------------------------------------------
    e_vec = np.cross(v_vec, h_vec) / mu_body - r_vec / rmag
    e_hat = e_vec / np.linalg.norm(e_vec) if e > 1e-8 else r_vec / rmag
    q_hat = np.cross(h_hat, e_hat)

    # --------------------------------------------------
    # true anomaly sweep
    # --------------------------------------------------
    if e < 1.0:
        nu_arr = np.linspace(0, 2 * np.pi, 500)
    elif abs(e - 1.0) < 1e-6:
        nu_arr = np.linspace(-np.radians(170), np.radians(170), 500)
    else:
        nu_inf = np.arccos(-1.0 / e)
        nu_lim = nu_inf - np.radians(5)
        nu_arr = np.linspace(-nu_lim, nu_lim, 500)

    # --------------------------------------------------
    # orbit curve
    # --------------------------------------------------
    p_orb = abs(a) * abs(1 - e**2) if abs(e - 1.0) > 1e-6 else hmag**2 / mu_body

    def radius_at(nu):
        return p_orb / (1 + e * np.cos(nu))

    pts = np.array([
        radius_at(nu) * (np.cos(nu) * e_hat + np.sin(nu) * q_hat)
        for nu in nu_arr
    ])

    asymptote = None
    if e > 1.0:
        nu_inf = np.arccos(-1.0 / e)
        asymptote = np.cos(-nu_inf) * e_hat + np.sin(-nu_inf) * q_hat

    plot_radius = rp * 1.8
    ca_pt = radius_at(nu_ca) * (np.cos(nu_ca) * e_hat + np.sin(nu_ca) * q_hat)

    # --------------------------------------------------
    # figure
    # --------------------------------------------------
    fig = plt.figure(figsize=(9, 8))
    ax  = fig.add_subplot(111, projection="3d")

    # --- body sphere ---
    u, v = np.mgrid[0:2*np.pi:60j, 0:np.pi:30j]
    R = body_radius_km
    bx = R * np.cos(u) * np.sin(v)
    by = R * np.sin(u) * np.sin(v)
    bz = R * np.cos(v)
    ax.plot_surface(bx, by, bz, color="steelblue", alpha=0.5,      # CHANGED: alpha up slightly for white bg
                    linewidth=0, zorder=1)
    th = np.linspace(0, 2*np.pi, 200)
    ax.plot(R*np.cos(th), R*np.sin(th), np.zeros_like(th),
            color="steelblue", lw=0.7, alpha=0.8)

    # --- orbit curve ---
    ax.plot(pts[:, 0], pts[:, 1], pts[:, 2],
            color="royalblue", lw=1.8, label="approach trajectory", zorder=3)  # CHANGED: gold→royalblue for white bg

    # --- periapsis marker ---
    peri_pt = rp * e_hat
    ax.scatter(*peri_pt, color="tomato", s=60, zorder=5,
               label=f"periapsis  {rp:.0f} km  (alt {elements['periapsis_altitude_km']:.0f} km)")

    # --- velocity vector at periapsis ---
    # at periapsis, velocity is purely tangential (perpendicular to e_hat)
    v_peri_dir = q_hat * np.sign(np.dot(v_vec, q_hat))
    vscale = plot_radius * 0.18 / np.linalg.norm(v_vec)
    ax.quiver(*peri_pt, *(v_peri_dir * np.linalg.norm(v_vec) * vscale),
              color="orange", linewidth=2, arrow_length_ratio=0.18,
              label=f"v periapsis  {np.sqrt(vinf ** 2 + 2 * mu_body / rp):.2f} km/s")

    # --- v∞ asymptote ---
    if asymptote is not None:
        ray_start = -asymptote * plot_radius * 0.9
        ray_end   = -asymptote * body_radius_km * 1.5
        ax.plot([ray_start[0], ray_end[0]],
                [ray_start[1], ray_end[1]],
                [ray_start[2], ray_end[2]],
                "--", color="darkorange", lw=1.2, alpha=0.8,
                label=f"v∞ = {vinf:.2f} km/s")

    # --- body-frame axes ---
    ax_len = body_radius_km * 1.6
    for vec, col, lbl in zip(
        [np.array([1,0,0]), np.array([0,1,0]), np.array([0,0,1])],
        ["red", "green", "royalblue"],
        ["x", "y", "z"]
    ):
        ax.quiver(0, 0, 0, *(vec*ax_len), color=col, linewidth=0.8,
                  arrow_length_ratio=0.25, alpha=0.6)
        ax.text(*(vec*ax_len*1.15), lbl, color=col, fontsize=8, alpha=0.8)

    # --- formatting ---
    lim = plot_radius
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_zlim(-lim, lim)
    ax.set_box_aspect([1, 1, 1])

    ax.set_xlabel("x (km)", labelpad=6, fontsize=9)
    ax.set_ylabel("y (km)", labelpad=6, fontsize=9)
    ax.set_zlabel("z (km)", labelpad=6, fontsize=9)

    ax.set_facecolor("white")                                          # CHANGED
    fig.patch.set_facecolor("white")                                   # CHANGED
    ax.tick_params(colors="black", labelsize=7)                        # CHANGED
    for pane in [ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane]:
        pane.fill = False
        pane.set_edgecolor("#cccccc")                                  # CHANGED

    ax.grid(True, color="#dddddd", linewidth=0.4)                      # CHANGED

    e_label = ("hyperbolic" if e > 1 else "elliptic" if e < 1 else "parabolic")
    title = (
        f"{body_name} approach  —  {e_label}\n"
        f"e = {e:.4f}   rp = {rp:.0f} km   "
        f"alt = {elements['periapsis_altitude_km']:.0f} km   "
        f"i = {elements['inclination_deg']:.1f}°   "
        f"ΔV capture (parabolic) = {delta_v_capture:.3f} km/s"                    # NEW
    )
    ax.set_title(title, color="black", fontsize=10, pad=10)            # CHANGED
    ax.legend(loc="upper left", fontsize=8,
              facecolor="white", edgecolor="#cccccc",
              labelcolor="black")                                       # CHANGED

    plt.tight_layout()
    return fig, ax