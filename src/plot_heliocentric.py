import numpy as np
import spiceypy as spice
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

AU = 1.495978707e8   # km per AU


def plot_heliocentric(sim, departure: dict, arrival: dict,
                      origin: str = "EARTH", destination: str = "MARS"):

    # --------------------------------------------------
    # spacecraft trajectory (heliocentric, in AU)
    # --------------------------------------------------
    traj = sim.int_props[0] / AU

    # --------------------------------------------------
    # origin and destination positions from SPICE
    # over the full transfer window
    # --------------------------------------------------
    times = np.array(sim.times)

    r_origin = np.array([
        spice.spkpos(origin, float(t), "ECLIPJ2000", "NONE", "SUN")[0]
        for t in times
    ]) / AU

    r_dest = np.array([
        spice.spkpos(destination, float(t), "ECLIPJ2000", "NONE", "SUN")[0]
        for t in times
    ]) / AU

    # --------------------------------------------------
    # departure and arrival points
    # --------------------------------------------------
    dep_et  = departure["departure_et"]
    arr_et  = arrival["epoch_et"]

    r_dep = np.array(
        spice.spkpos(origin, dep_et, "ECLIPJ2000", "NONE", "SUN")[0]
    ) / AU

    r_arr = np.array(
        spice.spkpos(destination, arr_et, "ECLIPJ2000", "NONE", "SUN")[0]
    ) / AU

    # spacecraft position at departure and arrival
    dep_idx = np.argmin(np.abs(times - dep_et))
    arr_idx = np.argmin(np.abs(times - arr_et))

    sc_dep = traj[dep_idx]
    sc_arr = traj[arr_idx]

    # --------------------------------------------------
    # Sun
    # --------------------------------------------------
    fig = plt.figure(figsize=(9, 8))
    ax  = fig.add_subplot(111, projection="3d")

    ax.scatter(0, 0, 0, color="gold", s=100, zorder=5, label="Sun")

    # --------------------------------------------------
    # origin orbit (full period approximation)
    # --------------------------------------------------
    ax.plot(r_origin[:, 0], r_origin[:, 1], r_origin[:, 2],
            color="steelblue", lw=0.8, alpha=0.5,
            label=origin.capitalize())

    # --------------------------------------------------
    # destination orbit
    # --------------------------------------------------
    ax.plot(r_dest[:, 0], r_dest[:, 1], r_dest[:, 2],
            color="tomato", lw=0.8, alpha=0.5,
            label=destination.capitalize())

    # --------------------------------------------------
    # transfer trajectory
    # --------------------------------------------------
    ax.plot(traj[:, 0], traj[:, 1], traj[:, 2],
            color="royalblue", lw=1.8,
            label="transfer trajectory", zorder=3)

    # --------------------------------------------------
    # departure point
    # --------------------------------------------------
    ax.scatter(*r_dep, color="steelblue", s=60, zorder=6)
    ax.scatter(*sc_dep, color="royalblue", s=60, zorder=6,
               label=f"departure  {departure['departure_utc']}")

    # --------------------------------------------------
    # arrival point
    # --------------------------------------------------
    ax.scatter(*r_arr, color="tomato", s=60, zorder=6)
    ax.scatter(*sc_arr, color="lime", s=80, zorder=6,
               label=f"arrival  {arrival['arrival_utc']}")

    # --------------------------------------------------
    # formatting
    # --------------------------------------------------
    rp = np.linalg.norm(r_dest[0])
    plot_radius = rp * 1.2
    lim = plot_radius
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_zlim(-lim, lim)

    ax.set_xlabel("x (AU)", labelpad=6, fontsize=9)
    ax.set_ylabel("y (AU)", labelpad=6, fontsize=9)
    ax.set_zlabel("z (AU)", labelpad=6, fontsize=9)

    ax.set_facecolor("white")
    fig.patch.set_facecolor("white")
    ax.tick_params(colors="black", labelsize=7)
    for pane in [ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane]:
        pane.fill = False
        pane.set_edgecolor("#cccccc")
    ax.grid(True, color="#dddddd", linewidth=0.4)
    ax.set_box_aspect([1, 1, 1])

    tof = (arr_et - dep_et) / 86400
    title = (
        f"{origin.capitalize()} → {destination.capitalize()} "
        f"heliocentric transfer\n"
        f"TOF = {tof:.1f} days   "
        f"C3 = {departure['c3_km2_s2']:.2f} km²/s²   "
        f"v∞ arr = {arrival['vinf_km_s']:.2f} km/s"
    )
    ax.set_title(title, color="black", fontsize=10, pad=10)
    ax.legend(loc="upper left", fontsize=8,
              facecolor="white", edgecolor="#cccccc", labelcolor="black")

    plt.tight_layout()
    return fig, ax