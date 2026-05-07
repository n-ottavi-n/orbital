import numpy as np
import plotly.graph_objects as go
import planetary_data as pd
import spiceypy as spice
from tools import load_solar_system_kernels

def get_semi_major_axis_km(body_name, mu_sun, epoch=0.0):
    """
    Compute mean semi-major axis of a body's heliocentric orbit from SPICE.
    Uses vis-viva at a single epoch — accurate enough for Tisserand analysis.
    """
    load_solar_system_kernels()
    state, _ = spice.spkezr(body_name, epoch, "ECLIPJ2000", "NONE", "SUN")
    r = np.linalg.norm(state[:3])
    v = np.linalg.norm(state[3:])
    a = 1.0 / (2.0/r - v**2/mu_sun)
    return a

AU = 1.495978707e8  # km per AU


def tisserand_graph(
    bodies,
    vinf_values,
    mu_sun=1.32712440018e11,
    rp_min_au=0.3,
    rp_max_au=10.0,
    n_points=500,
    save_html=None
):
    """
    Plot a Tisserand graph for a list of bodies and v∞ values.

    Parameters
    ----------
    bodies      : list of body name strings e.g. ["VENUS", "EARTH", "MARS", "JUPITER"]
    vinf_values : list of v∞ values in km/s e.g. [2, 4, 6, 8, 10]
    mu_sun      : gravitational parameter of the Sun (km³/s²)
    rp_min_au   : minimum periapsis to plot (AU)
    rp_max_au   : maximum periapsis to plot (AU)
    n_points    : number of theta points per curve
    save_html   : path to save interactive HTML (optional)
    """

    # --------------------------------------------------
    # color palette — one color per body
    # --------------------------------------------------
    colors = [
        '#e6194b', '#3cb44b', '#4363d8', '#f58231',
        '#911eb4', '#42d4f4', '#f032e6', '#bfef45',
    ]

    # --------------------------------------------------
    # line styles — one per v∞ value
    # --------------------------------------------------
    dash_styles = ['solid', 'dash', 'dot', 'dashdot', 'longdash',
                   'longdashdot', 'solid', 'dash']

    fig = go.Figure()

    # --------------------------------------------------
    # reference lines
    # --------------------------------------------------

    # diagonal: rp = ra (circular orbits)
    diag = np.linspace(rp_min_au, rp_max_au, 100)
    fig.add_trace(go.Scatter(
        x=diag, y=diag,
        mode='lines',
        line=dict(color='lightgrey', width=1, dash='dot'),
        name='circular orbit',
        showlegend=True,
        hoverinfo='skip',
    ))

    # --------------------------------------------------
    # body vertical reference lines + annotations
    # --------------------------------------------------
    for idx, body_name in enumerate(bodies):
        body     = pd.get_body(body_name)
        a_p_km = get_semi_major_axis_km(body_name, mu_sun)

        # if semi-major axis not in planetary_data, skip reference line
        # (we'll get it from the curve computation below)
        color = colors[idx % len(colors)]

        if a_p_km is not None:
            a_p_au = a_p_km / AU
            fig.add_vline(
                x=a_p_au,
                line=dict(color=color, width=1, dash='dot'),
            )
            fig.add_annotation(
                x=a_p_au, y=rp_max_au * 0.97,
                text=body_name.capitalize(),
                showarrow=False,
                font=dict(color=color, size=11, family='Arial Black'),
                xanchor='center',
            )

    # --------------------------------------------------
    # Tisserand curves
    # --------------------------------------------------
    # parametric sweep over theta (angle between v∞ and body velocity)
    theta_arr = np.linspace(0, 2 * np.pi, n_points)

    for idx, body_name in enumerate(bodies):
        body   = pd.get_body(body_name)
        color  = colors[idx % len(colors)]

        # get semi-major axis — from planetary_data or compute from mu
        a_p_km = get_semi_major_axis_km(body_name, mu_sun)
        if a_p_km is None:
            print(f"Warning: no semi-major axis for {body_name}, skipping")
            continue

        a_p    = a_p_km / AU          # in AU for plotting
        v_p    = np.sqrt(mu_sun / a_p_km)   # circular velocity km/s

        for v_idx, vinf in enumerate(vinf_values):

            rp_list = []
            ra_list = []

            for theta in theta_arr:

                # spacecraft heliocentric velocity components at flyby
                v_t = v_p + vinf * np.cos(theta)   # tangential
                v_r = vinf * np.sin(theta)          # radial

                # specific orbital energy
                v2  = v_t**2 + v_r**2
                E   = v2 / 2 - mu_sun / a_p_km

                # only elliptic orbits (E < 0)
                if E >= 0:
                    continue

                # semi-major axis
                a_km = -mu_sun / (2 * E)

                # specific angular momentum
                h = a_p_km * v_t

                # eccentricity
                e2 = 1 - h**2 / (mu_sun * a_km)
                if e2 < 0:
                    e2 = 0
                e = np.sqrt(e2)

                # periapsis and apoapsis in AU
                rp_au = a_km * (1 - e) / AU
                ra_au = a_km * (1 + e) / AU

                # filter to plot range
                if rp_au < rp_min_au or ra_au > rp_max_au * 3:
                    continue

                rp_list.append(rp_au)
                ra_list.append(ra_au)

            if len(rp_list) < 2:
                continue

            # sort by rp for clean curve
            sorted_pairs = sorted(zip(rp_list, ra_list))
            rp_arr = np.array([p[0] for p in sorted_pairs])
            ra_arr = np.array([p[1] for p in sorted_pairs])

            # clip ra to plot range
            ra_arr = np.clip(ra_arr, rp_min_au, rp_max_au)

            # only show legend entry for first v∞ of each body
            show_legend = (v_idx == 0)

            fig.add_trace(go.Scatter(
                x=rp_arr,
                y=ra_arr,
                mode='lines',
                line=dict(
                    color=color,
                    width=1.5,
                    dash=dash_styles[v_idx % len(dash_styles)]
                ),
                name=body_name.capitalize(),
                legendgroup=body_name,
                showlegend=show_legend,
                hovertemplate=(
                    f"<b>{body_name.capitalize()}</b><br>"
                    f"v∞ = {vinf:.1f} km/s<br>"
                    "rp = %{x:.2f} AU<br>"
                    "ra = %{y:.2f} AU<br>"
                    "<extra></extra>"
                ),
            ))

            # v∞ label at the right end of the curve
            if len(rp_arr) > 0:
                label_x = rp_arr[-1]
                label_y = ra_arr[-1]
                fig.add_annotation(
                    x=label_x, y=label_y,
                    text=f"{vinf:.0f}",
                    showarrow=False,
                    font=dict(color=color, size=8),
                    xanchor='left',
                )

    for v_idx, vinf in enumerate(vinf_values):
        fig.add_trace(go.Scatter(
            x=[None], y=[None],
            mode='lines',
            line=dict(
                color='grey',
                width=1.5,
                dash=dash_styles[v_idx % len(dash_styles)]
            ),
            name=f"v∞ = {vinf:.0f} km/s",
            legendgroup=f"vinf_{v_idx}",
            showlegend=True,
        ))

    # --------------------------------------------------
    # layout
    # --------------------------------------------------
    vinf_str = ", ".join([f"{v:.0f}" for v in vinf_values])

    fig.update_layout(
        title=dict(
            text=(f"Tisserand Graph — "
                  f"{', '.join([b.capitalize() for b in bodies])}<br>"
                  f"<sup>v∞ = {vinf_str} km/s   "
                  f"planar, circular orbit assumption</sup>"),
            font=dict(size=14, color='black'),
            x=0.5,
            xanchor='center',
        ),
        xaxis=dict(
            title=dict(text="Periapsis rp (AU)", font=dict(color='black')),
            range=[rp_min_au, rp_max_au],
            tickfont=dict(color='black'),
            gridcolor='lightgrey',
            showgrid=True,
            zeroline=False,
        ),
        yaxis=dict(
            title=dict(text="Apoapsis ra (AU)", font=dict(color='black')),
            range=[rp_min_au, rp_max_au],
            tickfont=dict(color='black'),
            gridcolor='lightgrey',
            showgrid=True,
            zeroline=False,
        ),
        paper_bgcolor='white',
        plot_bgcolor='white',
        font=dict(color='black'),
        width=900,
        height=800,
        legend=dict(
            x=1.02, y=1,
            bordercolor='lightgrey',
            borderwidth=1,
            font=dict(color='black'),
        ),
    )

    if save_html:
        fig.write_html(save_html)
        print(f"Saved to {save_html}")

    fig.show()
    return fig
