import numpy as np
import spiceypy as spice
from datetime import timedelta, datetime
from trajectory_utils import lambert_branches
from trajectory_utils import compute_c3
from trajectory_utils import hohmann_phase_angle


def optimize_c3(lambert_interface_class, origin, dest, base_start, base_arrival,
                days_range=3, step_hours=6):

    best = {
        "c3": float("inf"),
        "start": None,
        "arrival": None
    }

    base_start_et = spice.str2et(base_start)
    base_arrival_et = spice.str2et(base_arrival)

    step = step_hours * 3600

    mu = 1.32712440018e11  # Sun gravitational parameter (km^3/s^2)

    for dt0 in np.arange(-days_range*86400, days_range*86400 + step, step):
        for dt1 in np.arange(-days_range*86400, days_range*86400 + step, step):

            start_et = base_start_et + dt0
            arrival_et = base_arrival_et + dt1

            if arrival_et <= start_et:
                continue

            ##### PHASE ANGLE FILTERING speeds up everything #########
            # get positions first
            r1, _ = spice.spkpos(origin, start_et, 'ECLIPJ2000', 'NONE', 'SUN')
            r2, _ = spice.spkpos(dest, start_et, 'ECLIPJ2000', 'NONE', 'SUN')

            a1 = np.linalg.norm(r1)
            a2 = np.linalg.norm(r2)

            # compute phase angle
            dot = np.dot(r1, r2) / (np.linalg.norm(r1) * np.linalg.norm(r2))
            dot = np.clip(dot, -1, 1)
            phase_angle = np.degrees(np.arccos(dot))

            cross = np.cross(r1, r2)
            if cross[2] < 0:
                phase_angle = 360 - phase_angle

            # compute expected angle
            phi_nom = hohmann_phase_angle(a1, a2, mu)

            # filter
            diff = abs((phase_angle - phi_nom + 180) % 360 - 180)

            if diff > 40:
                print("phase angle outside range")
                continue

            ##### END PHASE ANGLE FILTERING ########

            # Convert back to UTC string for your interface
            start_str = spice.et2utc(start_et, "C", 0)
            arrival_str = spice.et2utc(arrival_et, "C", 0)

            try:
                sim = lambert_interface_class(
                    origin, dest,
                    start_str,
                    arrival_str,
                    arrival_str,   # end_date not important
                    steps=1000,
                    perturbations=[]
                )

                sim.solve()

                state, _ = spice.spkezr(origin, start_et, 'ECLIPJ2000', 'NONE', 'SUN')
                v_planet = state[3:]

                solutions = lambert_branches(
                    sim.start_r,
                    sim.end_r,
                    sim.tof,
                    mu,
                    prograde=True
                )

                if not solutions:
                    continue

                # evaluate all branches
                c3_candidates = []

                for sol in solutions:
                    v0 = sol[:3]
                    c3_candidates.append(compute_c3(v0, v_planet))

                c3 = min(c3_candidates)


                if c3 < best["c3"]:
                    best.update({
                        "c3": c3,
                        "start": start_str,
                        "arrival": arrival_str
                    })
                    print("New best:", best)


            except Exception as e:

                print("FAIL:", start_str, arrival_str, e)

                continue

    return best