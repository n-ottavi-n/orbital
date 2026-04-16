import numpy as np
import spiceypy as spice
from datetime import timedelta, datetime
from trajectory_utils import lambert_branches
from trajectory_utils import compute_c3


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