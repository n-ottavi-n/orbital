import numpy as np
import spiceypy as spice
from scipy.optimize import minimize
from arrival_orbit_elements import arrival_orbit_elements


def optimize_periapsis(
    lambert_interface_class,
    origin,
    dest,
    start_date,
    arrival_date,
    end_date,
    steps,
    perturbations,
    periapsis_altitude_km,
    mu_body,
    inc_target_deg=None,        # NEW: optional inclination target
    arrival_window_days=30,     # NEW: search window around nominal arrival
    max_iter=80,
    x0=None,       # NEW: warm start RTN vector
    dt0=0.0        # NEW: warm start arrival offset in days
):
    # --------------------------------------------------
    # Nominal Lambert solution
    # --------------------------------------------------
    sim = lambert_interface_class(
        origin, dest,
        start_date, arrival_date, end_date,
        steps, perturbations
    )
    sim.solve()

    v0_nominal = sim.v0.copy()
    r0         = sim.start_r.copy()

    # --------------------------------------------------
    # RTN basis at departure
    # --------------------------------------------------
    Rhat = r0 / np.linalg.norm(r0)
    That = v0_nominal / np.linalg.norm(v0_nominal)
    Nhat = np.cross(Rhat, That)
    Nhat /= np.linalg.norm(Nhat)
    That = np.cross(Nhat, Rhat)
    That /= np.linalg.norm(That)

    # --------------------------------------------------
    # Body radius and target periapsis
    # --------------------------------------------------
    try:
        radii = spice.bodvrd(dest, "RADII", 3)[1]
        body_radius = np.mean(radii)
    except:
        body_radius = 3396.0
    rp_target = body_radius + periapsis_altitude_km

    print(f"Target periapsis:   {rp_target:.1f} km")
    if inc_target_deg is not None:
        print(f"Target inclination: {inc_target_deg:.1f} deg")

    # --------------------------------------------------
    # Helper: run a trial simulation and return elements
    # --------------------------------------------------
    def run_trial(u):
        a, b, c, dt_arr = u

        dv = a * That + b * Rhat + c * Nhat
        arrival_et_offset = dt_arr * 86400.0    # days → seconds

        trial = lambert_interface_class(
            origin, dest,
            start_date, arrival_date, end_date,
            steps, perturbations,
            arrival_et_offset=arrival_et_offset
        )
        trial.solve()
        trial.v0 = v0_nominal + dv
        trial.propagate()

        elements = arrival_orbit_elements(trial, dest, mu_body)
        return elements

    # --------------------------------------------------
    # Objective: minimise |dv|
    # --------------------------------------------------
    def objective(u):
        a, b, c, dt_arr = u
        dv_mag = np.sqrt(a**2 + b**2 + c**2)
        return dv_mag

    # --------------------------------------------------
    # Constraints
    # --------------------------------------------------
    def constraint_periapsis(u):
        try:
            elements = run_trial(u)
            err = elements["periapsis_km"] - rp_target
            '''print(f"  rp err: {err:+.1f} km   "
                  f"rp: {elements['periapsis_km']:.1f} km   "
                  f"dt: {u[3]:+.2f} days   "
                  f"dv: {np.linalg.norm(u[:3])*1000:.1f} m/s")'''
            return err
        except Exception as e:
            print("constraint_periapsis failed:", e)
            return -1e10

    constraints = [
        {'type': 'eq', 'fun': constraint_periapsis},
    ]

    if inc_target_deg is not None:
        inc_tol_deg = 2.0  # accept anything within ±2°, adjust to taste

        def constraint_inclination(u):
            try:
                elements = run_trial(u)
                err = elements["inclination_deg"] - inc_target_deg
                '''print(f"  inc err: {err:+.2f} deg   "
                      f"inc: {elements['inclination_deg']:.2f} deg")'''
                return inc_tol_deg - abs(err)  # ineq: must be >= 0
            except Exception as e:
                print("constraint_inclination failed:", e)
                return -1e10

        constraints.append({'type': 'ineq', 'fun': constraint_inclination})



    # --------------------------------------------------
    # Bounds
    # --------------------------------------------------
    bounds = [
        (-0.5,  0.5),                             # a  km/s
        (-0.5,  0.5),                             # b  km/s
        (-0.5,  0.5),                             # c  km/s
        (-arrival_window_days, arrival_window_days)  # dt days
    ]

    # --------------------------------------------------
    # Optimize
    # --------------------------------------------------

    initial_guess = np.array([*(x0 if x0 is not None else [0., 0., 0.]), dt0])


    result = minimize(
        objective,
        x0=initial_guess,
        method='SLSQP',
        bounds=bounds,
        constraints=constraints,
        options={
            'maxiter': max_iter,
            'ftol': 1e-3,
            'eps': 1e-5,
            'disp': True
        }
    )

    # --------------------------------------------------
    # Final trajectory at converged solution
    # --------------------------------------------------
    a, b, c, dt_arr = result.x
    dv_best = a * That + b * Rhat + c * Nhat
    arrival_et_offset = dt_arr * 86400.0

    final = lambert_interface_class(
        origin, dest,
        start_date, arrival_date, end_date,
        steps, perturbations,
        arrival_et_offset=arrival_et_offset
    )
    final.solve()
    final.v0 = v0_nominal + dv_best
    final.propagate()

    elements_final = arrival_orbit_elements(final, dest, mu_body)

    print(f"\nResult:")
    print(f"  periapsis alt: {elements_final['periapsis_altitude_km']:.2f} km  "
          f"(target {periapsis_altitude_km} km)")
    print(f"  inclination:   {elements_final['inclination_deg']:.2f} deg"
          + (f"  (target {inc_target_deg} deg)" if inc_target_deg else ""))
    print(f"  dv:            {np.linalg.norm(dv_best)*1000:.2f} m/s")
    print(f"  arrival shift: {dt_arr:+.2f} days")

    return {
        "success":               result.success,
        "message":               result.message,
        "u_rtn":                 result.x[:3],
        "dt_arrival_days":       dt_arr,           # NEW
        "dv_xyz":                dv_best,
        "dv_magnitude_ms":       np.linalg.norm(dv_best) * 1000,
        "r0":                    final.int_props[0][0],
        "v0_corrected":          final.v0,
        "periapsis_km":          elements_final["periapsis_km"],
        "periapsis_altitude_km": elements_final["periapsis_altitude_km"],
        "inclination_deg":       elements_final["inclination_deg"],  # NEW
        "target_periapsis_km":   rp_target,
        "target_inc_deg":        inc_target_deg,                     # NEW
        "sim_object":            final,
    }