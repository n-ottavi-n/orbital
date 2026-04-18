import numpy as np
import spiceypy as spice
from scipy.optimize import minimize

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
    max_iter=80
):
    """
    Generic periapsis optimizer using RTN departure-frame correction.

    Optimizes:
        dv = a*T + b*R + c*N

    T = along-track
    R = radial
    N = normal to departure plane
    """

    # --------------------------------------------------
    # Nominal Lambert solution
    # --------------------------------------------------
    sim = lambert_interface_class(
        origin, dest,
        start_date,
        arrival_date,
        end_date,
        steps,
        perturbations
    )

    sim.solve()

    v0_nominal = sim.v0.copy()
    r0 = sim.start_r.copy()

    # --------------------------------------------------
    # Build RTN basis
    # --------------------------------------------------
    Rhat = r0 / np.linalg.norm(r0)

    That = v0_nominal / np.linalg.norm(v0_nominal)

    Nhat = np.cross(Rhat, That)
    Nhat /= np.linalg.norm(Nhat)

    # re-orthogonalize tangent
    That = np.cross(Nhat, Rhat)
    That /= np.linalg.norm(That)

    # --------------------------------------------------
    # Body radius
    # --------------------------------------------------
    try:
        radii = spice.bodvrd(dest, "RADII", 3)[1]
        body_radius = np.mean(radii)
    except:
        body_radius = 3396.0   # fallback Mars-like default

    rp_target = body_radius + periapsis_altitude_km

    print("Target periapsis:", rp_target, "km")

    # --------------------------------------------------
    # Helper: smooth periapsis estimate
    # --------------------------------------------------
    def estimate_periapsis(pos_sc, pos_tgt, vel_sc, times, mu_body):
        """
        Compute osculating periapsis radius at closest approach,
        using proper orbital mechanics.
        """
        # get Mars velocity from SPICE at each step to find true relative state
        rel = pos_sc - pos_tgt
        ranges = np.linalg.norm(rel, axis=1)
        i = np.argmin(ranges)
        epoch = times[i]

        # SPICE state for true relative velocity
        state_body, _ = spice.spkezr(dest, epoch, "ECLIPJ2000", "NONE", "SUN")
        r_body_spice = state_body[:3]
        v_body_spice = state_body[3:]

        r = pos_sc[i] - r_body_spice
        v = vel_sc[i] - v_body_spice

        rmag = np.linalg.norm(r)
        vmag = np.linalg.norm(v)
        energy = vmag ** 2 / 2 - mu_body / rmag
        h = np.cross(r, v)
        hmag = np.linalg.norm(h)
        e_vec = np.cross(v, h) / mu_body - r / rmag
        e = np.linalg.norm(e_vec)

        if abs(e - 1.0) < 1e-10:
            return hmag ** 2 / mu_body / 2
        elif abs(e - 1.0) < 1e-6:
            return hmag ** 2 / mu_body / (1 + e)
        else:
            a = -mu_body / (2 * energy)
            return a * (1 - e)

    # --------------------------------------------------
    # Objective
    # --------------------------------------------------
    def objective(u):
        """
        u = [a,b,c] in km/s
        """

        try:
            a, b, c = u

            dv = a * That + b * Rhat + c * Nhat

            trial = lambert_interface_class(
                origin, dest,
                start_date,
                arrival_date,
                end_date,
                steps,
                perturbations
            )

            trial.solve()

            trial.v0 = v0_nominal + dv

            trial.propagate()

            pos_sc = trial.int_props[0]
            pos_tgt = np.array([
                spice.spkpos(dest, t, "ECLIPJ2000", "NONE", "SUN")[0]
                for t in trial.times
            ])

            rel = pos_sc - pos_tgt
            vel_sc = trial.int_vels[0]

            rp = estimate_periapsis(pos_sc, pos_tgt, vel_sc, trial.times, mu_body)

            err = rp - rp_target

            # small penalty on huge maneuvers
            penalty = 10.0 * np.dot(u, u)

            J = err**2 + penalty

            print(
                "u =", np.round(u, 4),
                "rp =", round(rp, 1),
                "err =", round(err, 1)
            )

            return J

        except Exception as e:
            print("FAILED:", e)
            return 1e30

    # --------------------------------------------------
    # Initial simplex (50 m/s basis directions)
    # --------------------------------------------------
    simplex = np.array([
        [0.00, 0.00, 0.00],
        [0.05, 0.00, 0.00],
        [0.00, 0.05, 0.00],
        [0.00, 0.00, 0.05]
    ])

    result = minimize(
        objective,
        x0=np.zeros(3),
        method="Nelder-Mead",
        options={
            "maxiter": max_iter,
            "initial_simplex": simplex,
            "xatol": 1e-4,
            "fatol": 1e-2
        }
    )

    # --------------------------------------------------
    # Final trajectory
    # --------------------------------------------------
    a, b, c = result.x
    dv_best = a * That + b * Rhat + c * Nhat

    final = lambert_interface_class(
        origin, dest,
        start_date,
        arrival_date,
        end_date,
        steps,
        perturbations
    )

    final.solve()

    final.v0 = v0_nominal + dv_best
    final.propagate()

    pos_sc = final.int_props[0]
    pos_tgt = np.array([
        spice.spkpos(dest, t, "ECLIPJ2000", "NONE", "SUN")[0]
        for t in final.times
    ])


    rel = pos_sc - pos_tgt

    ranges = np.linalg.norm(rel, axis=1)

    vel_sc = final.int_vels[0]

    rp_final = estimate_periapsis(pos_sc, pos_tgt, vel_sc, final.times, mu_body)


    return {
        "success": result.success,
        "message": result.message,
        "u_rtn": result.x,
        "dv_xyz": dv_best,
        "r0": pos_sc[0],
        "v0_corrected": final.v0,
        "periapsis_km": rp_final,
        "target_km": rp_target,
        "error_km": rp_final - rp_target,
        "sim_object": final
    }