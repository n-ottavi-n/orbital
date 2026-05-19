from flyby_optimizer import FlybyOptimizer
import planetary_data
import tools as t

t.load_solar_system_kernels()

t0 = 'Jan 03, 2024, 11:20 UTC'
tf = 'Jan 03, 2025, 11:20 UTC'

res_orbit = t.get_resonant_orbit("EUROPA", [1, 5], t0, planetary_data.jupiter)

perturbations = [[
    planetary_data.io,
    planetary_data.europa,
    planetary_data.ganymede,
    planetary_data.callisto
]]

optimizer = FlybyOptimizer(
    target_body    = "EUROPA",
    central_body   = planetary_data.jupiter,
    start_date     = t0,
    end_date       = tf,
    resonant_orbit = res_orbit,
    perturbations  = perturbations,
    alt_min_km     = 100,
    alt_max_km     = 5000,
    a_bounds_pct   = 2.0,
    steps_opt      = 3000,
    steps_final    = 100000,
    popsize        = 15,
    maxiter        = 100,
)

optimizer.optimize()
optimizer.plot_results(save_html="europa_flyby_opt.html")