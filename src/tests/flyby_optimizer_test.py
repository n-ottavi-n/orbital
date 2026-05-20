from flyby_optimizer import FlybyOptimizer
import planetary_data
import tools as t

t.load_solar_system_kernels()

t0 = 'Jan 03, 2024, 11:20 UTC'
tf = 'Jan 03, 2025, 11:20 UTC'

res_orbit = t.get_resonant_orbit("EUROPA", [1, 6], t0, planetary_data.jupiter)

perturbations = [[
    planetary_data.io,
    planetary_data.europa,
    planetary_data.ganymede,
    planetary_data.callisto
]]

optimizer = FlybyOptimizer(
    target_body     = "EUROPA",
    central_body    = planetary_data.jupiter,
    start_date      = t0,
    end_date        = tf,
    resonant_orbit  = res_orbit,
    perturbations   = perturbations,
    alt_min_km      = 100,
    alt_max_km      = 5000,
    e_delta         = 0.05,
    steps_coarse    = 3000,
    steps_fine      = 10000,
    steps_final     = 100000,
    fine_delta_deg  = 20.0,
    fine_e_delta    = 0.01,
)

optimizer.optimize()
optimizer.plot_results(save_html="europa_flyby_opt.html")