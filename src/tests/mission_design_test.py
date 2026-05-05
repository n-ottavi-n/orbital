
from mission_design import MissionDesign
import planetary_data as pd

mission = MissionDesign(
    origin='EARTH',
    destination='JUPITER',
    launch_date='JAN 31, 2030, 00:00 UTC',
    arrival_date='AUG  21, 2032, 00:00 UTC',
    periapsis_alt_km=10000,
    inc_target_deg=20,
    parking_alt_km=400,
    arrival_window_days=60
)

mission.solve()
mission.report()
print(mission.arrival)
mission.plot(save_report=True, report_path="earth_jupiter_2030.png")