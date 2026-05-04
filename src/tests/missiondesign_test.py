
from mission_design import MissionDesign

mission = MissionDesign(
    origin='EARTH',
    destination='MARS',
    launch_date='JAN 30, 2031, 23:00 UTC',
    arrival_date='JUL 12, 2031, 00:00 UTC',
    periapsis_alt_km=300,
    inc_target_deg=25,
    parking_alt_km=400,
)

mission.solve()
mission.report()
mission.plot(save_report=True, report_path="earth_mars_2031.png")