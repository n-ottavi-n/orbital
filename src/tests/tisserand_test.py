from tisserand import tisserand_graph

bodies = ["VENUS", "EARTH", "JUPITER", "SATURN"]

tisserand_graph(
    bodies      = bodies,
    vinf_values = [4, 6, 8, 10, 12],
    rp_min_au=0.3,
    rp_max_au=20.0,
    save_html   = "reports/tisserand.html"
)