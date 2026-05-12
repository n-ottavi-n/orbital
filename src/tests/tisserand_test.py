from tisserand import tisserand_graph

bodies = ["EARTH", "JUPITER", "SATURN"]

tisserand_graph(
    bodies      = bodies,
    vinf_values = [4, 6, 8, 10, 12],
    rp_min_au=0.8,
    rp_max_au=15.0,
    save_html   = "reports/tisserand.html"
)