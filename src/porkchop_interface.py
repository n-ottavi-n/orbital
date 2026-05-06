import numpy as np
import spiceypy as spice
import matplotlib.pyplot as plt
from trajectory_utils import lambert_branches
from lambert_interface import lambert_interface
from datetime import datetime
from trajectory_utils import hohmann_phase_angle
import time
from trajectory_utils import compute_arrival_inclination


def et_to_datetime(et):
    utc_str = spice.et2utc(et, "ISOC", 0)  # e.g. 2028-12-12T00:00:00
    return datetime.strptime(utc_str, "%Y-%m-%dT%H:%M:%S")



class porkchop_interface:
    def __init__(self,origin_body, target_body, launch_date, arrival_date, max_c3, span=90, grid_res=5, mu=1.32712440018e11):

        base_launch = spice.str2et(launch_date)
        base_arrival = spice.str2et(arrival_date)
        '''
        max_tof # in days solutions longer than this will not be tested
        max_c3 # c3 values above this will be masked
        span # will look at launch and arrival dates ±span days around the base_launch and base_arrival
        grid_res #grid resolution in days
        '''
        launch_span = span * 86400
        arrival_span = span * 86400

        #grid_res = 5 #grid resolution in days
        step = grid_res * 3600 * 24

        self.mu = mu
        self.origin_body = origin_body
        self.target_body = target_body

        self.max_c3 = max_c3


        self.launches = np.arange(base_launch - launch_span,
                              base_launch + launch_span,
                              step)

        self.arrivals = np.arange(base_arrival - arrival_span,
                              base_arrival + arrival_span,
                              step)
        
        self.max_tof = (self.arrivals[-1] - self.launches[0]) / 86400



    def plot(self):

        launches = self.launches
        arrivals = self.arrivals

        mu = self.mu
        max_c3 = self.max_c3
        max_tof = self.max_tof

        C3 = np.zeros((len(launches), len(arrivals)))
        C3[:] = np.nan

        VINF_ARR = np.zeros_like(C3)
        VINF_ARR[:] = np.nan

        INC_MIN = np.zeros_like(C3)
        INC_MIN[:] = np.nan

        #for testing filter
        tested=0
        rejected=0
        passed=0
        lambert_time = 0
        spice_time = 0
        plot_time = 0

        for i, t0 in enumerate(launches):
            for j, t1 in enumerate(arrivals):

                tof_days = (t1-t0)/86400

                if tof_days > self.max_tof:
                    continue

                ##### PHASE ANGLE FILTERING speeds up everything #########
                # get positions first
                t=time.time()
                r1, _ = spice.spkpos(self.origin_body, t0, 'ECLIPJ2000', 'NONE', 'SUN')
                r2, _ = spice.spkpos(self.target_body, t0, 'ECLIPJ2000', 'NONE', 'SUN')
                spice_time += time.time() - t

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

                tested+=1
                if diff > 60:
                    rejected+=1
                    continue
                passed+=1

                ##### END PHASE ANGLE FILTERING ########

                try:
                    start_str = spice.et2utc(t0, "C", 0)
                    arr_str   = spice.et2utc(t1, "C", 0)

                    sim = lambert_interface(
                        self.origin_body, self.target_body,
                        start_str,
                        arr_str,
                        arr_str,
                        steps=500,
                        perturbations=[]
                    )

                    t = time.time()
                    sim.solve()
                    lambert_time += time.time() - t

                    state, _ = spice.spkezr(self.origin_body, t0,
                                            'ECLIPJ2000', 'NONE', 'SUN')
                    v_planet = state[3:]

                    state_arr, _ = spice.spkezr(self.target_body, t1,
                                                'ECLIPJ2000', 'NONE', 'SUN')
                    v_planet_arr = state_arr[3:]



                    sols = lambert_branches(sim.start_r,
                                            sim.end_r,
                                            sim.tof,
                                            mu=self.mu)

                    best_c3 = np.inf
                    best_vinf_arr = np.inf
                    best_vinf_arr_vec = None

                    for sol in sols:
                        v0 = sol[:3]
                        v1 = sol[3:]
                        v_inf_dep = v0 - v_planet
                        v_inf_arr_vec = v1 - v_planet_arr
                        v_inf_arr = np.linalg.norm(v_inf_arr_vec)

                        c3 = np.dot(v_inf_dep, v_inf_dep)
                        best_c3 = min(best_c3, c3)

                        if v_inf_arr < best_vinf_arr:
                            best_vinf_arr = v_inf_arr
                            best_vinf_arr_vec = v_inf_arr_vec  # track the vector, not just the magnitude

                    # after the loop:
                    C3[i, j] = best_c3
                    VINF_ARR[i, j] = best_vinf_arr

                    if best_vinf_arr_vec is not None:
                        INC_MIN[i, j] = compute_arrival_inclination(
                            best_vinf_arr_vec, self.target_body, t1
                        )


                except Exception as e:

                    print("FAIL:", i, j, e)

                    continue

        print(f"tested={tested}, passed={passed}, rejected={rejected}")

        # mask high C3
        C3_plot = np.copy(C3)
        C3_plot[C3_plot > max_c3] = np.nan
        C3_masked = np.ma.masked_invalid(C3_plot) # mask very large values for visibility

        t = time.time()
        plt.figure(figsize=(10, 6))

        # plot heatmap
        cf = plt.contourf(launches/86400,
                          arrivals/86400,
                          C3_masked.T,
                     levels=50)

        plt.colorbar(cf, label="C3 (km^2/s^2)")

        # plot C3 contours
        c3_min = np.nanmin(C3_masked)

        # aim for ~10-15 contour lines, rounded to clean intervals
        range_c3 = max_c3 - c3_min
        step = round(range_c3 / 12)   # ~12 levels
        step = max(step, 1)            # minimum step of 1

        c3_levels = np.arange(
            round(c3_min / step) * step,
            round(max_c3 / step) * step + step,
            step
        )

        # C3 contours — white labels with dark outline for visibility
        cs = plt.contour(launches/86400, arrivals/86400, C3_masked.T,
                        levels=c3_levels, linewidths=0.8, colors='yellow')
        plt.clabel(cs, inline=True, fontsize=8, fmt="%.0f",
                inline_spacing=5)
        for text in cs.labelTexts:
            text.set_color('white')
            text.set_backgroundcolor('none')
            text.set_fontweight('bold')

        # plot  arrival v_inf contours
        vinf_min = np.nanmin(VINF_ARR)
        vinf_max = np.nanmax(VINF_ARR)

        range_vinf = vinf_max - vinf_min
        step_vinf  = max(round(range_vinf / 10, 1), 0.5)  # ~10 levels, min step 0.5 km/s

        vinf_levels = np.arange(
            round(vinf_min / step_vinf) * step_vinf,
            round(vinf_max / step_vinf) * step_vinf + step_vinf,
            step_vinf
        )

        cs_vinf = plt.contour(launches/86400,
                              arrivals/86400,
                              VINF_ARR.T,
                              levels=vinf_levels,
                              linestyles='dotted',
                              linewidths=1.0,
                              colors='grey')

        plt.clabel(cs_vinf, inline=True, fontsize=8, fmt="%.1f km/s")

        # plot TOF lines
        TOF = np.zeros_like(C3)

        for i, t0 in enumerate(launches):
            for j, t1 in enumerate(arrivals):
                TOF[i, j] = (t1 - t0) / 86400  # days

        tof_min = np.nanmin(TOF[~np.isnan(C3_plot)])   # only where C3 is valid
        tof_max = self.max_tof

        step_tof = max(round((tof_max - tof_min) / 10 / 10) * 10, 10)  # rounded to nearest 10 days

        tof_levels = np.arange(
            round(tof_min / step_tof) * step_tof,
            tof_max,
            step_tof
        )

        cs_tof = plt.contour(launches/86400,
                             arrivals/86400,
                             TOF.T,
                             levels=tof_levels,
                             linestyles='dashed',
                             colors='red')

        plt.clabel(cs_tof, inline=True, fontsize=8, fmt="%d d", colors='red')

        # plot minimum inclination contours
        inc_levels = np.arange(0, 90, 5)  # every 5°

        cs_inc = plt.contour(
            launches / 86400,
            arrivals / 86400,
            INC_MIN.T,
            levels=inc_levels,
            linestyles='solid',
            linewidths=0.8,
            colors='orange'
        )
        plt.clabel(cs_inc, inline=True, fontsize=7, fmt="%d°")

        # define clean ticks for pretty graph
        tick_idx = np.linspace(0, len(launches)-1, 10, dtype=int)

        launch_ticks = launches[tick_idx] / 86400
        launch_labels = [
            et_to_datetime(launches[i]).strftime("%d %b")
            for i in tick_idx
        ]

        tick_idx_y = np.linspace(0, len(arrivals)-1, 10, dtype=int)

        arrival_ticks = arrivals[tick_idx_y] / 86400
        arrival_labels = [
            et_to_datetime(arrivals[i]).strftime("%d %b")
            for i in tick_idx_y
        ]

        plt.xticks(launch_ticks, launch_labels, rotation=30)
        plt.yticks(arrival_ticks, arrival_labels)

        #legend and labels
        start_year = et_to_datetime(launches[0]).year
        end_year   = et_to_datetime(arrivals[-1]).year

        plt.title(f"{self.origin_body} → {self.target_body} transfer\nLaunch: {start_year} | Arrival: {end_year}")

        plt.xlabel("Launch Date")
        plt.ylabel("Arrival Date")

        plot_time += time.time() - t

        print("SPICE:", spice_time)
        print("Lambert:", lambert_time)
        print("Plotting:", plot_time)

        plt.show()