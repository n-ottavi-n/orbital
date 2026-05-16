import numpy as np
import matplotlib.pyplot as plt
import planetary_data
from scipy.integrate import ode
import tools as t
import planetary_data as pd
import matplotlib.animation as animation
import spiceypy as spice
from plotly.subplots import make_subplots
import plotly.graph_objects as go

class Propagator:

    def __init__(self, state0, n_steps, start_date='Sep 16, 2023, 00:00 UTC', end_date='Sep 17, 2023, 00:00 UTC',spacecraft_data={}, coes=False, deg=False, cb=pd.earth, perts=[], integrator='lsoda'):
        '''
        @param state0: if cartesian: [x0,y0,z0,vx0,vy0,vz0]  if keplerian: [a,e,i,m,aop,raan]
        @param n_steps: number of timesteps for integration
        @param start_date: format:  '%b %d, %Y, %H:%M %Z'
        @param end_date: format:  '%b %d, %Y, %H:%M %Z'
        @param coes: bool: True if state0 is keplerian elements
        @param deg: bool: True if input keplerian elements in degrees
        @param cb: central body
        @param perts: perturbations to account for / ['j2','srp', [pd.mars,pd.venus,pd.moon,.....]]
        '''
        if coes:
            self.r0,self.v0=t.coes2rv(state0, deg=deg, mu=cb['mu'])
        else:
            self.r0=state0[:3]
            self.v0=state0[3:]

        #self.tspan = tspan
        self.n_steps=n_steps
        self.cb = cb
        self.perts=perts
        self.integrator=integrator

        #make a vector of times from start to end dates
        etOne = spice.str2et(start_date)
        etTwo = spice.str2et(end_date)
        self.times = [x * (etTwo - etOne) / n_steps + etOne for x in range(n_steps)]
        self.dt=self.times[1]-self.times[0]

        self.j2=False
        if "j2" in perts:
            self.j2=True

        self.srp=False
        if "srp" in perts:
            self.vectors_from_sun, lightTimes =spice.spkpos(self.cb['name'], self.times, 'J2000', 'NONE', "SUN")
            self.srp=True

        
        self.bodies=[]
        #check if list of perturbing bodies is present then store it.
        for pert in perts:
            if type(pert) is list: 
                self.pert_bodies=True 
                self.bodies=pert # list of pd objects  if type is list then it is a list of perturbing bodies, else it is a perturbation type like j2 or srp

        self.body_names=[]
        self.r_cb2nb_lst=[]

        if self.bodies: 
            for body in self.bodies: #get vector from central body to perturbing body for each perturbing body and store it for later use in acceleration calculation
                self.body_names.append(body['name']) #for future debugging
                r_cb2nb, lightTimes = spice.spkpos(body['name'], self.times, 'J2000', 'NONE', cb['name']) #vector from central body to perturbing body
                self.r_cb2nb_lst.append(r_cb2nb)

        self.sc_data=spacecraft_data


    def propagate(self):

        self.ts = np.zeros((self.n_steps,1))
        self.ys = np.zeros((self.n_steps,6))
        self.y0 = self.r0+self.v0
        self.ts[0]=0
        self.ys[0]=self.y0
        self.steps=1

        self.solver=ode(self.diffy_q)
        self.solver.set_integrator(self.integrator)
        self.solver.set_initial_value(self.y0, 0)

        while self.solver.successful() and self.steps < self.n_steps:
            self.solver.integrate((self.solver.t + self.dt))
            self.ts[self.steps]=self.solver.t
            new_ys=self.solver.y
            self.ys[self.steps]=new_ys
            self.steps += 1

        self.rs=self.ys[:,:3]
        self.vs = self.ys[:,3:]



    def diffy_q(self, t, y):
        rx,ry,rz,vx,vy,vz=y
        r=np.array([rx,ry,rz])
        v=np.array([vx,vy,vz])

        norm_r=np.linalg.norm(r)

        ax,ay,az=-r*self.cb['mu']/norm_r**3

        a=[ax,ay,az]

        if self.j2:
            j2=self.cb['j2']

            z2 = r[2] ** 2
            r2 = norm_r ** 2
            tx = r[0] / norm_r * (5 * z2 / r2 - 1)
            ty = r[1] / norm_r * (5 * z2 / r2 - 1)
            tz = r[2] / norm_r * (5 * z2 / r2 - 3)
            a += 1.5 * self.cb['j2'] * self.cb['mu'] * self.cb['radius'] ** 2 / norm_r ** 4 * np.array([tx, ty, tz])

        if self.srp:
            pressure=4.5*(10**(-6)) # sun pressure at earth in Pa
            force=pressure*self.sc_data['area']
            acc=force/self.sc_data['mass']
            d=self.vectors_from_sun[self.steps, :]+r #vector from sun to spacecraft
            d=d/np.linalg.norm(d)
            a+=acc*d
            #source: https://ntrs.nasa.gov/api/citations/20080012725/downloads/20080012725.pdf

        self.pert_bodies=True


        if self.pert_bodies:
            for i in range(len(self.bodies)):
                mu=self.bodies[i]['mu']  #gravitational parameter of perturbing body
                r_sat2body = self.r_cb2nb_lst[i][self.steps, :] - r #vector from satellite to perturbing body
                a += mu * ((r_sat2body / (np.linalg.norm(r_sat2body) ** 3)) - (self.r_cb2nb_lst[i][self.steps, :] / (np.linalg.norm(self.r_cb2nb_lst[i][self.steps, :]) ** 3))) #Cowell's method

        return [vx, vy, vz, a[0], a[1], a[2]]

    def calculate_coes(self):
        print("calculating...")
        coes=[]
        for time in range(len(self.rs)):
            coes.append(t.rv2coes(self.rs[time],self.vs[time], deg=True, mu=self.cb['mu']))
        self.coes=np.array(coes)

    def plot_coes_plotly(self, hours=False, days=False, title=None, save_html=None):
        """
        Plot orbital elements evolution in Plotly.
        Includes linear trend for argument of periapsis (apsidal rotation)
        and longitude of ascending node (nodal regression).
        """
        ts     = self.ts.flatten()
        xlabel = "Time (seconds)"

        if hours:
            ts     = ts / 3600
            xlabel = "Time (hours)"
        elif days:
            ts     = ts / (3600 * 24)
            xlabel = "Time (days)"

        # --------------------------------------------------
        # linear trends for aop and lan
        # --------------------------------------------------
        z_aop  = np.polyfit(ts, self.coes[:, 4], 1)
        p_aop  = np.poly1d(z_aop)
        aop_trend = p_aop(ts)

        z_lan  = np.polyfit(ts, self.coes[:, 5], 1)
        p_lan  = np.poly1d(z_lan)
        lan_trend = p_lan(ts)

        # trend rates
        aop_rate = z_aop[0]   # deg per time unit
        lan_rate = z_lan[0]

        # --------------------------------------------------
        # subplot layout — 2 rows x 3 cols
        # --------------------------------------------------
        fig = make_subplots(
            rows=2, cols=3,
            subplot_titles=(
                "Semi-major axis",
                "Eccentricity",
                "Inclination",
                "True anomaly",
                f"Apsidal rotation  ({aop_rate:+.4f} deg/{xlabel.split()[1]})",
                f"Nodal regression  ({lan_rate:+.4f} deg/{xlabel.split()[1]})",
            ),
            vertical_spacing=0.14,
            horizontal_spacing=0.08,
        )

        # common line style
        line = dict(color='royalblue', width=1.5)
        trend_line = dict(color='tomato', width=1.5, dash='dash')

        # --------------------------------------------------
        # semi-major axis
        # --------------------------------------------------
        fig.add_trace(go.Scatter(
            x=ts, y=self.coes[:, 0],
            mode='lines', line=line,
            name='a', showlegend=False,
            hovertemplate=f"{xlabel.split()[1]}: %{{x:.2f}}<br>a: %{{y:.2f}} km<extra></extra>",
        ), row=1, col=1)

        # --------------------------------------------------
        # eccentricity
        # --------------------------------------------------
        fig.add_trace(go.Scatter(
            x=ts, y=self.coes[:, 1],
            mode='lines', line=line,
            name='e', showlegend=False,
            hovertemplate=f"{xlabel.split()[1]}: %{{x:.2f}}<br>e: %{{y:.6f}}<extra></extra>",
        ), row=1, col=2)

        # --------------------------------------------------
        # inclination
        # --------------------------------------------------
        fig.add_trace(go.Scatter(
            x=ts, y=self.coes[:, 2],
            mode='lines', line=line,
            name='i', showlegend=False,
            hovertemplate=f"{xlabel.split()[1]}: %{{x:.2f}}<br>i: %{{y:.4f}} deg<extra></extra>",
        ), row=1, col=3)

        # --------------------------------------------------
        # true anomaly
        # --------------------------------------------------
        fig.add_trace(go.Scatter(
            x=ts, y=self.coes[:, 3],
            mode='lines', line=line,
            name='nu', showlegend=False,
            hovertemplate=f"{xlabel.split()[1]}: %{{x:.2f}}<br>ν: %{{y:.2f}} deg<extra></extra>",
        ), row=2, col=1)

        # --------------------------------------------------
        # argument of periapsis — raw + trend
        # --------------------------------------------------
        fig.add_trace(go.Scatter(
            x=ts, y=self.coes[:, 4],
            mode='lines',
            line=dict(color='royalblue', width=1, dash='dot'),
            name='aop raw', showlegend=True,
            opacity=0.4,
            hovertemplate=f"{xlabel.split()[1]}: %{{x:.2f}}<br>ω: %{{y:.2f}} deg<extra></extra>",
        ), row=2, col=2)

        fig.add_trace(go.Scatter(
            x=ts, y=aop_trend,
            mode='lines', line=trend_line,
            name='aop trend', showlegend=True,
            hovertemplate=f"{xlabel.split()[1]}: %{{x:.2f}}<br>ω trend: %{{y:.2f}} deg<extra></extra>",
        ), row=2, col=2)

        # --------------------------------------------------
        # longitude of ascending node — raw + trend
        # --------------------------------------------------
        fig.add_trace(go.Scatter(
            x=ts, y=self.coes[:, 5],
            mode='lines',
            line=dict(color='royalblue', width=1, dash='dot'),
            name='lan raw', showlegend=True,
            opacity=0.4,
            hovertemplate=f"{xlabel.split()[1]}: %{{x:.2f}}<br>Ω: %{{y:.2f}} deg<extra></extra>",
        ), row=2, col=3)

        fig.add_trace(go.Scatter(
            x=ts, y=lan_trend,
            mode='lines', line=trend_line,
            name='lan trend', showlegend=True,
            hovertemplate=f"{xlabel.split()[1]}: %{{x:.2f}}<br>Ω trend: %{{y:.2f}} deg<extra></extra>",
        ), row=2, col=3)

        # --------------------------------------------------
        # axis labels
        # --------------------------------------------------
        ylabels = {
            (1,1): "a (km)",
            (1,2): "e",
            (1,3): "i (deg)",
            (2,1): "ν (deg)",
            (2,2): "ω (deg)",
            (2,3): "Ω (deg)",
        }

        for (row, col), ylabel in ylabels.items():
            fig.update_yaxes(
                title=dict(text=ylabel, font=dict(color='black')),
                tickfont=dict(color='black'),
                gridcolor='lightgrey',
                showgrid=True,
                zeroline=False,
                row=row, col=col,
            )
            fig.update_xaxes(
                title=dict(text=xlabel, font=dict(color='black')),
                tickfont=dict(color='black'),
                gridcolor='lightgrey',
                showgrid=True,
                zeroline=False,
                row=row, col=col,
            )

        # --------------------------------------------------
        # layout
        # --------------------------------------------------
        plot_title = title or "Orbital Elements Evolution"

        fig.update_layout(
            title=dict(
                text=plot_title,
                font=dict(size=14, color='black'),
                x=0.5, xanchor='center',
            ),
            paper_bgcolor='white',
            plot_bgcolor='white',
            font=dict(color='black'),
            width=1200,
            height=700,
            legend=dict(
                x=1.02, y=1,
                bordercolor='lightgrey',
                borderwidth=1,
                font=dict(color='black', size=9),
            ),
        )

        if save_html:
            fig.write_html(save_html)
            print(f"Saved to {save_html}")

        fig.show()
        return fig
    
    def plot_distance_to_body(self, body_name, hours=False, days=False,
                           title=None, save_html=None):
        """
        Plot distance from spacecraft to a given body over time.

        Parameters
        ----------
        body_name : str
            SPICE name of the body (e.g. "EUROPA", "JUPITER")
        hours     : bool — plot time in hours
        days      : bool — plot time in days
        title     : str — plot title (optional)
        save_html : str — path to save HTML (optional)
        """
        import plotly.graph_objects as go

        ts     = self.ts.flatten()
        xlabel = "Time (seconds)"

        if hours:
            ts     = ts / 3600
            xlabel = "Time (hours)"
        elif days:
            ts     = ts / (3600 * 24)
            xlabel = "Time (days)"

        # --------------------------------------------------
        # get body positions from SPICE at each timestep
        # --------------------------------------------------
        body_radius = None
        try:
            radii       = spice.bodvrd(body_name, "RADII", 3)[1]
            body_radius = float(np.mean(radii))
        except Exception:
            pass

        r_body = np.array([
            spice.spkpos(body_name, float(t),
                        'J2000', 'NONE',
                        self.cb['name'].upper())[0]
            for t in self.times
        ])

        # --------------------------------------------------
        # distance and altitude
        # --------------------------------------------------
        rel      = self.rs - r_body
        dist_km  = np.linalg.norm(rel, axis=1)

        if body_radius is not None:
            alt_km = dist_km - body_radius
        else:
            alt_km = None

        # --------------------------------------------------
        # figure
        # --------------------------------------------------
        fig = go.Figure()

        # altitude trace (if body radius available)
        if alt_km is not None:
            fig.add_trace(go.Scatter(
                x=ts,
                y=alt_km,
                mode='lines',
                line=dict(color='royalblue', width=1.5),
                name=f"altitude above {body_name.capitalize()}",
                hovertemplate=(
                    f"{xlabel.split()[1]}: %{{x:.2f}}<br>"
                    f"altitude: %{{y:.2f}} km"
                    "<extra></extra>"
                ),
            ))

            # surface line
            fig.add_hline(
                y=0,
                line=dict(color='grey', width=1, dash='dot'),
                annotation_text=f"{body_name.capitalize()} surface",
                annotation_position="bottom right",
                annotation_font=dict(color='grey', size=9),
            )

            # body radius line
            fig.add_hline(
                y=body_radius,
                line=dict(color='steelblue', width=1, dash='dot'),
                annotation_text=f"R = {body_radius:.0f} km",
                annotation_position="top right",
                annotation_font=dict(color='steelblue', size=9),
            )

        # --------------------------------------------------
        # closest approach annotation
        # --------------------------------------------------
        i_min    = np.argmin(dist_km)
        t_min    = ts[i_min]
        d_min    = dist_km[i_min]

        fig.add_trace(go.Scatter(
            x=[t_min], y=[d_min],
            mode='markers+text',
            marker=dict(color='tomato', size=10, symbol='star'),
            text=[f"CA: {d_min:.0f} km"],
            textposition='top right',
            textfont=dict(color='tomato', size=10),
            name=f"closest approach",
            showlegend=True,
            hovertemplate=(
                f"closest approach<br>"
                f"{xlabel.split()[1]}: {t_min:.2f}<br>"
                f"distance: {d_min:.2f} km"
                + (f"<br>altitude: {d_min - body_radius:.2f} km"
                if body_radius else "")
                + "<extra></extra>"
            ),
        ))

        # --------------------------------------------------
        # layout
        # --------------------------------------------------
        plot_title = title or (
            f"Distance to {body_name.capitalize()} over time"
        )

        fig.update_layout(
            title=dict(
                text=plot_title,
                font=dict(size=14, color='black'),
                x=0.5, xanchor='center',
            ),
            xaxis=dict(
                title=dict(text=xlabel, font=dict(color='black')),
                tickfont=dict(color='black'),
                gridcolor='lightgrey',
                showgrid=True,
                zeroline=False,
            ),
            yaxis=dict(
                title=dict(text="Distance (km)", font=dict(color='black')),
                tickfont=dict(color='black'),
                gridcolor='lightgrey',
                showgrid=True,
                zeroline=False,
            ),
            paper_bgcolor='white',
            plot_bgcolor='white',
            font=dict(color='black'),
            width=1000,
            height=500,
            legend=dict(
                x=1.02, y=1,
                bordercolor='lightgrey',
                borderwidth=1,
                font=dict(color='black', size=10),
            ),
        )

        if save_html:
            fig.write_html(save_html)
            print(f"Saved to {save_html}")

        fig.show()
        return fig

    def plot_coes(self, hours=False, days=False):
        print("plotting...")

        ts=self.ts
        xlabel = "seconds"

        if hours:
            ts=ts/3600
            xlabel = "hours"

        if days:
            ts=ts/(3600*24)
            xlabel = "days"

        fig,axs=plt.subplots(2,3)

        # plot semi major axis
        axs[0,0].plot(ts, self.coes[:,0])
        axs[0,0].set_title('semi-major axis vs time')
        axs[0, 0].grid(True)
        axs[0, 0].set_ylabel("a (km)")
        axs[0, 0].set_xlabel(xlabel)

        # plot eccentricity
        axs[0, 1].plot(ts, self.coes[:, 1])
        axs[0, 1].set_title('eccentricity vs time')
        axs[0, 1].grid(True)
        axs[0, 1].set_ylabel("e")
        axs[0, 1].set_xlabel(xlabel)

        # plot inclination
        axs[0, 2].plot(ts, self.coes[:, 2])
        axs[0, 2].set_title('inclination vs time')
        axs[0, 2].grid(True)
        axs[0, 2].set_ylabel("i (deg)")
        axs[0, 2].set_xlabel(xlabel)

        # plot true anomaly
        axs[1, 0].plot(ts, self.coes[:, 3])
        axs[1, 0].set_title('true anomaly vs time')
        axs[1, 0].grid(True)
        axs[1, 0].set_ylabel("nu (deg)")
        axs[1, 0].set_xlabel(xlabel)

        # plot argument of periapsis
        #axs[1, 1].plot(ts, self.coes[:, 4])
        z = np.polyfit(ts.flatten(), self.coes[:, 4] , 1)
        p = np.poly1d(z)
        axs[1, 1].plot(ts, p(ts) )
        axs[1, 1].set_title('apsidal rotation')
        axs[1, 1].grid(True)
        axs[1, 1].set_ylabel("aop (deg)")
        axs[1, 1].set_xlabel(xlabel)

        # plot longitude of ascending node
        axs[1, 2].plot(ts, self.coes[:, 5])
        axs[1, 2].set_title('nodal regression')
        axs[1, 2].grid(True)
        axs[1, 2].set_ylabel("lan (deg)")
        axs[1, 2].set_xlabel(xlabel)

        plt.show()