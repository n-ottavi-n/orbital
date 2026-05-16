import matplotlib.pyplot as plt
import numpy as np
import planetary_data
import math
from scipy.optimize import fsolve
import re
import matplotlib.animation as animation
import datetime as date
from datetime import datetime
import os
import pandas as pd
import spiceypy as spice

def get_resonant_orbit(body, resonance_ratio, t0, cb=planetary_data.earth):
    '''
    returns the coes of a resonant orbit with a given body and resonance ratio
    @param body: spice body string
    @param resonance_ratio: list of two integers [p, q] where p is the number of orbits of the satellite and q is the number of orbits of the body in the same time period
    @param t0: epoch of date in string format 'Sep 17, 2023, 00:00 UTC'
    @param cb: central body, spice body string
    @return: coes of the resonant orbit
    '''
    mu = cb['mu']
    # get the orbital period of the body
    t0_et = spice.str2et(t0)
    sv, _lighttimes = spice.spkezr(body, t0_et, 'J2000', 'NONE', cb['name'])
    elems=spice.oscelt(sv, t0_et, mu)
    rp, e, i, raan, argp, m0 = elems[0:6]
    mu = elems[-1]
    ra = rp * (1 + e) / (1 - e) # apoapsis distance
    a = (rp + ra) / 2 # semi-major axis
    body_period = 2 * math.pi * math.sqrt((a**3) / mu)
    # calculate the orbital period of the satellite using the resonance ratio
    satellite_period = body_period * (resonance_ratio[1] / resonance_ratio[0])
    # calculate the semi-major axis of the satellite using Kepler's third law
    satellite_a = ((mu * (satellite_period**2)) / (4 * (math.pi**2)))**(1/3)
    satellite_rp = rp
    satellite_ra = 2 * satellite_a - satellite_rp
    satellite_e = (satellite_ra - satellite_rp) / (satellite_ra + satellite_rp)
    satellite_i = i*spice.dpr()
    satellite_raan = raan*spice.dpr()
    satellite_argp = argp*spice.dpr()
    satellite_m0 = 0
    sat_period = 2 * math.pi * math.sqrt((satellite_a**3) / mu)
    print("period ratio:", sat_period/body_period)
    return [satellite_a, satellite_e, satellite_i, satellite_m0, satellite_argp, satellite_raan]
    # return the coes of the resonant orbit

def plot_3d(rs, cb, show_plot=False, save_plot=False, au_units=False):
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_aspect('equal')

    # plot central body

    radius = cb['radius']

    if au_units:
        radius = radius / 1.495978707e8
        rs = rs / 1.495978707e8

    # _u, _v = np.mgrid[0:2 * np.pi:50j, 0:np.pi:50j]
    _u = np.linspace(0, 2 * np.pi, 100)
    _v = np.linspace(0, np.pi, 100)
    _x = radius * np.outer(np.cos(_u), np.sin(_v))
    _y = radius * np.outer(np.sin(_u), np.sin(_v))
    _z = radius * np.outer(np.ones(np.size(_u)), np.cos(_v))
    ax.plot_surface(_x, _y, _z, color='linen', alpha=0.5)

    # plot circular curves over the surface
    theta = np.linspace(0, 2 * np.pi, 100)
    z = np.zeros(100)
    x = radius * np.sin(theta)
    y = radius * np.cos(theta)

    ax.plot(x, y, z, color='black', alpha=0.75)
    ax.plot(z, x, y, color='black', alpha=0.75)

    ## add axis lines
    zeros = np.zeros(1000)
    line = np.linspace(-10, 10, 1000)

    ax.plot(line, zeros, zeros, color='black', alpha=0.75)
    ax.plot(zeros, line, zeros, color='black', alpha=0.75)
    ax.plot(zeros, zeros, line, color='black', alpha=0.75)

    # plot trajectory
    ax.plot(rs[:, 0], rs[:, 1], rs[:, 2], 'g', label='trajectory', zorder=4, linewidth=0.25)
    ax.plot([rs[0, 0]], [rs[0, 1]], [rs[0, 2]], 'wo', label='initial position', zorder=5, linewidth=0.25)
    ax.plot([rs[-1, 0]], [rs[-1, 1]], [rs[-1, 2]], 'yo', label='final position', zorder=6, linewidth=0.25)

    # plot x y z vectors
    l = radius * 2

    max_val = np.max(np.abs(rs))

    ax.set_xlim([-max_val, max_val])
    ax.set_ylim([-max_val, max_val])
    ax.set_zlim([-max_val, max_val])

    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')

    if au_units:
        ax.set_xlabel('X (au)')
        ax.set_ylabel('Y (au)')
        ax.set_zlabel('Z (au)')

    plt.legend()

    if show_plot:
        plt.show()
    if save_plot:
        plt.savefig('title' + '.png')

def plot_n_orbits(rs, step_t, labels, cb, show_plot=False, save_plot=False, au_units=False, equal_axes=True, save_file='matplot003.png'):
    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_aspect('equal')

    # plot central body

    radius=cb['radius']

    if au_units:
        radius =radius /1.495978707e8
        rs =rs /1.495978707e8

    #_u, _v = np.mgrid[0:2 * np.pi:50j, 0:np.pi:50j]
    _u = np.linspace(0, 2 * np.pi, 100)
    _v = np.linspace(0, np.pi, 100)
    _x = radius * np.outer(np.cos(_u), np.sin(_v))
    _y = radius * np.outer(np.sin(_u), np.sin(_v))
    _z = radius * np.outer(np.ones(np.size(_u)), np.cos(_v))
    ax.plot_surface(_x, _y, _z, color='linen', alpha=0.5)

    # plot circular curves over the surface
    theta = np.linspace(0, 2 * np.pi, 100)
    z = np.zeros(100)
    x = radius * np.sin(theta)
    y = radius * np.cos(theta)

    ax.plot(x, y, z, color='black', alpha=0.75)
    ax.plot(z, x, y, color='black', alpha=0.75)

    ## add axis lines
    zeros = np.zeros(1000)
    line = np.linspace(-radius, radius, 1000)

    ax.plot(line, zeros, zeros, color='black', alpha=0.75)
    ax.plot(zeros, line, zeros, color='black', alpha=0.75)
    ax.plot(zeros, zeros, line, color='black', alpha=0.75)

    # plot trajectory
    n=0
    for r in rs:
        ax.plot(r[:, 0], r[:, 1], r[:, 2],  label=labels[n], zorder=4, linewidth=0.5)
        ax.plot([r[0, 0]], [r[0, 1]], [r[0, 2]],'wo', zorder=5, linewidth=0.5)
        ax.plot([r[-1, 0]], [r[-1, 1]], [r[-1, 2]], 'yo', zorder=6, linewidth=0.5)
        n+=1

    max_val = np.max(np.abs(rs))

    ax.set_xlim([-max_val, max_val])
    ax.set_ylim([-max_val, max_val])
    ax.set_zlim([-max_val, max_val])

    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')

    if au_units:
        ax.set_xlabel('X (au)')
        ax.set_ylabel('Y (au)')
        ax.set_zlabel('Z (au)')

    plt.legend()

    #return fig

    if show_plot:
        plt.show()
    if save_plot:
        plt.savefig(save_file+'.png', dpi=400)


def coes2rv(coes, deg=False, mu=planetary_data.earth['mu']):
    a, e, i, m, aop, raan = coes
    if deg:
        i    = np.radians(i)
        m    = np.radians(m)
        aop  = np.radians(aop)
        raan = np.radians(raan)
    E = ecc_anom(m, e)  # eccentric anomaly
    beta   = e / (1 + np.sqrt(1 - e**2))
    nu     = E + 2 * np.arctan((beta * np.sin(E)) / (1 - beta * np.cos(E)))  # true anomaly
    radius = a * (1 - e * np.cos(E))
    h      = np.sqrt(mu * a * (1 - e**2))
    p      = radius * (1 + e * np.cos(nu))
    x      = radius * (np.cos(raan) * np.cos(aop + nu) - np.sin(raan) * np.sin(aop + nu) * np.cos(i))
    y      = radius * (np.sin(raan) * np.cos(aop + nu) + np.cos(raan) * np.sin(aop + nu) * np.cos(i))
    z      = radius * (np.sin(i) * np.sin(aop + nu))
    r      = [x, y, z]
    x_dot  = (((x * h * e) / (radius * p)) * np.sin(nu)) - (h / radius) * (np.cos(raan) * np.sin(aop + nu) + np.sin(raan) * np.cos(aop + nu) * np.cos(i))
    y_dot  = (((y * h * e) / (radius * p)) * np.sin(nu)) - (h / radius) * (np.sin(raan) * np.sin(aop + nu) - np.cos(raan) * np.cos(aop + nu) * np.cos(i))
    z_dot  = (((z * h * e) / (radius * p)) * np.sin(nu)) + (h / radius) * (np.sin(i) * np.cos(aop + nu))
    v      = [x_dot, y_dot, z_dot]
    return r, v

def rv2coes(r,v, deg=False, mu=planetary_data.earth['mu']):
    h=np.cross(r,v) #angular momentum
    e=(np.cross(v,h)/mu)-(r/np.linalg.norm(r)) #eccentricity vector
    n=np.cross(np.array([0,0,1]).transpose(),h) #vector toward ascending node
    norm_n=np.linalg.norm(n)
    nu = math.acos(np.dot(e, r) / (np.linalg.norm(e) * np.linalg.norm(r))) # true anomaly
    if np.dot(r,v)<0:
        nu = (2*math.pi)-nu
    i=math.acos(h[2]/np.linalg.norm(h)) #inclination
    e_norm=np.linalg.norm(e) #eccentricity
    E=2*math.atan2(math.tan(nu/2),math.sqrt((1+e_norm)/(1-e_norm))) #eccentric anomaly
    lan=math.acos(n[0]/norm_n) # longitude of ascending node
    if n[1]<0:
        lan=(2*math.pi)-lan
    aop=math.acos(np.dot(n,e)/(norm_n*e_norm)) #argument of periapsis
    if e[2]<0:
        aop=(2*math.pi)-aop
    M=E-e*math.sin(E) #mean anomaly
    a=1/((2/np.linalg.norm(r))-((np.linalg.norm(v)**2)/mu)) #semi major axis

    if deg:
        nu=math.degrees(nu)
        i=math.degrees(i)
        lan=math.degrees(lan)
        aop=math.degrees(aop)


    return [a,e_norm,i,nu,aop,lan]

def ecc_anom(m,e):
    m = float(m)
    e = float(e)
    res=m
    if e!=0:
        f=lambda x: x-e*np.sin(x)-m
        res=fsolve(f, [50])
        return float(res[0])
    return res

def tle2coes(filename, mu=planetary_data.earth['mu']):
    f = open(filename, "r")
    #get name
    line0=f.readline().strip().split(' ')
    name=' '.join(line0[1:])
    # get elements
    line1=f.readline()
    line2=f.readline().strip().split(' ')
    i=float(line2[3])
    raan=float(line2[4])
    e=float('0.'+line2[5])
    aop=float(line2[6])
    ma=float(line2[8])
    mm=float(line2[10]) #mean motion
    T=(1/mm)*24*3600
    a=((mu*(T**2))/(4*(math.pi**2)))**(1/3.0)
    f.close()
    return [a,e,i,ma,aop,raan, name]

def n_tle2coes(filename, n_objects, t0, mu=planetary_data.earth['mu']):
    f = open(filename, "r")
    #get name
    names=[]
    elements=[]
    n_object=0
    while n_object<n_objects:
        line0 = f.readline().strip().split(' ')
        name = ' '.join(line0[1:]) #get full name

        # get elements
        line1 = f.readline()
        line2 = re.split(r"\s+", f.readline().strip())
        i = float(line2[2])
        raan = float(line2[3])
        e = float('0.' + line2[4])
        aop = float(line2[5])
        ma = float(line2[6]) #mean anomaly at epoch
        mm = float(line2[7])  # mean motion rev/day

        #days between t0 and epoch
        es_epoch=float(line1[3]) #element set epoch
        t0_epoch=get_epoch(t0) # epoch at t0

        ma_t0=(ma+(t0_epoch-es_epoch)*mm*360)%360 #mean anomaly at t0

        T = (1 / mm) * 24 * 3600
        a = ((mu * (T ** 2)) / (4 * (math.pi ** 2))) ** (1 / 3.0)
        elements.append([a,e,i,ma_t0,aop,raan])
        names.append(name)
        n_object+=1
    f.close()
    #print((a * (1 - e))-6378.0, (a * (1 + e))-6378.0)
    return elements, names

def perif2eq(rv, coes):
    return 0

def plot_n_orbits_animate(rs, step_t, labels, cb, show_plot=False, save=False, au_units=False, equal_axes=True, interval=10, frame_step=5, save_file='matplot003.mp4'):
    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_aspect('equal')

    # plot central body

    radius=cb['radius']
    if au_units:
        radius=radius/1.495978707e8
        rs=rs/1.495978707e8

    #_u, _v = np.mgrid[0:2 * np.pi:50j, 0:np.pi:50j]
    _u = np.linspace(0, 2 * np.pi, 100)
    _v = np.linspace(0, np.pi, 100)
    _x = radius * np.outer(np.cos(_u), np.sin(_v))
    _y = radius * np.outer(np.sin(_u), np.sin(_v))
    _z = radius * np.outer(np.ones(np.size(_u)), np.cos(_v))
    ax.plot_surface(_x, _y, _z, color='linen', alpha=0.5)

    # plot circular curves over the surface
    theta = np.linspace(0, 2 * np.pi, 100)
    z = np.zeros(100)
    x = radius * np.sin(theta)
    y = radius * np.cos(theta)

    ax.plot(x, y, z, color='black', alpha=0.75)
    ax.plot(z, x, y, color='black', alpha=0.75)

    ## add axis lines
    zeros = np.zeros(1000)
    line = np.linspace(-radius, radius, 1000)

    ax.plot(line, zeros, zeros, color='black', alpha=0.75)
    ax.plot(zeros, line, zeros, color='black', alpha=0.75)
    ax.plot(zeros, zeros, line, color='black', alpha=0.75)

    def update(num, data, line2_lst, current_poses):
        day = (step_t * num) // (3600*24)
        for i in range(len(data)):
            ax.set_title("time: {} days".format(day))
            line2_lst[i].set_data(data[i][:num, 0], data[i][:num, 1])
            line2_lst[i].set_3d_properties(data[i][:num, 2])
            current_poses[i].set_data([data[i][num, 0]], [data[i][num, 1]])
            current_poses[i].set_3d_properties([data[i][num, 2]])
        if num == N:
            return 0


    data = rs
    N = len(data[0])

    line2_lst = []
    current_poses=[]
    n=0
    for r in data:
        ax.plot([r[0, 0]], [r[0, 1]], [r[0, 2]], 'wo',  zorder=5) # initial position
        r_current, = ax.plot([r[0, 0]], [r[0, 1]], [r[0, 2]], 'yo', zorder=5) #final position
        line2, = ax.plot([r[0, 0]], [r[0, 1]], [r[0, 2]], label=labels[n], zorder=4)  # initialize
        current_poses.append(r_current)
        line2_lst.append(line2)
        n+=1


    plt.legend()

    print("animating...")
    frames = range(0, N, frame_step)
    ani = animation.FuncAnimation(fig, update, frames, fargs=(data, line2_lst, current_poses), interval=interval, blit=False, repeat=True)
    print("done!")


    max_val = np.max(np.abs(rs))

    if equal_axes:
        ax.set_xlim([-max_val, max_val])
        ax.set_ylim([-max_val, max_val])
        ax.set_zlim([-max_val, max_val])

    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')

    if au_units:
        ax.set_xlabel('X (au)')
        ax.set_ylabel('Y (au)')
        ax.set_zlabel('Z (au)')

    plt.legend()

    if show_plot:
        plt.show()
    if save:
        print("saving...")
        ani.save(save_file, writer='pillow', fps=30)
        print("saved!")

def get_epoch(t0):
    '''
    returns the tle formatted epoch of date t0
    @param t0: str format='Sep 17, 2023, 00:00 UTC'
    @return: float time in format= yyddd.fraction of days
    exemple: Jan 14 2022 12:00 UTC = 22014.5
    '''
    start = datetime.strptime(t0, '%b %d, %Y, %H:%M %Z')
    yr = start.year % 100
    day_of_year = start.timetuple().tm_yday  # returns 1 for January 1st

    dt = date.timedelta(hours=start.hour, minutes=start.minute, seconds=start.second)
    secs_per_day = 24 * 60 * 60  # hours * mins * secs
    frac_day = dt.total_seconds() / secs_per_day
    res = yr * 1000 + day_of_year + frac_day
    return res

def plot_pert_coes(coes, ts, labels,  hours=False, days=False):
    print("plotting...")

    xlabel = "seconds"

    if hours:
        ts = ts / 3600
        xlabel = "hours"

    if days:
        ts = ts / (3600 * 24)
        xlabel = "days"

    fig, axs = plt.subplots(2, 3)
    fig.suptitle("Evolution of orbital elements")


    for i in range(len(coes)):
        c=coes[i]
        # plot semi major axis
        axs[0, 0].plot(ts, c[:, 0], label=labels[i])
        axs[0, 0].legend(loc="upper right")
        # plot eccentricity
        axs[0, 1].plot(ts, c[:, 1], label=labels[i])
        axs[0, 1].legend(loc="upper right")
        # plot inclination
        axs[0, 2].plot(ts, c[:, 2], label=labels[i])
        axs[0, 2].legend(loc="upper right")
        # plot true anomaly
        axs[1, 0].plot(ts, c[:, 3], label=labels[i])
        axs[1, 0].legend(loc="upper right")
        # plot argument of periapsis trendline
        #z = np.polyfit(ts.flatten(), c[:, 4], 1)
        #p = np.poly1d(z)
        #axs[1, 1].plot(ts, p(ts), label=labels[i])
        #plot argument of periapsis
        axs[1, 1].plot(ts, c[:, 4], label=labels[i])
        axs[1, 1].legend(loc="upper right")
        # plot longitude of ascending node
        axs[1, 2].plot(ts, c[:, 5], label=labels[i])
        axs[1, 2].legend(loc="upper right")

    axs[0, 0].set_title('semi-major axis vs time')
    axs[0, 0].grid(True)
    axs[0, 0].set_ylabel("a (km)")
    axs[0, 0].set_xlabel(xlabel)

    axs[0, 1].set_title('eccentricity vs time')
    axs[0, 1].grid(True)
    axs[0, 1].set_ylabel("e")
    axs[0, 1].set_xlabel(xlabel)

    axs[0, 2].set_title('inclination vs time')
    axs[0, 2].grid(True)
    axs[0, 2].set_ylabel("i (deg)")
    axs[0, 2].set_xlabel(xlabel)

    axs[1, 0].set_title('true anomaly vs time')
    axs[1, 0].grid(True)
    axs[1, 0].set_ylabel("nu (deg)")
    axs[1, 0].set_xlabel(xlabel)

    axs[1, 1].set_title('argument of periapsis vs time')
    axs[1, 1].grid(True)
    axs[1, 1].set_ylabel("aop (deg)")
    axs[1, 1].set_xlabel(xlabel)

    axs[1, 2].set_title('longitude of ascending node vs time')
    axs[1, 2].grid(True)
    axs[1, 2].set_ylabel("lan (deg)")
    axs[1, 2].set_xlabel(xlabel)

    plt.show()

def get_sats_from_file(sat_names, t0, mu=planetary_data.earth['mu']):
    print(sat_names)
    directory = 'data'
    files=[]
    labels = []
    coes_lst = []
    for filename in os.scandir(directory):
        if filename.is_file():
            f = open(filename, "r")
            for line in f:
                #try to get a name from line
                line0 = f.readline().strip().split(' ')
                name = ' '.join(line0[1:])  # get full name
                if name in sat_names:
                    labels.append(name)
                    #read the next two lines
                    line1 = f.readline()
                    line2 = re.split(r"\s+", f.readline().strip())
                    coes=get_coes_from_tle(line1, line2, t0, mu)
                    coes_lst.append(coes)
            f.close()
            files.append(filename.path)
    return labels, coes_lst

def get_coes_from_tle(line1, line2, t0, mu):
    i = float(line2[2])
    raan = float(line2[3])
    e = float('0.' + line2[4])
    aop = float(line2[5])
    ma = float(line2[6])  # mean anomaly at epoch
    mm = float(line2[7])  # mean motion rev/day

    # days between t0 and epoch
    es_epoch = float(line1[3])  # element set epoch
    t0_epoch = get_epoch(t0)  # epoch at t0

    ma_t0 = (ma + (t0_epoch - es_epoch) * mm * 360) % 360  # mean anomaly at t0

    T = (1 / mm) * 24 * 3600
    a = ((mu * (T ** 2)) / (4 * (math.pi ** 2))) ** (1 / 3.0)
    return [a, e, i, ma_t0, aop, raan]

def read_tle_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Extract satellite name from the first line
    satellite_name = lines[0][2:].strip()

    return satellite_name

def read_tle_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Extract satellite name from the first line
    satellite_names = [line[2:].strip() for line in lines if line.startswith('0')]
    return satellite_names

def get_satellite_names(folder_path):
    satellite_names = []

    # Iterate over all files in the folder
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)

        # Check if the item is a file and not a subdirectory
        if os.path.isfile(file_path):
            satellite_names.extend(read_tle_file(file_path))

    return satellite_names

def create_dataframe_and_write_to_csv(satellite_names):
    # Create a DataFrame with the satellite names
    df = pd.DataFrame({'Satellite Names': satellite_names})

    return df

def load_solar_system_kernels():
    SRC_DIR = os.path.dirname(os.path.abspath(__file__)).replace('\\', '/')
    # Generate metakernel content dynamically
    mk_content = f"""\\begindata

    PATH_VALUES = ('{SRC_DIR}')
    PATH_SYMBOLS = ('SRCDIR')

    KERNELS_TO_LOAD=(
    '$SRCDIR/spice_solar_system/codes_300ast_20100725.bsp',
    '$SRCDIR/spice_solar_system/codes_300ast_20100725.tf',
    '$SRCDIR/spice_lunar/de440.bsp',
    '$SRCDIR/spice_solar_system/latest_leapseconds.tls',
    '$SRCDIR/spice_solar_system/mars_iau2000_v1.tpc',
    '$SRCDIR/spice_solar_system/mar099s.bsp',
    '$SRCDIR/spice_solar_system/jup348.bsp',
    '$SRCDIR/spice_solar_system/jup365.bsp',
    '$SRCDIR/spice_solar_system/sat459.bsp',
    '$SRCDIR/spice_solar_system/54404085.bsp',
    '$SRCDIR/spice_solar_system/pck00011.tpc',
    '$SRCDIR/spice_solar_system/naif0012.tls'
    )
    \\begintext
    """

    # Write temporary metakernel and load it
    mk_path = os.path.join(SRC_DIR, 'temp_kernel.tm')
    with open(mk_path, 'w') as f:
        f.write(mk_content)

    spice.kclear()
    spice.furnsh(mk_path)

    # Clean up temp file
    os.remove(mk_path)

def load_earth_moon_kernels():
    SRC_DIR = os.path.dirname(os.path.abspath(__file__)).replace('\\', '/')
    # Generate metakernel content dynamically
    mk_content = f"""\\begindata

    PATH_VALUES = ('{SRC_DIR}')
    PATH_SYMBOLS = ('SRCDIR')

    KERNELS_TO_LOAD=(
    '$SRCDIR/spice_lunar/de440.bsp',
    '$SRCDIR/spice_lunar/moon_de440_220930.tf',
    '$SRCDIR/spice_lunar/moon_pa_de440_200625.bpc',
    '$SRCDIR/spice_lunar/clem_nrl.bsp',
    '$SRCDIR/spice_lunar/clem_v20.tf'
    '$SRCDIR/spice_solar_system/latest_leapseconds.tls',
    '$SRCDIR/spice_solar_system/naif0012.tls'
    )
    \\begintext
    """

    # Write temporary metakernel and load it
    mk_path = os.path.join(SRC_DIR, 'temp_kernel.tm')
    with open(mk_path, 'w') as f:
        f.write(mk_content)

    spice.kclear()
    spice.furnsh(mk_path)

    # Clean up temp file
    os.remove(mk_path)
