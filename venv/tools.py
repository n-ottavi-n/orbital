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
    ax.plot(rs[:, 0], rs[:, 1], rs[:, 2], 'g', label='trajectory', zorder=4)
    ax.plot([rs[0, 0]], [rs[0, 1]], [rs[0, 2]], 'wo', label='initial position', zorder=5)
    ax.plot([rs[-1, 0]], [rs[-1, 1]], [rs[-1, 2]], 'yo', label='final position', zorder=6)

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
        plt.savefig(title + '.png')

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
        ax.plot(r[:, 0], r[:, 1], r[:, 2],  label=labels[n], zorder=4)
        ax.plot([r[0, 0]], [r[0, 1]], [r[0, 2]],'wo', zorder=5)
        ax.plot([r[-1, 0]], [r[-1, 1]], [r[-1, 2]], 'yo', zorder=6)
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

    if show_plot:
        plt.show()
    if save_plot:
        plt.savefig(save_file+'.png', dpi=400)

def coes2rv(coes, deg=False, mu=planetary_data.earth['mu']):
    '''
    :param coes: a,e,i,m,aop,raan
    :param deg: true if degrees
    :param mu: of central body
    :return: [x,y,z],[vx,vy,vz]
    '''
    a,e,i,m,aop,raan=coes
    E=ecc_anom(m,e) #eccentric anomaly
    if deg:
        i = math.radians(i)
        m = math.radians(m)
        aop = math.radians(aop)
        raan = math.radians(raan)


    beta=e/(1+math.sqrt(1-e**2))
    nu=E+2*math.atan((beta*math.sin(E))/(1-beta*math.cos(E))) #true anomaly
    radius=a*(1-e*math.cos(E)) #radius
    h=math.sqrt(mu*a*(1-e**2)) #specific angular momentum
    p=radius*(1+e*math.cos(nu)) #semi-latus rectum :)

    x = radius * (math.cos(raan) * math.cos(aop + nu) - math.sin(raan) * math.sin(aop + nu) * math.cos(i))
    y = radius * (math.sin(raan) * math.cos(aop + nu) + math.cos(raan) * math.sin(aop + nu) * math.cos(i))
    z = radius * (math.sin(i) * math.sin(aop + nu))

    r = [x,y,z]

    x_dot = (((x*h*e)/(radius*p))*math.sin(nu))-(h/radius)*(math.cos(raan) * math.sin(aop + nu) + math.sin(raan) * math.cos(aop + nu) * math.cos(i))
    y_dot = (((y * h * e) / (radius * p)) * math.sin(nu)) - (h / radius) * (math.sin(raan) * math.sin(aop + nu) - math.cos(raan) * math.cos(aop + nu) * math.cos(i))
    z_dot = (((z * h * e) / (radius * p)) * math.sin(nu))+(h/radius)*(math.sin(i) * math.cos(aop + nu))

    v = [x_dot,y_dot,z_dot]

    return r,v

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
    res=m
    if e!=0:
        f=lambda x: x-e*math.sin(x)-m
        res=fsolve(f, [50])
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
        line2 = re.split("\s+", f.readline().strip())
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

def plot_n_orbits_animate(rs, step_t, labels, cb, show_plot=False, save=False, au_units=False, equal_axes=True, interval=0.1, save_file='matplot003.gif'):
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
        hr=(step_t*num)//3600
        for i in range(len(data)):
            ax.set_title("time: {} hrs".format(hr))
            line2_lst[i].set_data(data[i][:num, :2].transpose())
            line2_lst[i].set_3d_properties(data[i][:num, 2])
            # plot current position
            current_poses[i].set_data(data[i][num, :2].transpose())
            current_poses[i].set_3d_properties(data[i][num, 2])
        if num==N:
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
    ani = animation.FuncAnimation(fig, update, N, fargs=(data, line2_lst, current_poses), interval=interval, blit=False, repeat=True)
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
        ani.save(save_file, writer='imagemagick', fps=30)
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
    directory = '../data'
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
                    line2 = re.split("\s+", f.readline().strip())
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