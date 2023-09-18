from matplotlib import pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

radius=1

_u = np.linspace(0, 2 * np.pi, 100)
_v = np.linspace(0, np.pi, 100)
_x = radius * np.outer(np.cos(_u) , np.sin(_v))
_y = radius * np.outer(np.sin(_u) , np.sin(_v))
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

def gen(n):
    phi = 0
    while phi < 2*np.pi:
        yield np.array([np.cos(phi), np.sin(phi), phi])
        phi += 2*np.pi/n

def update(num, data, line):
    print("updating...")
    print(data[:2, :num])
    print(data[:2, :num].shape)
    line.set_data(data[:2, :num])
    line.set_3d_properties(data[2, :num])

N = 100
data = np.array(list(gen(N))).T
line, = ax.plot(data[0, 0:1], data[1, 0:1], data[2, 0:1])


# Setting the axes properties
ax.set_xlim3d([-1.0, 1.0])
ax.set_xlabel('X')

ax.set_ylim3d([-1.0, 1.0])
ax.set_ylabel('Y')

ax.set_zlim3d([0.0, 10.0])
ax.set_zlabel('Z')

ani = animation.FuncAnimation(fig, update, N, fargs=(data, line), interval=10000/N, blit=False)
#ani.save('matplot003.gif', writer='imagemagick')
plt.show()