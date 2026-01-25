import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import G, c

def dudphi(phi, u, w):
    return w 

def dwdphi(phi, u, w, M: float, h: float, b):
    """b = 1 for matter, b = 0 for light"""
    return -u + b*M/h**2 + 3*M*u**2

def RK4_2nd(f, g, IC, x_final, n):
    x = [IC[0]]
    y = [IC[1]]
    dy = [IC[2]]
    dx = (x_final - x[0])/n
    dx_2 = dx/2
    for i in range(n):
        xi, yi, dyi = x[i], y[i], dy[i]
        k0 = dx*f(xi, yi, dyi)
        l0 = dx*g(x[i], y[i], dy[i])
        k1 = dx*f(x[i] + dx_2, y[i] + k0/2, dy[i] + l0/2)
        l1 = dx*g(x[i] + dx_2, y[i] + k0/2, dy[i] + l0/2)
        k2 = dx*f(x[i] + dx_2, y[i] + k1/2, dy[i] + l1/2)
        l2 = dx*g(x[i] + dx_2, y[i] + k1/2, dy[i] + l1/2)
        k3 = dx*f(x[i] + dx, y[i] + k2, dy[i] + l2)
        l3 = dx*g(x[i] + dx, y[i] + k2, dy[i] + l2)
        x.append(x[i] + dx)
        y.append(y[i] + (k0 + 2*k1 + 2*k2 + k3)/6)
        dy.append(dy[i] + (l0 + 2*l1 + 2*l2 + l3)/6)
    return x, y, dy

#Parameters
M = 1.0
r0 = 20
u0 = 1/r0
#Circular orbit specific momentum
h = np.sqrt(M/(u0 - 3*M*u0**2))
b = 1.0
orbits = 10

wrapped_dwdphi = lambda phi, u, w: dwdphi(phi, u, w, M, h, b)
initial_conditions = [0.0, u0 - 0.01, 0.0]
phi, u, du = RK4_2nd(dudphi, wrapped_dwdphi, initial_conditions, x_final = 2*orbits*np.pi, n = 5000)
phi = np.array(phi)
u = np.array(u)
du = np.array(du)

r = 1/u
x = r*np.cos(phi)
y = r*np.sin(phi)

plt.figure(figsize = [6,6])
plt.plot(x, y, color = 'green', label = 'Matter particle trajectory')
plt.scatter(x[0], y[0], color = 'red', label = 'Initial position', zorder = 10)
plt.scatter(x[-1], y[-1], color = 'blue', label = 'Final position', zorder = 11)
plt.plot(0, 0, 'ko', label = 'Central Mass at (0,0)')
plt.xlabel('x', fontsize = 20)
plt.ylabel('y', fontsize = 20)
plt.axis("equal")
plt.legend()
plt.grid()
plt.show()