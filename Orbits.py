import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from einsteinpy.geodesic import Timelike
from scipy.constants import G, c
from scipy.integrate import solve_ivp

def system(s, y, M, h, object):

    if object == 'matter':
        b = 1
    elif object == 'photon':
        b = 0
    else:
        raise ValueError("Object specification must be either 'matter' or 'photon'.")
    u, dudphi = y
    q = dudphi
    dq = -u + b*M/h**2 +3*M*u**2
    return [q, dq]

u0 = 1/10
du0 = 0
y0 = [u0, du0]

phi_span = [0, 10*np.pi]
phi_eval = np.linspace(phi_span[0], phi_span[1], 2000)

sol = solve_ivp(system, phi_span, y0, 'DOP853', t_eval = phi_eval, args= (1, 4, 'matter'))

phi = sol.t
r = 1/sol.y[0]
dudphi = sol.y[1]

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