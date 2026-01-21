import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import G, c
from scipy.integrate import solve_ivp

def system(s, y, M):
    #Photon equations
    u, dudphi = y
    q = dudphi
    dq = -u +3*M*u**2
    return [q, dq]

#Constants
AU = 1.496e11 #m
M_sun_kg = 1.989e30 #kg
M_geo = G*M_sun_kg/c**2 #m, solar mass in geometrised units
R_sun = 6.96e8  # m
r0 = R_sun*1.1


#Integration parameters
phi_span = [0, 8*np.pi]
phi_eval = np.linspace(phi_span[0], phi_span[1], 800000)
u0 = [1/r0, 0.0]

photon_solution = solve_ivp(
    system, 
    t_span = phi_span, 
    y0 = u0, 
    method = 'RK45', 
    t_eval = phi_eval,
    args = (M_geo,), 
    rtol = 1e-12, 
    atol = 1e-15)

phi = photon_solution.t
u = photon_solution.y[0]
r = 1/u
x = r*np.cos(phi)
y = r*np.sin(phi)
dudphi = photon_solution.y[1]

plt.figure(figsize = (6, 6))
plt.title('Photon Trajectory with GR Effects', fontsize = 20)
plt.plot(x/R_sun, y/R_sun, color = 'magenta', label = 'Photon Trajectory')
plt.scatter(0, 0, marker = 'o', zorder  = 10, s = 100, color = 'orange', label = 'Sun')
plt.scatter(x[0]/R_sun, y[0]/R_sun, zorder = 10, s = 20, color = 'red', label = 'Initial')
plt.scatter(x[-1]/R_sun, y[-1]/R_sun, zorder = 9, s = 60, color = 'blue', label = 'Final')
plt.xlabel('x (solar radii)', fontsize = 18)
plt.ylabel('y (solar radii)', fontsize = 18)
plt.grid()
plt.legend(loc = 'upper right')
plt.show()