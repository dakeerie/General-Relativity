import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from einsteinpy.geodesic import Timelike
from scipy.constants import G, c
from scipy.integrate import solve_ivp
from scipy.signal import argrelextrema

def system(s, y, M, h):
    #Matter particle equations
    u, dudphi = y
    q = dudphi
    dq = -u + M/h**2 +3*M*u**2
    return [q, dq]

def specific_am(M, a, e):
    #Newtonian specific angular momentum calculation
    #Works fine for Mercury as r >> 2M for Mercury-Sun system 
    h_sq = M*a*(1-e**2)
    return np.sqrt(h_sq)

#Constants
AU = 1.496e11 #m
M_sun_kg = 1.989e30  # kg
M_sol = G*M_sun_kg/c**2 #m, solar mass in geometrised units
e_merc = 0.2056 #CHECK SOURCES
a_merc = 0.387*AU #m CHECK SOURCES
h_merc = specific_am(M_sol, a_merc, e_merc)
peri_merc = a_merc*(1-e_merc)

#Integration parameters
phi_span = [0, 1000*np.pi]
phi_eval = np.linspace(phi_span[0], phi_span[1], 50000)
u0 = [1/peri_merc, 0] #start at perihelion where du/dphi = 0

#Solve
mercury_solution = solve_ivp(
    system, 
    t_span = phi_span, 
    y0 = u0, 
    method = 'RK45', 
    t_eval = phi_eval, 
    args = (M_sol, h_merc), 
    rtol = 1e-12, 
    atol = 1e-15)

phi = mercury_solution.t
u = mercury_solution.y[0]
r = 1/u
x = r*np.cos(phi)
y = r*np.sin(phi)
dudphi = mercury_solution.y[1]

peri_idx = argrelextrema(r, np.less)
print('Local minima indices:', peri_idx)
# print('Number of local minima:', peri_idx.shape()[1])
print(type(peri_idx))
print(peri_idx[0][0])
print(r[0], r[peri_idx[0][0]], r[peri_idx[0][-1]])

plt.figure(figsize = (6, 6))
plt.title('Mercury Orbit with GR Effects', fontsize = 20)
plt.plot(x/AU, y/AU, color = 'gray', label = 'Mercury Trajectory')
plt.scatter(0, 0, marker = 'o', zorder  = 10, s = 100, color = 'orange', label = 'Sun')
plt.scatter(x[peri_idx[0][0]]/AU, y[peri_idx[0][0]]/AU, zorder = 10, s = 60, color = 'black', label = 'First index')
plt.scatter(x[peri_idx[0][-1]]/AU, y[peri_idx[0][-1]]/AU, zorder = 10, s = 60, color = 'green', label = 'Last index')
plt.scatter(x[0]/AU, y[0]/AU, zorder = 10, s = 20, color = 'red', label = 'Initial')
# plt.scatter(x[-1]/AU, y[-1]/AU, zorder = 9, s = 60, color = 'blue', label = 'Final')
plt.xlabel('x (AU)', fontsize = 18)
plt.ylabel('y (AU)', fontsize = 18)
plt.grid()
plt.legend(loc = 'upper right')
plt.show()

plt.figure(figsize = (10, 5))
plt.plot(phi, r)
plt.show()