import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import G, c
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks

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
M_geo = G*M_sun_kg/c**2 #m, solar mass in geometrised units
e_merc = 0.2056 #CHECK SOURCES
a_merc = 0.387*AU #m CHECK SOURCES
h_merc = specific_am(M_geo, a_merc, e_merc)
peri_merc = a_merc*(1-e_merc)

#Integration parameters
phi_span = [0, 8*np.pi]
phi_eval = np.linspace(phi_span[0], phi_span[1], 800000)
u0 = [1/peri_merc, 0] #start at perihelion where du/dphi = 0

#Solve
mercury_solution = solve_ivp(
    system, 
    t_span = phi_span, 
    y0 = u0, 
    method = 'RK45', 
    t_eval = phi_eval, 
    dense_output = True,
    args = (M_geo, h_merc), 
    rtol = 1e-12, 
    atol = 1e-15)

from scipy.optimize import brentq

sol = mercury_solution  # your solve_ivp result (dense_output=True)
sol_cont = sol.sol      # callable: sol_cont(phi) -> y (2-vector)

def dudphi_at(phi):
    return sol_cont(phi)[1]   # y[1] is dudphi

peri_phi = []
for k in range(1, 4):
    guess = k * 2*np.pi
    left = guess - 5e-7
    right = guess + 5e-7

    # Check that dudphi changes sign inside the bracket
    if np.sign(dudphi_at(left)) == np.sign(dudphi_at(right)):
        # If no sign change, slightly shift bracket
        left -= 1e-7
        right += 1e-7

    phi_root = brentq(dudphi_at, left, right, xtol=1e-14, rtol=1e-12)
    peri_phi.append(phi_root)



# precession per orbit (radians)
delta_phi_meas = abs(peri_phi[1] - peri_phi[0] - 2*np.pi)
delta_phi_arcsec = delta_phi_meas*180*3600/np.pi
orbits_per_century = 100*365.25/87.969
delta_phi_per_century = delta_phi_arcsec*orbits_per_century
print(delta_phi_per_century)


phi = mercury_solution.t
u = mercury_solution.y[0]
r = 1/u
x = r*np.cos(phi)
y = r*np.sin(phi)
dudphi = mercury_solution.y[1]

#Determination of stationary points
minima_indices = find_peaks(1/r)[0]
stationary_phi = []
for i in minima_indices:
    angle = phi[i]
    stationary_phi.append(angle)

#Precession calculation
precession_per_orbit_arcsec = np.abs(stationary_phi[1] - stationary_phi[0] - 2*np.pi)*(180*3600/np.pi)
precession_per_century_arcsec = precession_per_orbit_arcsec*100*365.25/87.969
print(precession_per_century_arcsec)

plt.figure(figsize = (6, 6))
plt.title('Mercury Orbit with GR Effects', fontsize = 20)
plt.plot(x/AU, y/AU, color = 'gray', label = 'Mercury Trajectory')
plt.scatter(0, 0, marker = 'o', zorder  = 10, s = 100, color = 'orange', label = 'Sun')
plt.scatter(x[0]/AU, y[0]/AU, zorder = 10, s = 20, color = 'red', label = 'Initial')
plt.scatter(x[-1]/AU, y[-1]/AU, zorder = 9, s = 60, color = 'blue', label = 'Final')
plt.xlabel('x (AU)', fontsize = 18)
plt.ylabel('y (AU)', fontsize = 18)
plt.grid()
plt.legend(loc = 'upper right')
plt.show()

