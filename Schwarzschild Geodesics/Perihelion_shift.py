import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import G, c
from scipy.signal import find_peaks

#Mercury-Sun system isn't relativistic enough to accurately calculate the precession
#angle using RK4 so I'll make up a test system to show that the simulation works
def dudphi(phi, u, w):
    return w 

def dwdphi(phi, u, w, M: float, h: float, b):
    """b = 1 for matter, b = 0 for light"""
    return -u + b*M/h**2 + 3*M*u**2

#RK4 solver for 2nd order differential equations
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
    return x, y, dy, dx

#Simulation parameters
M = 1.0
r0 = 50 #far enough away that the analytical phi prediction should match simulation value
u0 = 1/r0 #but close enough that precession is resolvable by RK4 method
h = 8 #stable orbit 
orbits = 10
b = 1
#Analytical precession
dphi_rad_analytical = 2*np.pi/(1 - 3*M**2/h**2) - 2*np.pi

wrapped_dwdphi = lambda phi, u, w: dwdphi(phi, u, w, M, h, b)
initial_conditions = [0.0, u0, 0.0]
phi, u, du, dx = RK4_2nd(dudphi, wrapped_dwdphi, initial_conditions, x_final = 2*orbits*np.pi, n = 100000)
phi = np.array(phi)
u = np.array(u)
du = np.array(du)

r = 1/u
x = r*np.cos(phi)
y = r*np.sin(phi)

# Maxima in u are minima in r ie perihelia
peaks, _ = find_peaks(u)
peaks = np.insert(peaks, 0, 0)
if len(peaks) >= 2:
    phi_peaks = phi[peaks]
    dphi = np.diff(phi_peaks) - 2*np.pi
    mean_precession = np.mean(dphi)
    precession_std = np.std(dphi)
    print('Precession')
    print(f"Numerical precession: {mean_precession:.3f} rad +- {precession_std:.5f} rad")
    print(f"Analytical precession: {dphi_rad_analytical:.3f} rad")
    print(f"Difference: {np.abs(mean_precession - dphi_rad_analytical)*100/dphi_rad_analytical:.2f} %")
else:
    raise ValueError("Not enough perihelia found. Increase 'orbits' or make sure motion is bounded.")

#From the Schwarzschild metric and Lagrangian, the specific energy, E, should be conserved
def specific_energy(u, du, M, h):
    kinetic = h**2*(du)**2
    potential = (1 + (h*u)**2)*(1 - 2*M*u)
    return np.sqrt(kinetic + potential)

E = specific_energy(u, du, M, h)
print()
print('Specific Energy')
print(f'Initial specific energy: {E[0]:.3f}')
print(f'Final specific energy: {E[-1]:.3f}')
print(f'Maximum variation: {np.max(np.abs(E - E[0])):.2e}')

#Convergence test
def convergence_test(n_values):
    precessions = []
    for n_step in n_values:
        phi_test, u_test, du_test, dx_test = RK4_2nd(dudphi, wrapped_dwdphi, initial_conditions, x_final=2*orbits*np.pi, n = n_step)
        indices, _ = find_peaks(u_test)
        indices = np.insert(indices, 0, 0)
        perihelia_shift = np.diff(np.array(phi_test)[indices]) - 2*np.pi
        precessions.append(np.mean(perihelia_shift))
    return np.array([precessions])

ns = np.array([5000, 10000, 20000, 40000, 80000, 100000])
test = convergence_test(ns)


plt.figure(figsize = [6,6])
plt.plot(x, y, color = 'green', label = 'Trajectory')
plt.scatter(x[0], y[0], color = 'red', label = 'Initial position', zorder = 10)
plt.scatter(x[-1], y[-1], color = 'blue', label = 'Final position', zorder = 11)
for i in peaks:
    plt.scatter(x[i], y[i], color = 'orange', zorder = 9)
plt.scatter(x[peaks[-1]], y[peaks[-1]], color = 'orange', label = 'Local minima')
plt.plot([0, 50*np.cos(mean_precession)], [0, 50*np.sin(mean_precession)], 'r--', label = 'Test line')
plt.plot(0, 0, 'ko', label = 'Central Mass at (0,0)')
plt.title('Matter particle trajectory'
        '\n'
        'in Schwarzschild Metric', fontsize = 20)
plt.xlabel('x', fontsize = 20)
plt.ylabel('y', fontsize = 20)
plt.axis("equal")
plt.legend(loc = 'lower left')
plt.tight_layout()
plt.grid()

plt.figure()
plt.plot(phi, E, '.', color = 'orange', label = r'$E(\phi)$')
plt.xlabel(r'$\phi$', fontsize = 20)
plt.ylabel(r'$E$', fontsize = 20)
plt.title('Specific Energy', fontsize = 20)
plt.legend()
plt.tight_layout()
plt.grid()

plt.figure()
plt.scatter(phi, np.log10(E), color = 'orange', label = r'$E(\phi)$')
plt.xlabel(r'$\phi$', fontsize = 20)
plt.ylabel(r'$\log_{10}(E)$', fontsize = 20)
plt.title('Specific Energy, log scale', fontsize = 20)
plt.legend()
plt.tight_layout()
plt.grid()
plt.show()

