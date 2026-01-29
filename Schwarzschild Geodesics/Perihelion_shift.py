import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import G, c

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
r0 = 50 #far enough away that the analytical phi prediction should be close to simulation value
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
    dx_values = []
    for n_step in n_values:
        phi_test, u_test, du_test, dx_test = RK4_2nd(dudphi, wrapped_dwdphi, initial_conditions, x_final=2*orbits*np.pi, n=n_step)
        dx_values.append(dx_test)

        indices, _ = find_peaks(u_test)
        indices = np.insert(indices, 0, 0)
        
        current_peaks = [] 
        for i in indices:
            # Only interpolate if there are neighbours
            if 0 < i < len(u_test) - 1:
                y1, y2, y3 = u_test[i - 1], u_test[i], u_test[i + 1]
                denom = 2*(y1 - 2*y2 + y3)
                if denom != 0:
                    offset = dx_test*(y1 - y3)/denom
                    current_peaks.append(phi_test[i] + offset)
                else:
                    current_peaks.append(phi_test[i])
            else:
                current_peaks.append(phi_test[i])
        precessions.append(np.mean(np.diff(current_peaks) - 2*np.pi))  
    return np.array(precessions), np.array(dx_values)

ns = np.array([5000, 10000, 20000, 40000, 80000, 100000])
test, dxs = convergence_test(ns)
errors = np.abs(test[:-1] - test[-1])
dxs = dxs[:-1]

print()
print("Convergence Tests")
slope = np.polyfit(np.log(dxs), np.log(errors), 1)[0]
print(f"Precession measurement convergence rate (slope): {slope:.2f}")
print("Convergence test uses quadratic interpolation which is second order so this affects the " \
"slope and brings it towards 2, away from 4.")

def energy_drift_test(n_values):
    drifts = []
    dx_values = []
    for n_step in n_values:
        phi_t, u_t, du_t, dx_t = RK4_2nd(dudphi, wrapped_dwdphi, initial_conditions, x_final=2*orbits*np.pi, n=n_step)
        
        # Calculate Energy at the very first and last point
        E_start = specific_energy(u_t[0], du_t[0], M, h)
        E_end = specific_energy(u_t[-1], du_t[-1], M, h)
        
        # The 'Error' is how much Energy was "lost" or "gained" due to numerical drift
        drifts.append(np.abs(E_end - E_start))
        dx_values.append(dx_t)
        
    return np.array(drifts), np.array(dx_values)

# Use these for the plot
n = [500, 1000, 1500, 2000, 2500]
e_errors, e_dxs = energy_drift_test(n)
e_slope = np.polyfit(np.log(e_dxs), np.log(e_errors), 1)[0]
print(f"Energy measurement convergence rate (slope): {e_slope:.2f}")

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

plt.figure(figsize=(8, 5))
plt.loglog(dxs, errors, 'o-', label='Numerical Error')
plt.loglog(dxs, (errors[0]/dxs[0]**4)*dxs**4, 'k--', alpha=0.5, label=r'Theoretical RK4 Slope $(\mathcal{O}(dx^4))$')
plt.xlabel(r'Step Size ($\Delta \phi$)', fontsize=14)
plt.ylabel(r'Precession Error $(rad)$', fontsize=14)
plt.title('RK4 Convergence Test', fontsize=16)
plt.legend()
plt.grid()

plt.figure(figsize=(8, 5))
plt.loglog(e_dxs, e_errors, 's-', label='Energy Drift')
plt.loglog(e_dxs, (e_errors[0]/e_dxs[0]**4) * e_dxs**4, 'k--', alpha=0.5, label='Theoretical $O(dx^4)$')
plt.title('Energy Conservation Convergence', fontsize=16)
plt.xlabel(r'Step Size ($\Delta \phi$)', fontsize=14)
plt.ylabel('Absolute Energy Drift', fontsize=14)
plt.legend()
plt.grid(True, which="both", alpha=0.3)
plt.show()

# plt.figure()
# plt.scatter(phi, np.log10(E), color = 'orange', label = r'$E(\phi)$')
# plt.xlabel(r'$\phi$', fontsize = 20)
# plt.ylabel(r'$\log_{10}(E)$', fontsize = 20)
# plt.title('Specific Energy, log scale', fontsize = 20)
# plt.legend()
# plt.tight_layout()
# plt.grid()
# plt.show()

