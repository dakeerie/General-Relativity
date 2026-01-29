import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import G, c

#Light deflection by massive object
def dudphi(phi, u, w):
    return w 

def dwdphi(phi, u, w, M: float):
    return -u + 3*M*u**2

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
b_imp = 20.0
phi0 = -np.pi/2
u0 = np.sin(phi0)/b_imp
w0 = np.cos(phi0)/b_imp
r0 = 1/u0
analytical = 4*M/b_imp
print("Analytical deflection:", analytical)

wrapped_dwdphi = lambda phi, u, w: dwdphi(phi, u, w, M)
initial_conditions = [phi0, u0, w0]
phi, u, du, dx = RK4_2nd(dudphi, wrapped_dwdphi, initial_conditions, x_final = 1.5*np.pi, n = 100000)
phi = np.array(phi)
u = np.array(u)
du = np.array(du)
r = 1/u

phif = phi[-1]
rf = r[-1]
wf = np.cos(phi[-1])/b_imp
drdphif = -rf**2*wf
num = drdphif*np.sin(phif) + rf*np.cos(phif)
den = drdphif*np.cos(phif) - rf*np.sin(phif)
cart_grad = num/den
xf = rf*np.cos(phif)
yf = rf*np.sin(phif)
x_coords = np.linspace(xf, xf + 2000, 2000)
y_coords = yf + cart_grad*(x_coords - xf)
print(xf)

u_cutoff = 1e-3
mask = (u > u_cutoff)
x = r[mask]*np.cos(phi[mask])
y = r[mask]*np.sin(phi[mask])

# mask_far = r > 10  # adjust depending on scale
# x_far = x[mask_far]
# y_far = y[mask_far]
# coeffs = np.polyfit(x_far[-100:], y_far[-100:], 1)
# grad_out = coeffs[0]

# theta_in = np.arctan(cart_grad)
# theta_out = np.arctan(grad_out)
# delta_phi = theta_in - theta_out
# print("Numerical deflection:", delta_phi)

plt.figure(figsize = [6,6])
plt.plot(x, y, color = 'green', label = 'Trajectory')
plt.scatter(x[0], y[0], color = 'red', label = 'Initial position', zorder = 10)
plt.scatter(x[-1], y[-1], color = 'blue', label = 'Final position', zorder = 11)
plt.plot(0, 0, 'ko', ms = 2.5, label = 'Central Mass at (0,0)')
plt.plot(x_coords, y_coords, '--', color='grey', label='Flat-space trajectory')
plt.title('Photon trajectory'
        '\n'
        'in Schwarzschild Metric', fontsize = 20)
plt.xlabel('x', fontsize = 20)
plt.ylabel('y', fontsize = 20)
plt.axis("equal")
plt.legend(loc = 'lower left')
plt.tight_layout()
plt.grid()
plt.show()