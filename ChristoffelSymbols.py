import sympy as sp
from sympy import *
init_printing(use_unicode=True)
# from sympy import MutableDenseNDimArray

r, theta = sp.symbols('r theta')
coords = [r, theta]

g_polar = sp.Matrix([[1, 0], [0, r**2]])

def christoffel_symbols(g, coords):
    n = len(coords)
    g_inv = g.inv()
    Gamma = MutableDenseNDimArray.zeros(n, n, n)
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for l in range(n):
                    Gamma[i, j, k] += 0.5*g_inv[i, l]*(sp.diff(g[l,j], coords[k]) + sp.diff(g[l,k], coords[l]) - sp.diff(g[j,k], coords[l])) #+= is required as we're summing over different l contributions
                    Gamma[i,j,k] = sp.simplify(sp.trigsimp(Gamma[i,j,k]))
    return Gamma

Gamma = christoffel_symbols(g_polar, coords)
print(Gamma[0, 1, 1])           

def display_christoffel(Gamma, coords):
    n = len(coords)
    print(n)
    for i in range(n):
        for j in range(n):
            for k in range(n):
                if Gamma[i,j,k] != 0:
                    print(f"Γ^{i}_[{j}{k}] = {Gamma[i,j,k]}")
    #tried doing coords[i, j, k] but it looks better with just 0s and 1s in the printout
    print("Rest are zero.")
    print()

display_christoffel(Gamma, coords)

#Now do same for spherical coordinates
r, theta, phi = sp.symbols('r theta phi')
coords_spherical = [r, theta, phi]

g_spherical = sp.Matrix([[1, 0, 0], [0, r**2, 0], [0, 0, r**2 * sp.sin(theta)**2]])

Gamma_spherical = christoffel_symbols(g_spherical, coords_spherical)
display_christoffel(Gamma_spherical, coords_spherical)

#Schwarzschild Metric

t, r, theta, phi = sp.symbols('t r theta phi')
M = sp.symbols('M')
# G=c=1
coords_schwarz = [t, r, theta, phi]
#using [-, +, +, +] signature
g_schwarz = sp.Matrix([[-(1-2*M/r), 0, 0, 0], [0, (1-2*M/r)**(-1), 0, 0], [0, 0, r**2, 0], [0, 0, 0, (r**2)*(sp.sin(theta))**2]])

Gamma_schwarz = christoffel_symbols(g_schwarz, coords_schwarz)
display_christoffel(Gamma_schwarz, coords_schwarz)
