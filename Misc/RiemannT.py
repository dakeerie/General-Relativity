from Week1.test import function
from Misc.ChristoffelSymbols import *
# IMPORTING FUNCTIONS NEEDS WORK!!

import sympy as sp
from sympy import *
init_printing(use_unicode=True)
# from sympy import MutableDenseNDimArray

#Redefine Christoffel Function
def christoffel_symbols(g, coords):
    n = len(coords)
    g_inv = g.inv()
    Gamma = MutableDenseNDimArray.zeros(n, n, n)
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for l in range(n):
                    term1 = sp.diff(g[l,j], coords[k])
                    term2 = sp.diff(g[l,k], coords[j])
                    term3 = sp.diff(g[j,k], coords[l])
                    Gamma[i,j,k] += 0.5*g_inv[i, l]*(term1+term2-term3)
                    Gamma[i,j,k] = sp.simplify(sp.trigsimp(Gamma[i,j,k]))
    return Gamma

#printer function
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

#Riemann tensor
def Riemann_Curvature(g, coords, Gamma):
    n = len(coords)
    R = Gamma = MutableDenseNDimArray.zeros(n, n, n, n)
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for l in range(n):
                    for m in range(n):
                        term1 = sp.diff(Gamma[i, j, l], coords[k])
                        term2 = sp.diff(Gamma[i, j, k], coords[l])
                        term3 = Gamma[i, m, k]*Gamma[m, j, l]
                        term4 = Gamma[i, m, l]*Gamma[m, j, k]
                        R[i, j, k, l] += term1 - term2 + term3 - term4
    return R

#Ricci tensor
# def Ricci_Tensor(g, R):


