import sympy as sp
from sympy import *
init_printing(use_unicode=True)

# M = sp.Matrix([[4,2], [1, 0]])
# print(M.shape)
# print(M[0,1])
# print(M)

r, theta = sp.symbols('r theta')

g = sp.Matrix([[1, 0], [0, r**2]])

def christoffel_symbols(g):
    n = g.shape[0]
    Gamma = sp.zeros(n, n, n)

n = g.shape[0]
zeroes = sp.zeros(n, n, n)
print(zeroes)