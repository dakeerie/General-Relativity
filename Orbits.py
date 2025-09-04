import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from einsteinpy.geodesic import Timelike
from scipy.constants import G, c
from scipy.integrate import solve_ivp


def system(s, y, M, h, object):

    if object == 'matter':
        b = 1
    if object == 'photon':
        b = 0
    else:
        raise ValueError("Object specification must be either 'matter' or 'photon'.")
    
    u, dudphi = y
    q = dudphi
    dq = -u + b*M/h**2 +3*M*u**2

    return [q, dq]





    