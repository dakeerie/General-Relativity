import numpy as np
from LorentzBoost import lorentz_factor, lorentz_boost

def four_momentum(m, v):
    #c = 1
    gamma = lorentz_factor(np.linalg.norm(v))
    E = gamma*m
    if type(v) == float:
        p = gamma*m*v
        P = np.array([E, p])
        return P
    else:
        velocity = np.array(v)
        print(gamma*m)
        px, py, pz = gamma*m*velocity[0:3]
        P = np.array([E, px, py, pz])
        return P

#Two particle momenta

vel = [0.1, 0.2, 0.3]
# P1 = four_momentum(0.1, 0.5)
P2 = four_momentum(0.2, vel)
print(P1, P2)