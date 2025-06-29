import numpy as np
from LorentzBoost import lorentz_factor, lorentz_boost

def four_momentum(m, v):
    #c = 1
    v_array = np.array(v) if not isinstance(v, (int, float)) else np.array([v])

    gamma = lorentz_factor(np.linalg.norm(v_array))
    E = gamma*m

    if isinstance(v, (int, float)):
        p_vec = np.array([gamma*m*v, 0.0, 0.0]) #motion in x
    else:
        p_vec = gamma*m*v_array[:3]        
    P = np.concatenate(([E], p_vec))
        
    return P


#Two particle momenta

vel = [0, 0, 0]
P1 = four_momentum(0.1, 0.5)
P2 = four_momentum(0.2, vel)
print(P1, P2)