import numpy as np
from LorentzBoost import lorentz_factor, lorentz_boost, lorentz_boost_arb

# ms = [0.2, 0.3]
# print(len(ms))

# def four_momentum(m, v):
#     #c = 1
#     v_array = np.array(v) if not isinstance(v, (int, float)) else np.array([v])
#     gamma = lorentz_factor(np.linalg.norm(v_array))
#     E = gamma*m
#     if isinstance(v, (int, float)):
#         p_vec = np.array([gamma*m*v, 0.0, 0.0]) #motion in x
#     else:
#         p_vec = gamma*m*v_array[:3]        
#     P = np.concatenate(([E], p_vec))
#     return P

#Two particle momenta
# vel = [0.5, 0.6, 0.01]
# P1 = four_momentum(0.1, 0.5)
# P2 = four_momentum(0.2, 0.6)
# print(P1, P2)

#Elastic Collision- Energy and momentum conserved
#Total energy is conserved: E1 + E2 = E3 + E4
#Total momentum is conserved: p1 + p2 = p3 + p4
#Restrict motion to x direction only so input four momenta need to have scalar v term
def elastic(masses, speeds):
    if len(masses) == 2:
        four_momenta = [] 
        for n in range(len(masses)):
            m = masses[n]
            v = speeds[n]
            v_array = np.array(v) if not isinstance(v, (int, float)) else np.array([v])

            gamma = lorentz_factor(np.linalg.norm(v_array))
            E = gamma*m

            if isinstance(v, (int, float)):
                p_vec = np.array([gamma*m*v, 0.0, 0.0]) #motion in x
            else:
                p_vec = gamma*m*v_array[:3]        
            P = np.concatenate(([E], p_vec))
            four_momenta.append(P)
    else: raise ValueError("Function only works for two particle collisions.")
    
    P1, P2 = four_momenta[0:2]
    total_P = P1 + P2
    v_CoM = []
    v_CoM[0:3] = total_P[1:3]/total_P[0]
    print(v_CoM)
    #CoM components
    total_P_CoM = lorentz_boost_arb(v_CoM, total_P)
    P1_CoM = lorentz_boost_arb(v_CoM, P1)
    P2_CoM = lorentz_boost_arb(v_CoM, P2)
    #After collision
    final_P1_CoM = np.diag([1, -1, -1, -1]) @ P1_CoM
    final_P2_CoM = np.diag([1, -1, -1, -1]) @ P2_CoM
    #Transform back to Lab frame
    final_P1 = lorentz_boost_arb(-v_CoM, final_P1_CoM)
    final_P2 = lorentz_boost_arb(-v_CoM, final_P2_CoM)
    final_v1 = final_P1[1]/final_P1[0]
    final_v2 = final_P2[1]/final_P2[0]

    print(f"The total four-momentum in the CoM frame is {total_P_CoM}.")

    print(f"The final velocity and four-momentum of particle 1 are {final_v1} and {final_P1} in the lab frame.")
    print(f"The final velocity and four-momentum of particle 2 are {final_v2} and {final_P2} in the lab frame.")

    return final_v1, final_v2, final_P1, final_P2

ms = [0.2, 0.2]
vs = [[0.1, 0.1, 0.1], [0.1, 0.1, 0.1]]

elastic(ms, vs)

#Try now with arbitrary 3-velocities

# def elastic_3d(masses, velocities):
#     if len(masses) == 2:
#         four_momenta = [] 
#         for n in range(len(masses)):
#             m = masses[n]
#             v = speeds[n]




#Perfectly Inelastic Collision- KE not conserved eg. two particles merge into one

# def inelastic(P1, P2):
