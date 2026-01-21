import scipy.constants as const
import numpy as np
c = const.speed_of_light

def lorentz_factor(v):
    #v in natural units
    if np.linalg.norm(v) <= 1:
        gamma = 1/np.sqrt(1-v**2)
    else:
        raise ValueError(f"velocity v={v} must satisfy |v| <= 1.")
    return gamma

def lorentz_boost(v, event):
    #v in natural units
    #Lorentz boost of velocity v in the x direction
    t, x, y, z = event[0:4]
    gamma = lorentz_factor(v)

    #Matrix Equation
    lorentz_matrix = np.array([[gamma, -gamma*v, 0, 0], [-gamma*v, gamma, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    barred_event = np.matmul(lorentz_matrix, event.transpose())
    #Explicit equations to check correct matrix multiplication
    # barred_coordinates = np.zeros(4)
    # barred_coordinates[0] = gamma*(t - v*x)
    # barred_coordinates[1] = gamma*(-v*t + x)
    # barred_coordinates[2] = y #y and z coordinates the same for boost in x direction alone
    # barred_coordinates[3] = z
    # t_prime, x_prime, y_prime, z_prime = barred_coordinates[0:4]
    
    return barred_event

if __name__ == "__main__":
    test = np.array([1, 1, 1, 1])
    barred = lorentz_boost(0.5, test)
    print(barred)
else:
    print("lorentz_boost, lorentz_factor  and lorentz_boost_arb were imported.")



#Arbitrary direction velocity boost
def lorentz_factor_arb(v_vec):
    mag = np.linalg.norm(v_vec)
    return 1/np.sqrt(1- mag**2)

def lorentz_boost_arb(v_vec, event):
    #Rapidity
    beta = np.array(v_vec)
    beta_mag = np.linalg.norm(beta)
    if beta_mag >= 1:
        raise ValueError("Speed must be less than the speed of light ie |v|<1 when expressed in natural units")

    gamma = lorentz_factor_arb(beta)
    t, x, y, z = event
    position_vec = np.array([x, y, z])

    boost_matrix = np.eye(4)
    boost_matrix[0, 0] = gamma
    boost_matrix[0, 1:4] = (-gamma * beta).reshape(3,)
    boost_matrix[1:4, 0] = (-gamma * beta).reshape(3,)

    outer = np.outer(beta, beta)
    if beta_mag > 0:
        boost_matrix[1:4, 1:4] += ((gamma-1)/beta_mag**2)*outer
    # print(boost_matrix)

    return boost_matrix @ event

v = [0.1, 0.1, 0.1]
e = [1, 1, 1, 1]
if __name__ == "__main__":
    print(lorentz_boost_arb(v, e))





#Add boost in any direction. Velocity 3-vector?
#Look into rapidity 
