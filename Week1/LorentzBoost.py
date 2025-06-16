import scipy.constants as const
import numpy as np
c = const.speed_of_light

def lorentz_factor(v):
    #v in natural units
    gamma = (1-v**2)**(-1/2)
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

test = np.array([1, 1, 1, 1])
barred = lorentz_boost(0.5, test)
print(barred)



#Add boost in any direction. Velocity 3-vector?
#Look into rapidity 
