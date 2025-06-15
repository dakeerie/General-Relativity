import scipy.constants as const
import numpy as np
c = const.speed_of_light

def lorentz_factor(v):
    #v in natural units
    gamma = (1-v**2)**(-1/2)
    return gamma


array = np.array([1,2,3,4])
a,b,c,d = array[0:4]
print(c)

def lorentz_boost(v, event):
    #v in natural units
    #lorentz boost of velocity v in the x direction
    # t, x, y, z = event[0:3]
    barred_coordinates = np.zeros(4)
    gamma = lorentz_factor(v)

    # t,x,y,z = event[0:3]

    barred_coordinates[0] = gamma*event[0] - gamma*v*event[1]
    barred_coordinates[1] = -gamma*v*event[0] + gamma*event[1]
    barred_coordinates[2] = event[2] #y and z coordinates the same for boost in x direction alone
    barred_coordinates[3] = event[3]
    #CHECK FORMULAE AND T, X, Y, Z LINE!!
    return barred_coordinates

#Add boost in any direction. Velocity 3-vector?
#Look into rapidity 
