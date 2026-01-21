import numpy as np 
from einsteinpy.geodesic import Nulllike
from einsteinpy.plotting.geodesic import GeodesicPlotter, StaticGeodesicPlotter, InteractiveGeodesicPlotter

X = [10, np.pi/2 , 0] #position vector, set initial position radial value as 5.86, unstable circular orbit
# at r = 4, stable circular orbit at r = 12
P = [0., 0., 4] #momentum set p_phi = 4

a = 0.
steps = 100 #good plot set steps = 1000, 26 steps for unstable circular orbit and 227 for stable
delta = 1.

geod = Nulllike(metric = 'Schwarzschild', metric_params=(a,), position = X, momentum = P, steps = steps,
                delta = delta, return_cartesian = True)

gpl = GeodesicPlotter(ax=None, bh_colors=('#000', '#FFC'), draw_ergosphere=False)
#%%

# plot
gpl.plot2D(geod, coordinates=(1, 2), color= 'blue')
gpl.show()
# gpl.plot(geod, color= 'blue')
# gpl.clear()
#%%
#parameters plot
gpl.parametric_plot(geod) #this will plot the variables values
gpl.show()

#%%
import matplotlib.pyplot as plt
trajectory = geod.trajectory[1]
x1_list = []
x2_list = []
iterations = []
for i in range(0,steps):
    x1 = trajectory[i][1] # X1 values
    x2 = trajectory[i][2] # X2 values 
    ite = i # keep the iteartions
    x1_list.append(x1)
    x2_list.append(x2)
    iterations.append(ite)
# plotting the results
plt.plot(iterations, x1_list, color = 'red', label = r'$X_1$ (cartesian)')
plt.plot(iterations, x2_list, color = 'blue', label = r'$X_2$ (cartesian)')
plt.legend(loc = 'upper right', bbox_to_anchor = (1.5, 1))
plt.title(r'$X_1$ and $X_2$ in cartesian')
plt.xlabel(r'Affine parameter $\lambda$')
plt.ylabel('Coordinates')
plt.show()
#%%
rs = 2 # since G=M=c=1 
x2_negative = []
for i in x2_list:
    d = -i
    x2_negative.append(d)

circle = plt.Circle((0, 0), rs, color='black', alpha=0.5, label = r'Singularity of radius $R_s$')
# Plot the trajectory in the xy-plane
plt.plot(x1_list, x2_list, label = 'orbit', color = 'orange')
plt.plot(x1_list, x2_negative, color = 'orange')
plt.axhline(0, color='black',linewidth=0.5)
plt.axvline(0, color='black',linewidth=0.5)
plt.gca().set_aspect('equal', adjustable='box')
plt.gcf().gca().add_artist(circle)

lim = 6 # set 20 or 5 for unstable circular orbit
plt.xlim(-lim,lim)
plt.ylim(-lim,lim)

#plt.scatter(x1_list[-1], x2_list[-1], label = "Particle's current position", color = 'black', s = 15 )
#plt.scatter(x1_list[0], x2_list[0], label = "Photon", color = 'blue', s = 15)

plt.ylabel(r'$\frac{y}{R_s}$')
plt.xlabel(r'$\frac{x}{R_s}$')
plt.title('Gravitational deflection of light')
#plt.title('Unstable circular orbit')
#plt.title('Stable circular orbit')
plt.legend(loc='upper right', bbox_to_anchor=(1.8, 1))
plt.show()

#%%
# =============================================================================
# forcasting values
# =============================================================================

# is important to notice that the minimal radius is 5.831, once this is reached 
# from solving the geodesic equation, the iterations will stop leaving a gap in the plot 
# if we try to provide the full plot using the symmetries of the Schwarzschilds black hole. Hence
# we require to forecast the values for the missing gap produced from the above code
print(x2_list[0], x2_list[1], x2_list[2], x2_list[3], x2_list[4])
# note from the above it is quite clear to see that the values of the x2 coordinate (y) is 
# doubling, therefore, if x1 values reached its minimum, we can forecast the missing values for x2. 

value1 = []
value2 = []
xval = x1_list[0]
yval = x2_list[0]

for i in range(0,25):# generate iterative loop
    yval = yval/2
    value1.append(xval)#append constant values to a list
    value2.append(yval)# append forecasted values in a list
    
print(value1)
print(value2)
print(x1_list[0])
#%%
xs = value1 + x1_list 
ys = value2 + x2_list
print(xs)
rs = 2 # since G=M=c=1 
ys_negative = []
for i in ys:
    d = -i
    ys_negative.append(d)

circle = plt.Circle((0, 0), rs, color='black', alpha=0.5, label = r'Singularity of radius $R_s$')
# Plot the trajectory in the xy-plane
plt.plot(xs, ys, label = 'light deflection', color = 'orange')
plt.plot(xs, ys_negative, color = 'orange')
plt.axhline(0, color='black',linewidth=0.5)
plt.axvline(0, color='black',linewidth=0.5)
plt.gca().set_aspect('equal', adjustable='box')
plt.gcf().gca().add_artist(circle)

lim = 20 # set 20 or 5 for unstable circular orbit
plt.xlim(-lim,lim)
plt.ylim(-lim,lim)

#plt.scatter(x1_list[-1], x2_list[-1], label = "Particle's current position", color = 'black', s = 15 )
#plt.scatter(x1_list[0], x2_list[0], label = "Photon", color = 'blue', s = 15)

plt.ylabel(r'$\frac{y}{R_s}$')
plt.xlabel(r'$\frac{x}{R_s}$')
plt.title('Gravitational deflection of light')
#plt.title('Unstable circular orbit')
#plt.title('Stable circular orbit')
plt.legend(loc='upper right', bbox_to_anchor=(1.8, 1))
plt.show()
