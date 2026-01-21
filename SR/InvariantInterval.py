import numpy as np
from LorentzBoost import lorentz_boost, lorentz_factor

#Minkowski metric
eta = np.diag([-1, 1, 1, 1])

def interval(event1, event2):
    coordinate_diff  = event2 - event1
    dual = eta @ coordinate_diff
    s_sq = dual @ coordinate_diff
    return s_sq

if __name__ == "__main__":
    test1 = np.array([1, 1, 1, 1])
    test2 = np.array([0, 0, 0, 0])
    variable = interval(test1, test2)
    print(variable)
else:
    print("interval was imported.")

origin = np.zeros(4)
unprimed = np.array([0, 1, 0, 0])
primed = lorentz_boost(0.5, unprimed)

#Commented out specific example
# interval_unprime = interval(origin, unprimed)
# interval_prime = interval(origin, primed)
# print(interval_prime, interval_unprime)

# if np.isclose(interval_unprime, interval_prime):
#     #used isclose instead of == as if statement wouldn't produce the correct print out due to 
#     #floating point error.
#     print('The interval is invariant under the Lorentz transformation.')
# else:
#     print("Something's gone wrong.")

#Put into function

def check(event1, event2, v):
    primed1 = lorentz_boost(v, event1)
    primed2 = lorentz_boost(v, event2)

    interval_unprime = interval(event1, event2)
    interval_prime = interval(primed1, primed2)

    if np.isclose(interval_prime, interval_unprime):
        print('The interval is invariant under the Lorentz transformation.')
    else:
        print("Something's gone wrong!")

if __name__ == "__main__":
    check(origin, unprimed, 0.5)

#go on to start with a set of events in one frame, calculate the interval,
#then transform to barred basis and calculate again to verify the interval is conserved between frames.


