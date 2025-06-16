import numpy as np
from LorentzBoost import lorentz_boost, lorentz_factor

#Minkowski metric
eta = np.diag([-1, 1, 1, 1])

def interval(event1, event2):
    coordinate_diff  = event2 - event1
    dual = np.matmul(eta, coordinate_diff)
    s_sq = np.dot(dual, coordinate_diff)
    return s_sq

test1 = np.array([1, 1, 1, 1])
test2 = np.array([0, 0, 0, 0])
print(interval(test1, test2))

#go on to start with a set of events in one frame, calculate the interval,
#then transform to barred basis and calculate again to verify the interval is conserved between frames.


