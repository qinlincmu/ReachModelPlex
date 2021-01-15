#!/usr/bin/env python
#author: qinlin@andrew.cmu.edu
import numpy as np
import math

def positive(x):
    return max(x, 0)
def range_by_step(start, end, step):
    range = []
    if end >= start:
        range = [item for item in np.arange(start, end, step)]
    else:
        range = [item for item in np.arange(start, end, -1*step)]
    return range
def pi_2_pi(angle):
    #to[0, pi] or [-pi, 0]
    return (angle + math.pi) % (2 * math.pi) - math.pi
