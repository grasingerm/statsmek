#!/usr/bin/env python3

import random
import numpy as np

from time import clock

t0 = clock()

trials = 60000000
sides  = 6

print ("\nNumber of trials = {:d}\n".format(trials))

histogram = np.zeros(sides,int)   #initialize the histogram array

for j in range(trials):
    r = int(random.random() * sides )   # this goes from 0 to sides-1
    histogram[r] +=  1

print("eyes    counts     abs. dev.     rel. dev.")
for j in range(sides):
    print ("{:3d}    {:7d} {:14.6f}   {: 12.7f}"\
           .format(j+1,   # add 1 to go back from 1 to sides
                   histogram[j],
                   histogram[j] - trials/sides,
                  (histogram[j] - trials/sides) / trials))

t_end = clock()
print ("\nSimulation time = {:1.3f} seconds".format(t_end - t0))
