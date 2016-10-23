from __future__ import print_function
import random
import sys

def direct_pi(nsamples):
    nhits = 0
    for i in range(nsamples):
        x = random.uniform(-1, 1)
        y = random.uniform(-1, 1)
        if x**2 + y**2 < 1:
            nhits = nhits + 1
    return 4.0 * nhits / float(nsamples)

assert(len(sys.argv) == 2)
try:
    print(direct_pi(int(sys.argv[1])))
except:
    print("invalid argument ", sys.argv[1], file=sys.stderr)
    raise
