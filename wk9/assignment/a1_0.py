import random, math

def prob(x):
    s1 = math.exp(-(x + 1.2) ** 2 / 0.72)
    s2 = math.exp(-(x - 1.5) ** 2 / 0.08)
    return (s1 + 2.0 * s2) / math.sqrt(2.0 * math.pi)

delta = 10.0
nsteps = 10000000
acc_tot = 0
x = 0.0
xsum = 0.0
for step in xrange(nsteps):
    x = random.uniform(-1e3, 1e3)
    if random.uniform(0.0, 1.0) < prob(x):
        xsum += x
        acc_tot += 1

print '<x> =', xsum / float(acc_tot)
