import random, math, pylab

def prob(x):
    s1 = math.exp(-(x + 1.2) ** 2 / 0.72)
    s2 = math.exp(-(x - 1.5) ** 2 / 0.08)
    return (s1 + 2.0 * s2) / math.sqrt(2.0 * math.pi)

delta = 3.0
nsteps = 1000000
acc_tot = 0
x = 0.0
x_av = 0.0
data = []
for step in xrange(nsteps):
    xnew = x + random.uniform(-delta, delta)
    if random.uniform(0.0, 1.0) < prob(xnew) / prob(x):
        x = xnew
        acc_tot += 1
    x_av += x
    data.append(x)

print 'global acceptance ratio:', acc_tot / float(nsteps)
print '<x> =', x_av / float(nsteps)
pylab.hist(data, bins=100, normed=True)
pylab.grid()
pylab.show()
