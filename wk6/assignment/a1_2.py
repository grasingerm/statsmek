import random, math, pylab

alpha = 0.5
nsamples = 1000000
samples_x = [0]*nsamples
samples_y = [0]*nsamples
for k in range(nsamples):
    while True:
        phi = random.uniform(0.0, 2.0*math.pi)
        upsilon = random.uniform(0.0, 1.0)
        psi = -math.log(upsilon)
        r = math.sqrt(2*psi)
        x = r * math.cos(phi)
        y = r * math.sin(phi)
        p = math.exp(-alpha * (x**4 + y**4))
        if random.uniform(0.0, 1.0) < p:
            break
    samples_x[k] = x
    samples_y[k] = y

pylab.hexbin(samples_x, samples_y, gridsize=50, bins=1000)
pylab.axis([-1.0, 1.0, -1.0, 1.0])
cb = pylab.colorbar()
pylab.xlabel('x')
pylab.ylabel('y')
pylab.title('A1_2')
pylab.savefig('plot_A1_2.png')
pylab.show()
