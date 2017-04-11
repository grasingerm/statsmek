import random, math, pylab

def gauss_cut():
    while True:
        x = random.gauss(0.0, 1.0)
        if abs(x) <= 1.0:
            return x

alpha = 0.5
nsamples = 1000000
samples_x = [0]*nsamples
samples_y = [0]*nsamples
for k in range(nsamples):
    while True:
        x = gauss_cut()
        y = gauss_cut()
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
pylab.title('A1_3')
pylab.savefig('plot_A1_3.png')
pylab.show()
