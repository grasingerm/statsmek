import random, math, pylab

alpha = 0.5
delta = 0.10
nsamples = 10000000
samples_xy = [[0]*2 for i in range(nsamples)]
samples_xy[0][:] = [0.0, 0.0] # initial conditions

x, y = samples_xy[0][:]
old_weight = math.exp(-0.5 * (x**2 + y**2) - alpha * (x**4 + y**4))

for k in range(nsamples-1):
    c = random.choice([0, 1])
    x_trial = samples_xy[k][:]
    x_trial[c] += random.uniform(-delta, delta)
    x, y = x_trial[:]
    new_weight = math.exp(-0.5 * (x **2 + y**2) - alpha * (x**4 + y**4))
    if random.uniform(0.0, 1.0) < new_weight / old_weight:
        samples_xy[k+1][:] = x_trial[:]
    else:
        samples_xy[k+1][:] = samples_xy[k][:]

samples_x = [xy[0] for xy in samples_xy]
samples_y = [xy[1] for xy in samples_xy]
pylab.hexbin(samples_x, samples_y, gridsize=50, bins=1000)
pylab.axis([-1.0, 1.0, -1.0, 1.0])
cb = pylab.colorbar()
pylab.xlabel('x')
pylab.ylabel('y')
pylab.title('A2_0')
pylab.savefig('plot_A2_0.png')
pylab.show()
