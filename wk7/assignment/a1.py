import math, random, pylab

def pi_x(x, beta):
    sigma = 1.0 / math.sqrt(2.0 * math.tanh(beta / 2.0))
    return math.exp(-x ** 2 / (2.0 * sigma ** 2)) / math.sqrt(2.0 * math.pi) / sigma

def levy_harmonic_path(k):
    x = [random.gauss(0.0, 1.0 / math.sqrt(2.0 * math.tanh(k * beta / 2.0)))]
    if k == 2:
        Ups1 = 2.0 / math.tanh(beta)
        Ups2 = 2.0 * x[0] / math.sinh(beta)
        x.append(random.gauss(Ups2 / Ups1, 1.0 / math.sqrt(Ups1)))
    return x[:]

beta = 2.0
nsteps = 1000000
low = levy_harmonic_path(2)
high = low[:]
data = []
nsamples = 0
for step in xrange(nsteps):
    k = random.choice([0, 1])
    low[k] = levy_harmonic_path(1)[0]
    high[k] = low[k]
    data.append(low[:])
    nsamples += 1
    if step % 100000 == 0:
        print 'step ', step

pylab.hist([data[k][0] for k in range(nsamples)], bins=50, normed=True, label='Particle 1')
pylab.hist([data[k][1] for k in range(nsamples)], bins=50, normed=True, label='Particle 2')
x_min = math.floor(min(min([data[k][0] for k in range(nsamples)]), min([data[k][1] for k in range(nsamples)])))
x_max = math.ceil(max(max([data[k][0] for k in range(nsamples)]), max([data[k][1] for k in range(nsamples)])))
nx = 200
dx = (x_max - x_min) / nx
list_x = [dx * i + x_min for i in range(nx)]
list_y = [pi_x(x, beta) for x in list_x]
pylab.plot(list_x, list_y, label="analytical")
pylab.legend()
pylab.show()
