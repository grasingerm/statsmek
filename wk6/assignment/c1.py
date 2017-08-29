import math, random, pylab

cubic = -1
quartic = 1

def V(x, cubic, quartic):
    pot = x ** 2 / 2.0 + cubic * x ** 3 + quartic * x ** 4
    return pot

def levy_free_path(xstart, xend, dtau, N):
    x = [xstart]
    for k in range(1, N):
        dtau_prime = (N - k) * dtau
        x_mean = (dtau_prime * x[k - 1] + dtau * xend) / (dtau + dtau_prime)
        sigma = math.sqrt(1.0 / (1.0 / dtau + 1.0 / dtau_prime))
        x.append(random.gauss(x_mean, sigma))
    return x

def weight(x):
    global cubic
    global quartic
    return math.exp(sum(-V(a, cubic, quartic) * dtau for a in x))

beta = 20.0
N = 100
Ncut = 15
dtau = beta / N
delta = 1.0
n_steps = 10000000
x = [0.0] * N
old_weight = weight(x)
data = []

for step in range(n_steps):
    x_new = levy_free_path(x[0], x[Ncut], dtau, Ncut) + x[Ncut:]
    new_weight = weight(x_new)
    if random.uniform(0.0, 1.0) < new_weight / old_weight:
        x = x_new[:]
        x = x[1:] + x[:1]
        old_weight = new_weight
    if step % N == 0:
        k = random.randint(0, N-1)
        data.append(x[k])
    if step % 100000 == 0:
        print 'step ', step

pylab.hist(data, normed=True, bins=100)
pylab.legend()
pylab.xlabel('$x$')
pylab.ylabel('$\\pi(x)$ (normalized)')
pylab.title('levy_anharmonic_path ($\\beta=%s, N=%i$)' % (beta, N))
pylab.savefig('c1.png');
pylab.show()
pylab.clf()
