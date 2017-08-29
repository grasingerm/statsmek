import math, random, pylab

cubic = -1
quartic = 1

def V_anharm(x, cubic, quartic):
    pot = cubic * x ** 3 + quartic * x ** 4
    return pot

def levy_harmonic_path(xstart, xend, dtau, N):
    x = [xstart]
    for k in range(1, N):
        dtau_prime = (N - k) * dtau
        Ups1 = 1.0 / math.tanh(dtau) + 1.0 / math.tanh(dtau_prime)
        Ups2 = x[k-1] / math.sinh(dtau) + xend / math.sinh(dtau_prime)
        x.append(random.gauss(Ups2 / Ups1, 1.0 / math.sqrt(Ups1)))
    return x

def weight(x):
    global cubic
    global quartic
    return math.exp(sum(-V_anharm(a, cubic, quartic) * dtau for a in x))

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
    x_new = levy_harmonic_path(x[0], x[Ncut], dtau, Ncut) + x[Ncut:]
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

final_path = x[:]
pylab.hist(data, normed=True, bins=100, label='QMC')
pylab.legend()
pylab.xlabel('$x$')
pylab.ylabel('$\\pi(x)$ (normalized)')
pylab.title('levy_anharmonic_path ($\\beta=%s, N=%i$)' % (beta, N))
pylab.savefig('c2.png');
pylab.show()
pylab.clf()
