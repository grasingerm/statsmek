import math, random, pylab

def rho_free(x, y, beta):
    return math.exp(-(x - y)**2 / (2.0 * beta))

def levy_free_path(xstart, xend, dtau, N):
    x = [xstart]
    for k in range(1, N):
        dtau_prime = (N - k) * dtau
        x_mean = (dtau_prime * x[k - 1] + dtau * end) / (dtau + dtau_prime)
        sigma = math.sqrt(1.0 / (1.0 / dtau + 1.0  dtau_prime))
        x.append(random.gauss(Ups2 / Ups1, 1.0 / math.sqrt(Ups1)))
    return x

beta = 20.0
N = 2
Ncut = N / 2
dtau = beta / N
delta = 1.0
n_steps = 1000000
x = [5.0] * N
data = []

for step in range(n_steps):
    x = levy_harmonic_path(x[0], x[0], dtau, N)
    x = x[Ncut:] + x[:Ncut]
    if step % N == 0:
        k = random.randint(0, N-1)
        data.append(x[k])
    if step % 100000 == 0:
        print 'step ', step

final_path = x[:]
pylab.hist(data, normed=True, bins=100, label='QMC')
list_x = [0.1 * a for a in range(-30, 31)]
list_y = [math.sqrt(math.tanh(beta / 2.0)) / math.sqrt(math.pi) * \
          math.exp(-x ** 2 * math.tanh(beta / 2.0)) for x in list_x]
pylab.plot(list_x, list_y, label='analytic')
pylab.legend()
pylab.xlabel('$x$')
pylab.ylabel('$\\pi(x)$ (normalized)')
pylab.title('levy_harmonic_path ($\\beta=%s, N=%i$)' % (beta, N))
pylab.xlim(-2, 2)
pylab.savefig('plot_B2_beta%s.png' % beta)
pylab.show()
pylab.clf()

pylab.plot(final_path, [dtau * n for n in range(N)])
pylab.xlabel('$x$')
pylab.ylabel('$\\tau$')
pylab.title('levy_harmonic_path ($\\beta=%s, N=%i$)' % (beta, N))
pylab.savefig('plot_B2_beta%s_final-path.png' % beta)
pylab.show()
