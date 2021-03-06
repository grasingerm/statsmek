import math, random, pylab

def rho_free(x, y, beta):
    return math.exp(-(x - y)**2 / (2.0 * beta))

def levy_harmonic_path(xstart, xend, dtau, N):
    x = [xstart]
    for k in range(1, N):
        dtau_prime = (N - k) * dtau
        Ups1 = 1.0 / math.tanh(dtau) + 1.0 / math.tanh(dtau_prime)
        Ups2 = x[k-1] / math.sinh(dtau) + xend / math.sinh(dtau_prime)
        x.append(random.gauss(Ups2 / Ups1, 1.0 / math.sqrt(Ups1)))
    return x

beta = 20.0
sigma = 1.0 / math.sqrt(2.0 * math.tanh(beta / 2.0))
N = 5
dtau = beta / N
delta = 1.0
n_steps = 1000000
x = [5.0] * N
data = []

for step in range(n_steps):
    x[0] = random.gauss(0.0, sigma)
    x = levy_harmonic_path(x[0], x[0], dtau, N)
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
pylab.savefig('plot_B3_beta%s.png' % beta)
pylab.show()
pylab.clf()

pylab.plot(final_path, [dtau * n for n in range(N)])
pylab.xlabel('$x$')
pylab.ylabel('$\\tau$')
pylab.title('levy_harmonic_path ($\\beta=%s, N=%i$)' % (beta, N))
pylab.savefig('plot_B3_beta%s_final-path.png' % beta)
pylab.show()
