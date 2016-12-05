import math, random, pylab

def V_sph(dim):
    return math.pi ** (dim / 2.0) / math.gamma(dim / 2.0 + 1.0)

def Q(dim):
    return V_sph(dim) / V_sph(dim-1)

def markov_chain_qd(d, n_trials, delta=0.1):
    n_hits = 0
    xs = [0.0] * (d - 1)
    rsq = 0.0
    for i in range(n_trials):

        # trial move
        k = random.randint(0, d - 2)
        xk = xs[k]
        new_xk = xk + random.uniform(-delta, delta)
        new_rsq = rsq + new_xk*new_xk - xk*xk 
        if new_rsq < 1.0:
            xs[k] = new_xk
            rsq = new_rsq

        alpha = random.uniform(-1.0, 1.0)
        if rsq + alpha*alpha < 1: n_hits += 1

    return 2 * n_hits / float(n_trials)

Qapproxs = []
for d in range(1, 200):
    Qapproxs.append(markov_chain_qd(d+1, 1000000))

V_sph_approxs = []
V_sph_analyts = []
V_sph_d = 2.0
ds = []
for d in range(0, 200-1):
    V_sph_d *= Qapproxs[d]
    V_sph_approxs.append(V_sph_d)
    V_sph_analyts.append(V_sph(d+2))
    ds.append(d+2)

pylab.plot(ds, V_sph_analyts, 'k-', ds, V_sph_approxs, 'kx')
pylab.yscale('log')
pylab.xlabel('d')
pylab.ylabel(r'$V_{sph}(d)$')
pylab.title('Volume of unit sphere in d dimensions')
pylab.legend(['analytical', 'approximate'])
pylab.grid()
pylab.savefig('markov_sphere.png')
pylab.show()
