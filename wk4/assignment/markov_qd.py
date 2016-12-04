import math, random

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

print "Q(4), analytical:  ", Q(4)
print "Q(4), approximate: ", markov_chain_qd(4, 10000000)
print ""
print "Q(200), analytical:  ", Q(200)
print "Q(200), approximate: ", markov_chain_qd(200, 10000000)
