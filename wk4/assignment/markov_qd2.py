import math, random, pymp

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

def V_sph_approx(dim, n_trials):
    retval = 2.0
    for d in range(1, dim):
        retval *= markov_chain_qd(d+1, n_trials)
    return retval

print "{:<10} | {:<22} | {:<22} | {:<22} | {:<22}".format('n_trials', '<V_sph(20)>',
        'V_sph(20) (exact)', 'error', 'difference')
n_runs = 20
for n_trials in [1, 10, 100, 1000, 10000, 100000, 1000000]:
    V_sph_sum = 0.0
    V_sphsq_sum = 0.0
    V_sph_runs = pymp.shared.array((n_runs), dtype='float64')
    with pymp.Parallel(8) as p:
        for i in p.range(0, n_runs):
            V_sph_runs[i] = V_sph_approx(20, n_trials)
    
    for V_sph_approx_run in V_sph_runs:
        V_sph_sum += V_sph_approx_run
        V_sphsq_sum += V_sph_approx_run**2
    
    V_sph_avg = V_sph_sum / n_runs
    V_sphsq_avg = V_sphsq_sum / n_runs
    error_est = math.sqrt(abs(V_sphsq_avg - V_sph_avg**2)) / math.sqrt(n_runs)

    print "{:<10} | {:<22} | {:<22} | {:<22} | {:<22}".format(n_trials, 
            V_sph_avg, V_sph(20), error_est, abs(V_sph_avg - V_sph(20)))
