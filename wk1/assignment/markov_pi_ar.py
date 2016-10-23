import random, math, pylab

def markov_pi_ar(N, delta):
    x, y = 1.0, 1.0
    nac = 0
    for i in range(N):
        del_x, del_y = random.uniform(-delta, delta), random.uniform(-delta, delta)
        if abs(x + del_x) < 1.0 and abs(y + del_y) < 1.0:
            x, y = x + del_x, y + del_y
            nac = nac + 1
    return nac / float(N)

print "{:10} | {:20}".format("delta", "acceptance rate")
n_trials = 2 ** 12
for delta in [0.062, 0.125, 0.25, 0.5, 1.0, 2.0, 4.0]:
    print "{:10} | {:20}".format(delta, markov_pi_ar(n_trials, delta))

