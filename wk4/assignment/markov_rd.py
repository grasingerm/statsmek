import math, random
import numpy as np
import matplotlib.pyplot as plt

def markov_chain_sphere(d, n_trials, delta=0.1):
    rs = []
    xs = [0.0] * d
    rsq = 0.0
    for i in range(n_trials):

        # trial move
        k = random.randint(0, d - 1)
        xk = xs[k]
        new_xk = xk + random.uniform(-delta, delta)
        new_rsq = rsq + new_xk*new_xk - xk*xk 
        if new_rsq < 1.0:
            xs[k] = new_xk
            rsq = new_rsq

        rs.append(math.sqrt(rsq))

    return rs

rs = markov_chain_sphere(20, 10000000)
plt.hist(rs, 100, normed=True)
xs = np.linspace(0, 1, 100)
#ys = [4 * r**3 for r in xs]
ys = [20 * r**19 for r in xs]
plt.plot(xs, ys)
plt.xlabel('r')
plt.ylabel('frequency')
plt.title('Markov sampling: radius of sphere histogram (20d)')
plt.grid()
plt.savefig('markov_20d.png')
plt.show()
