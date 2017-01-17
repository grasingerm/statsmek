import random, math, sys
import numpy as np
import matplotlib.pyplot as plt

"""
Analytical solution
"""
def pi_quant(x, beta):
    return (math.sqrt(math.tanh(beta / 2) / math.pi) * 
            math.exp(-(x ** 2) * math.tanh(beta / 2)))

"""
Normalized probability

    \param   x           Position
    \param   n           Index of eigenstate, e.g. groundstate is 0
    \return              Probability of being at x, in state n
"""
def psi_n_sq(x, n):
    if n == -1:
        return 0.0
    else:
        psi = [math.exp(-x ** 2 / 2.0) / math.pi ** 0.25]
        psi.append(math.sqrt(2.0) * x * psi[0])
        for k in range(2, n + 1):
            psi.append(math.sqrt(2.0 / k) * x * psi[k - 1] -
                       math.sqrt((k - 1.0) / k) * psi[k - 2])
        return psi[n] ** 2

"""
Acceptance probability of trial move for a quantum harmonic oscillator

    \param   x           Current position
    \param   n           Current eigenstate
    \param   x_prime     Trial position
    \param   n_prime     Trial eigenstate
    \param   beta        1 / (kB * T)
    \return              Acceptance probability
"""
def p_accept(x, n, x_prime, n_prime, beta):
    E = n + 0.5
    E_prime = n_prime + 0.5
    return min(1, psi_n_sq(x_prime, n_prime) / psi_n_sq(x, n) * 
                  math.exp(-beta * (E_prime - E)))

import argparse

parser = argparse.ArgumentParser(description="Markov chain quantum harmonic oscillator")
parser.add_argument('--nsteps', default=100000, type=int, 
                    help='Number of steps to run')
parser.add_argument('--delta', default=0.5, type=float, help='Maximum step size')
parser.add_argument('--beta', default=1.0, type=float, help='1 / (kB * T)')
parser.add_argument('--filename', default='markov_ho.png', type=str,
                    help='File name of plot')

# Main program
args = parser.parse_args()
nsteps = args.nsteps
beta = args.beta
delta = args.delta

x = 0.0
n = 0
data = [0] * nsteps # Allocate list for histogram data

# Main loop
for k in range(nsteps):

    # Notify user every 10% of steps completed
    if (k + 1) % (nsteps / 10) == 0:
        print k + 1, " / ", nsteps

    # Sample position
    x_prime = x + random.uniform(-delta, delta)
    if random.uniform(0.0, 1.0) < p_accept(x, n, x_prime, n, beta): 
        x = x_prime 

    n_prime = n + random.randint(-1, 1)
    if random.uniform(0.0, 1.0) < p_accept(x, n, x, n_prime, beta): 
        n = n_prime 

    data[k] = x

# Plot data
plt.hist(data, 100, normed=True)
max_abs = max([abs(min(data)), max(data)])
xs = np.linspace(-max_abs, max_abs, 100)
y1s = [pi_quant(x, beta) for x in xs]
y2s = ([math.sqrt(beta / (2 * math.pi)) * math.exp(-beta * x ** 2 / 2.0)
        for x in xs])
plt.plot(xs, y1s, linewidth=4.0, linestyle='-')
plt.plot(xs, y2s, linewidth=4.0, linestyle='--')
plt.xlabel('x')
plt.ylabel('p(x)')
plt.legend(['Analytical, Quantum', 'Analytical, Classical', 'Markov chain'])
plt.title(r'Markov chain sampling: 1D harmonic oscillator, $\beta = ' + str(beta)
           + r'$')
plt.grid()
plt.savefig(args.filename)
plt.show()
