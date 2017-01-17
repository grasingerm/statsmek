import random, math, sys
import numpy as np
import matplotlib.pyplot as plt

"""
Acceptance probability of trial move for a quantum harmonic oscillator

    \param   x           Current position
    \param   x_prime     Trial position
    \return              Acceptance probability
"""
def p_accept(x, x_prime):
    if abs(x) >= abs(x_prime):
        return 1.0
    else: 
        return math.exp(x**2 - x_prime**2)

"""
Ground-state normalized probability, i.e. wave function squared

    \param   x           Position
    \return              Probability of being at x
"""
__SQRT_PI__ = math.sqrt(math.pi)
def psi_0_sq(x):
    return math.exp(-(x**2)) / __SQRT_PI__

# Main program
nsteps = 100000
if len(sys.argv) >= 2: # Let number of steps be a (optional) command-line argument
    nsteps = int(sys.argv[1])

delta = 0.5
if len(sys.argv) == 3: # Let the step size be the second (optional) argument
    delta = float(sys.argv[2])

x = 0.0
data = [0] * nsteps # Allocate list for histogram data

# Main loop
for k in range(nsteps):

    # Notify user every 10% of steps completed
    if (k + 1) % (nsteps / 10) == 0:
        print k + 1, " / ", nsteps

    # Sample
    x_new = x + random.uniform(-delta, delta)
    if random.uniform(0.0, 1.0) < p_accept(x, x_new): 
        x = x_new 
    data[k] = x

# Plot data
plt.hist(data, 100, normed=True)
max_abs = max([abs(min(data)), max(data)])
xs = np.linspace(-max_abs, max_abs, 100)
ys = [psi_0_sq(x) for x in xs]
plt.plot(xs, ys, linewidth=4.0)
plt.xlabel('x')
plt.ylabel('p(x)')
plt.legend(['Analytical', 'Markov chain'])
plt.title('Markov chain sampling: 1D harmonic oscillator')
plt.grid()
plt.savefig('markov_ho.png')
plt.show()
