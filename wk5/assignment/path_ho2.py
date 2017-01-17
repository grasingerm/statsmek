import math, random, numpy
import matplotlib.pyplot as plt

def rho_free(x, y, beta):    # free off-diagonal density matrix
    return math.exp(-(x - y) ** 2 / (2.0 * beta))

def read_file(filename):
    list_x = []
    list_y = []
    with open(filename) as f:
        for line in f:
            x, y = line.split()
            list_x.append(x)
            list_y.append(y)
    return list_x, list_y

import argparse

parser = argparse.ArgumentParser(description="Path-integral Monte Carlo")
parser.add_argument('--n_steps', default=1000000, type=int, 
                    help='Number of steps to run')
parser.add_argument('--delta', default=1.0, type=float, help='Maximum step size')
parser.add_argument('--beta', default=4.0, type=float, help='1 / (kB * T)')
parser.add_argument('--N', default=10, type=int, help='Number of slices')
parser.add_argument('--xlim', type=str, help='Limits of x-axis')

args = parser.parse_args()
beta = args.beta
N = args.N                                        # number of slices
dtau = beta / N
delta = args.delta                                # maximum displacement on one slice
n_steps = args.n_steps                            # number of Monte Carlo steps
x = [0.0] * N                                     # initial path
hist_data = []
for step in range(n_steps):
    k = random.randint(0, N - 1)                  # random slice
    knext, kprev = (k + 1) % N, (k - 1) % N       # next/previous slices
    x_new = x[k] + random.uniform(-delta, delta)  # new position at slice k
    old_weight  = (rho_free(x[knext], x[k], dtau) *
                   rho_free(x[k], x[kprev], dtau) *
                   math.exp(-0.5 * dtau * x[k] ** 2))
    new_weight  = (rho_free(x[knext], x_new, dtau) *
                   rho_free(x_new, x[kprev], dtau) *
                   math.exp(-0.5 * dtau * x_new ** 2))
    if random.uniform(0.0, 1.0) < new_weight / old_weight:
        x[k] = x_new
    if step % 10 == 0:
        for xk in x:
            hist_data.append(xk)

plt.hist(hist_data, 100, normed=True)
matrix_sq_data = read_file('data_harm_matrixsquaring_beta4.0.dat')
plt.plot(matrix_sq_data[0], matrix_sq_data[1], 'ro')
if args.xlim != None:
    xlim = eval(args.xlim)
    plt.xlim(*xlim)
plt.legend(['Matrix squaring', 'Path integral'])
plt.xlabel('x')
plt.ylabel(r'$\rho(x)$')
plt.title(r'Path integral compared to matrix squaring, $\beta = ' + str(beta) +
          r'$')
plt.show()
