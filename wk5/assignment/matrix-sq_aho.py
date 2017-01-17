from __future__ import print_function
import math, numpy
import matplotlib.pyplot as plt
import argparse

def_x_max = 5.0
def_nx = 100

parser = argparse.ArgumentParser(description='Matrix squaring quantum harmonic oscillator')
parser.add_argument('--cubic', default=-1.0, type=float, help='Potential property')
parser.add_argument('--quartic', default=1.0, type=float, help='Potential property')
parser.add_argument('--x_max', default=5.0, type=float, help='Maximum x to discretize')
parser.add_argument('--nx', default=100, type=int, 
                    help='Number of discrete points')
parser.add_argument('--beta', default=4.0, type=float, help='1 / (kB * T)')
parser.add_argument('--beta_init', default=2.0 ** (-6), 
                    type=float, help='Initial beta, 1 / (kB * T)')
parser.add_argument('--filename', default='matrix_aho.png', type=str,
                    help='File name of plot')

def Energy_pert(n, cubic, quartic):
    return n + 0.5 - 15.0 / 4.0 * cubic ** 2 * (n ** 2 + n + 11.0 / 30.0) \
            + 3.0 / 2.0 * quartic * (n ** 2 + n + 1.0 / 2.0)

def Z_pert(cubic, quartic, beta, n_max):
    return sum(math.exp(-beta * Energy_pert(n, cubic, quartic)) \
               for n in range(n_max + 1))

def V(x, cubic, quartic):
    return x ** 2 / 2.0 + cubic * x ** 3 + quartic * x ** 4

# Free off-diagonal density matrix
def rho_free(x, xp, beta):
    return (math.exp(-(x - xp) ** 2 / (2.0 * beta)) /
            math.sqrt(2.0 * math.pi * beta))

# Harmonic density matrix in the Trotter approximation (returns the full matrix)
def rho_anharmonic_trotter(grid, beta, cubic, quartic):
    return numpy.array([[rho_free(x, xp, beta) * \
                         numpy.exp(-0.5 * beta * (V(x, cubic, quartic) + \
                         V(xp, cubic, quartic))) for x in grid] for xp in grid])

args = parser.parse_args()
cubic = args.cubic
quartic = args.quartic
x_max = args.x_max
nx = args.nx
dx = 2.0 * x_max / (nx - 1)
x = [i * dx for i in range(-(nx - 1) / 2, nx / 2 + 1)]
beta_tmp = args.beta_init            # initial value of beta (power of 2)
beta = args.beta                     # actual value of beta (power of 2)
rho = rho_anharmonic_trotter(x, beta_tmp, cubic, quartic)

# Print initial density matrix to see that it is approximately diagonal
print('density matrix at beta = ', beta_tmp)
for i in range(nx + 1):
    for j in range(nx + 1):
        print("{0:.3f}".format(rho[i, j]), end=' ')
    print(end='\n')

while beta_tmp < beta:
    rho = numpy.dot(rho, rho)
    rho *= dx
    beta_tmp *= 2.0
    print('beta: %s -> %s' % (beta_tmp / 2.0, beta_tmp))

Z = sum(rho[j, j] for j in range(nx + 1)) * dx
pi_of_x = [rho[j, j] / Z for j in range(nx + 1)]
f = open('data_aharm_matrixsquaring_beta' + str(beta) + '.dat', 'w')
for j in range(nx + 1):
    f.write(str(x[j]) + ' ' + str(rho[j, j] / Z) + '\n')
f.close()

plt.plot(x, pi_of_x, 'x')
plt.xlabel('x')
plt.ylabel(r'$\rho(x)$')
plt.title(r'Matrix squaring: anharmonic, $\beta = ' + str(beta)
           + r'$, cubic = ' + str(cubic) + ', quartic = ' + str(quartic))
plt.grid()
plt.savefig(args.filename)
plt.show()
