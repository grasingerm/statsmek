import math, numpy
import matplotlib.pyplot as plt
import argparse

def_x_max = 5.0
def_nx = 100

parser = argparse.ArgumentParser(description='Matrix squaring quantum harmonic oscillator')
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
x_max = args.x_max
nx = args.nx
dx = 2.0 * x_max / (nx - 1)
x = [i * dx for i in range(-(nx - 1) / 2, nx / 2 + 1)]
beta = args.beta                     # actual value of beta (power of 2)

print '{0:12} {1:12} {2:12} {3}'.format('quartic', 'Z matrix', 'Z pert', \
                                        '|Z_matrix - Z_pert|')
for quartic in [0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5]:
    beta_tmp = args.beta_init            # initial value of beta (power of 2)
    cubic = -quartic
    rho = rho_anharmonic_trotter(x, beta_tmp, cubic, quartic)
    while beta_tmp < beta:
        rho = numpy.dot(rho, rho)
        rho *= dx
        beta_tmp *= 2.0

    Z_q = sum(rho[j, j] for j in range(nx + 1)) * dx
    Z_pert_q = 0
    try:
        Z_pert_q = Z_pert(cubic, quartic, beta, 50)
        print '{0:12.8f} {1:12.8f} {2:12.8f} {3:12.8f}'.format(quartic, Z_q, \
                Z_pert_q, abs(Z_q - Z_pert_q)) 
    except OverflowError:
        print '{0:12.8f} {1:12.8f} {2:12} {3:12}'.format(quartic, Z_q, 'Inf', \
                'Inf') 

