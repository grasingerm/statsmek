import os, pylab, random, math, cmath

def delx_dely(x, y):
    d_x = (x[0] - y[0]) % 1.0
    if d_x > 0.5: d_x -= 1.0
    d_y = (x[1] - y[1]) % 1.0
    if d_y > 0.5: d_y -= 1.0
    return d_x, d_y

def Psi_6(L, sigma):
    sum_vector = 0j
    for i in range(N):
        vector  = 0j
        n_neighbor = 0
        for j in range(N):
            if dist(L[i], L[j]) < 2.8 * sigma and i != j:
                n_neighbor += 1
                dx, dy = delx_dely(L[j], L[i])
                angle = cmath.phase(complex(dx, dy))
                vector += cmath.exp(6.0j * angle)
        if n_neighbor > 0:
            vector /= n_neighbor
        sum_vector += vector
    return sum_vector / float(N)

def dist(x,y):
    d_x = abs(x[0] - y[0]) % 1.0
    d_x = min(d_x, 1.0 - d_x)
    d_y = abs(x[1] - y[1]) % 1.0
    d_y = min(d_y, 1.0 - d_y)
    return math.sqrt(d_x**2 + d_y**2)

def is_valid_conf(L, sigma):
    N = len(L)
    for i in range(N-1):
        for j in range(i+1, N):
            if dist(L[i], L[j]) < 2.0 * sigma:
                return False
    return True

def rand_conf(N):
    L = []
    for k in range(N):
        L.append([random.uniform(0.0, 1.0), random.uniform(0.0, 1.0)])
    return L

def init_conf(N, eta, filename):
    sigma_sq = eta / (math.pi * N)
    sigma = math.sqrt(sigma_sq)
    if os.path.isfile(filename):
        f = open(filename, 'r')
        L = []
        for line in f:
            a, b = line.split()
            L.append([float(a), float(b)])
        f.close()
        print 'starting from file', filename
    else:
        N_sqrt = int(math.sqrt(N))
        two_delxy = 1.0 / N_sqrt
        delxy = 0.5 * two_delxy
        L = [[delxy + i * two_delxy, delxy + j * two_delxy] for i in range(N_sqrt) for j in range(N_sqrt)]
    return L, sigma, sigma_sq

def show_conf(L, sigma, title, fname):
    pylab.axes()
    for [x, y] in L:
        for ix in range(-1, 2):
            for iy in range(-1, 2):
                cir = pylab.Circle((x + ix, y + iy), radius=sigma,  fc='r')
                pylab.gca().add_patch(cir)
    pylab.axis('scaled')
    pylab.title(title)
    pylab.axis([0.0, 1.0, 0.0, 1.0])
    pylab.savefig(fname)
    pylab.show()
    pylab.close()

def save_conf(L, filename):
    f = open(filename, 'w')
    for a in L:
       f.write(str(a[0]) + ' ' + str(a[1]) + '\n')
    f.close()

N = 64
eta = 0.10 
filename = 'disk_configuration_order_study.txt'
L, sigma, sigma_sq = init_conf(N, eta, filename)
delta = 0.3 * sigma
n_steps = 10000
Psi_sum = 0.0
samples = 0
for steps in range(n_steps):
    a = random.choice(L)
    b = [a[0] + random.uniform(-delta, delta), a[1] + random.uniform(-delta, delta)]
    min_dist = min(dist(b, c) for c in L if c != a)
    if min_dist >= 2.0 * sigma:
        a[:] = [x % 1.0 for x in b]
    if steps % 100 == 0:
        Psi_sum += abs(Psi_6(L, sigma))
        samples += 1

#show_conf(L, sigma, 'Markov disks, periodic BCs', 'disk_configuration_N%i_eta%.2f.png' % (N, eta))
with open("order_study.csv",'a') as w:
    w.write(str(eta) + "," + str(Psi_sum / float(samples)) + '\n')
save_conf(L, filename)
