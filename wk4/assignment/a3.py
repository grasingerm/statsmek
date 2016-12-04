import math, random

x, y = 0.0, 0.0
delta = 0.1
n_trials = 10000000
n_hits = 0
for i in range(n_trials):
    del_x, del_y = random.uniform(-delta, delta), random.uniform(-delta, delta)
    if (x + del_x)**2 + (y + del_y)**2 < 1.0:
        x, y = x + del_x, y + del_y
    z = random.uniform(-1.0, 1.0)
    if x**2 + y**2 + z**2 < 1.0: n_hits += 1

q3 = 2 * n_hits / float(n_trials)

x, y, z = 0.0, 0.0, 0.0
delta = 0.1
n_trials = 10000000
n_hits = 0
for i in range(n_trials):
    del_x, del_y, del_z = random.uniform(-delta, delta), random.uniform(-delta, delta), random.uniform(-delta, delta)
    if (x + del_x)**2 + (y + del_y)**2 + (z + del_z)**2 < 1.0:
        x, y, z = x + del_x, y + del_y, z + del_z
    alpha = random.uniform(-1.0, 1.0)
    if x**2 + y**2 + z**2 + alpha**2 < 1.0: n_hits += 1

q4 = 2 * n_hits / float(n_trials)

print "pi^2 / 2 = ", math.pi**2 / 2.0
print "V_sph(3) * Q_4 = ", 4 * math.pi / 3.0 * q4
print "V_sph(2) * Q_3 * Q_4 = ", math.pi * q3 * q4
