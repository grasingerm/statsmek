import math, random, pylab

def levy_harmonic_path(k):
    x = [random.gauss(0.0, 1.0 / math.sqrt(2.0 * math.tanh(k * beta / 2.0)))]
    if k == 2:
        Ups1 = 2.0 / math.tanh(beta)
        Ups2 = 2.0 * x[0] / math.sinh(beta)
        x.append(random.gauss(Ups2 / Ups1, 1.0 / math.sqrt(Ups1)))
    return x[:]

def rho_harm_1d(x, xp, beta):
    Upsilon_1 = (x + xp) ** 2 / 4.0 * math.tanh(beta / 2.0)
    Upsilon_2 = (x - xp) ** 2 / 4.0 / math.tanh(beta / 2.0)
    return math.exp(- Upsilon_1 - Upsilon_2)

def z(beta):
    return 1.0 / (1.0 - math.exp(- beta))

def pi_two_bosons(x, beta):
    pi_x_1 = math.sqrt(math.tanh(beta / 2.0)) / math.sqrt(math.pi) *\
             math.exp(-x ** 2 * math.tanh(beta / 2.0))
    pi_x_2 = math.sqrt(math.tanh(beta)) / math.sqrt(math.pi) *\
             math.exp(-x ** 2 * math.tanh(beta))
    weight_1 = z(beta) ** 2 / (z(beta) ** 2 + z(2.0 * beta))
    weight_2 = z(2.0 * beta) / (z(beta) ** 2 + z(2.0 * beta))
    pi_x = pi_x_1 * weight_1 + pi_x_2 * weight_2
    return pi_x

betas = []
p1s = []
p2s = []
nsteps = 500000

for beta in [0.1, 0.5, 1.25, 2.5, 5.0]:
    n1cycle = 0
    n2cycle = 0
    low = levy_harmonic_path(2)
    high = low[:]
    for step in xrange(nsteps):
        # move 1
        if low[0] == high[0]:
            k = random.choice([0, 1])
            low[k] = levy_harmonic_path(1)[0]
            high[k] = low[k]
        else:
            low[0], low[1] = levy_harmonic_path(2)
            high[1] = low[0]
            high[0] = low[1]
        # move 2
        weight_old = (rho_harm_1d(low[0], high[0], beta) *
                      rho_harm_1d(low[1], high[1], beta))
        weight_new = (rho_harm_1d(low[0], high[1], beta) *
                      rho_harm_1d(low[1], high[0], beta))
        if random.uniform(0.0, 1.0) < weight_new / weight_old:
            high[0], high[1] = high[1], high[0]

        if low[0] == high[1] and high[0] == low[1]:
            n1cycle += 1
        else:
            n2cycle += 1
   
    betas.append(beta)
    p1s.append(float(n1cycle) / float(nsteps))
    p2s.append(float(n2cycle) / float(nsteps))

nx = 200
x_min, x_max = min(betas), max(betas)
dx = (x_max - x_min) / nx
list_x = [dx * i + x_min for i in range(nx)]
list_p2 = [z(beta) ** 2 / (z(beta) ** 2 + z(2.0 * beta)) for beta in list_x]
list_p1 = [z(2.0 * beta) / (z(beta) ** 2 + z(2.0 * beta)) for beta in list_x]
pylab.scatter(betas, p1s, marker="<", label="p1 cycle (MC)")
pylab.scatter(betas, p2s, marker="o", label="p2 cycle (MC)")
pylab.plot(list_x, list_p1, label="p1 analytical")
pylab.plot(list_x, list_p2, label="p2 analytical")
pylab.legend()
pylab.show()
