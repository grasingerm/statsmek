import math, random

nsamples = 1000
x_samples = [0]*nsamples
y_samples = [0]*nsamples
for k in range(nsamples):
    upsilon = random.uniform(0.0, 1.0)
    phi = random.uniform(0.0, 2.0*math.pi)
    psi = -math.log(upsilon)
    r = math.sqrt(2*psi)
    x_samples[k] = r * math.cos(phi)
    y_samples[k] = r * math.sin(phi)
    print x_samples[k], y_samples[k]
