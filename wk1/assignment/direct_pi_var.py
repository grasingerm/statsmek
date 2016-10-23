import random, math
n_trials = 400000
n_hits = 0
sum_obs = 0.0
sum_obs_sqr = 0.0
for iter in range(n_trials):
    x, y = random.uniform(-1.0, 1.0), random.uniform(-1.0, 1.0)
    if x**2 + y**2 < 1.0:
        n_hits += 1
        sum_obs = sum_obs + 4.0
        sum_obs_sqr = sum_obs_sqr + 16.0

exp_obs = sum_obs / n_trials
exp_obs_sqr = sum_obs_sqr / n_trials
var = exp_obs_sqr - exp_obs**2

print 4.0 * n_hits / float(n_trials), math.sqrt(var)
