from math import ceil

import numpy as np

from utils.vec import max_abs_diff
from utils.one_d.viscosity import linear_profile

from writers.one_d.data import write_data, write_statistics

qs = 1.

nus = 1.
nub = 0.01

h = 1.
ht = 0.3

dz_min = 0.01

r = 0.3

end_time = 0.3
out_time_step = 0.003

out_dir = f'.out/1d/couette_flow/wind_induced_current/uniform_dz/qs={qs}, nus={nus}, nub={nub}, ht={ht}, H={h}/dz={dz}, r={r}'

dt = r*dz_min*dz_min/nus

nz = ceil(h / dz_min)

dz = dz_min * np.ones(nz)

nu = np.array([linear_profile(i*dz, ht, nus, nub) for i in range(nz + 1)])

u = np.zeros(nz)
u1 = np.zeros(nz)

up = u.copy()

t = 0.
tn = out_time_step

n = 1
m = 0

stats = []

write_data(u, t, dy, m, out_dir)

m += 1

while t <= end_time:
    for i in range(1, ny-1):
        u1[i] = u[i] + r*(u[i+1] - 2*u[i] + u[i-1])

    u1[0] = u_bot
    u1[ny-1] = u[ny-1] + r*(-2*u[ny-1] + 2*u[ny-2])

    u, u1 = u1, u

    t += dt

    if t >= tn:
        write_data(u, t, dy, m, out_dir)

        max_diff = max_abs_diff(up, u)

        print(f'Write data in file t={t:.3f}, convergence of u={max_diff:.5f}')

        stats.append([n, t, max_diff])

        up = u.copy()

        m += 1

        tn += out_time_step

    n += 1

write_statistics(stats, out_dir)
