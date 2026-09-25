from math import ceil

import numpy as np

from utils.vec import max_abs_diff
from utils.one_d.viscosity import linear_profile

from writers.one_d.data import write_non_uniform_data, write_statistics

qs = 1.

kb = 0.002

nus = 1.
nub = 0.5

h = 1.
ht = 0.3

dz_min = 0.1

r = 0.3

end_time = 0.3
out_time_step = 0.003

out_dir = f'.out/1d/couette_flow/wind_induced_current/uniform_dz/qs={qs}, nus={nus}, nub={nub}, ht={ht}, H={h}/dz={dz_min}, r={r}'

dt = r*dz_min*dz_min/nus

nz = ceil(h / dz_min)

dz = dz_min * np.ones(nz)

nu = np.array([linear_profile(i*dz_min, ht, nus, nub) for i in range(nz + 1)])

u = np.zeros(nz)
u1 = np.zeros(nz)

up = u.copy()

t = 0.
tn = out_time_step

n = 1
m = 0

stats = []

write_non_uniform_data(u, t, dz, m, out_dir)

m += 1

while t <= end_time:
    for i in range(1, nz-1):
        u1[i] = u[i] + 2*dt*(
            nu[i+1]*(u[i+1] - u[i])/(dz[i+1] + dz[i]) -
            nu[i]*(u[i] - u[i-1])/(dz[i] + dz[i-1])
        )/dz[i]

    u1[0] = u[0] + 2*nu[1]*dt*(u[1] - u[0])/dz[0]/(dz[1] + dz[0]) - qs*dt/dz[0]
    u1[nz-1] = u[nz-1] + kb*u[nz-1]*dt/dz[nz-1] -\
                2*nu[nz-1]*dt*(u[nz-1] - u[nz-2])/dz[nz-1]  /(dz[nz-1] + dz[nz-2])

    u, u1 = u1, u

    t += dt

    if t >= tn:
        write_non_uniform_data(u, t, dz, m, out_dir)

        max_diff = max_abs_diff(up, u)

        print(f'Write data in file t={t:.3f}, convergence of u={max_diff:.5f}')

        stats.append([n, t, max_diff])

        up = u.copy()

        m += 1

        tn += out_time_step

    n += 1

write_statistics(stats, out_dir)
