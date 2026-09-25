from math import ceil

import numpy as np

from slae.direct.tridiagonal import solve_uniform

from utils.vec import max_abs_diff
from writers.one_d.data import write_data, write_statistics

u_top = 1.
nu = 1.

h = 1.
dy = 0.01

r = 30

end_time = 0.3
out_time_step = 0.003

out_dir = f'.out/1d/couette_flow/implicit_euler/nu={nu}, U={u_top}, H={h}/dy={dy}, r={r}'

dt = r*dy*dy/nu

ny = ceil(h / dy) + 1

u = np.zeros(ny)
u1 = np.zeros(ny)

u[0] = 0
u[ny-1] = u_top

up = u.copy()

t = 0.
tn = out_time_step

n = 1
m = 0

stats = []

write_data(u, t, dy, m, out_dir)

m += 1

al = -r
al1 = 0

ac = 1 + 2*r
ac0 = 1
ac1 = 1

ar = -r
ar0 = 0

while t <= end_time:
    solve_uniform(al, al1, ac, ac0, ac1, ar, ar0, u, u1)

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
