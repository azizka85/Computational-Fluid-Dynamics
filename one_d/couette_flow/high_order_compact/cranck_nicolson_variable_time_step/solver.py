from math import ceil

import numpy as np

from slae.direct.tridiagonal import solve_uniform

from utils.vec import max_abs_diff
from utils.time import adjust_time_step

from writers.one_d.data import write_data, write_statistics

b = 1.1

u_top = 1.
nu = 1.

h = 1.
dy = 0.01

r = 0.3

end_time = 0.3
out_time_step = 0.003

out_dir = f'.out/1d/couette_flow/high_order_compact/cranck_nicolson_variable_time_step/nu={nu}, U={u_top}, H={h}/dy={dy}, r={r}'

cr = r
dt = cr*dy*dy/nu

ny = ceil(h / dy) + 1

u = np.zeros(ny)
d = np.zeros(ny)

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

while t <= end_time:
    al = 1 - 6*cr
    al1 = 0

    ac = 10 + 12*cr
    ac0 = 1
    ac1 = 1

    ar = 1 - 6*cr
    ar0 = 0

    d[0] = u[0]
    d[ny-1] = u[ny-1]

    for i in range(1, ny-1):
        d[i] = (1 + 6*cr)*u[i+1] + (10 - 12*cr)*u[i] + (1 + 6*cr)*u[i-1]

    solve_uniform(al, al1, ac, ac0, ac1, ar, ar0, d, u)

    t += dt

    cdt = adjust_time_step(b, t, dt, out_time_step)

    cr = cr * cdt / dt
    dt = cdt

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
