from math import ceil

import numpy as np

from utils.vec import max_abs_diff
from writers.one_d.data import write_data, write_statistics

u_bot = 1.
nu = 1.

h = 1.
dy = 0.01

r = 0.3

end_time = 0.3
out_time_step = 0.003

out_dir = f'.out/1d/couette_flow/free_surface/second_order_bc/nu={nu}, U={u_bot}, H={h}/dy={dy}, r={r}'

dt = r*dy*dy/nu

ny = ceil(h / dy) + 1

u = np.zeros(ny)
u1 = np.zeros(ny)

u[0] = u_bot

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
