from math import ceil, sin, pi

import numpy as np

from utils.vec import max_abs_diff

from writers.one_d.data import write_statistics
from writers.two_d.data import write_data

def set_boundary_condition(w, nx, ny):
    for i in range(nx):
        w[i*ny] = 0.
        w[ny - 1 + i*ny] = 0

    for j in range(ny):
        w[j] = 0
        w[j + (nx - 1)*ny] = 0

nu = 1.

dx = 0.01
dy = 0.01

r = 0.15

end_time = 0.3
out_time_step = 0.003

out_dir = f'.out/2d/poiseuille_flow/explicit_euler/nu={nu}/dx={dx}, dy={dy}, r={r}'

dt = min(r*dx*dx/nu, r*dy*dy/nu)

rx = nu*dt/dx/dx
ry = nu*dt/dy/dy

nx = ceil(1 / dx) + 1
ny = ceil(1 / dy) + 1

w = np.zeros(nx*ny)
w1 = np.zeros(nx*ny)

for i in range(nx):
    for j in range(ny):
        x = i*dx
        y = j*dy

        p = j + i*ny

        w[p] = sin(pi*x) * sin(pi*y)

set_boundary_condition(w, nx, ny)

wp = w.copy()

t = 0
tn = out_time_step

n = 1
m = 0

stats = []

write_data(w, t, dx, dy, m, nx, ny, out_dir)

m += 1

while t <= end_time:
    for i in range(1, nx - 1):
        for j in range(1, ny - 1):
            p = j + i*ny

            pr = j + (i + 1)*ny
            pl = j + (i - 1)*ny

            pu = j + 1 + i*ny
            pd = j - 1 + i*ny

            w1[p] = w[p] + rx*(w[pr] - 2*w[p] + w[pl]) + ry*(w[pu] - 2*w[p] + w[pd])

    set_boundary_condition(w, nx, ny)

    w, w1 = w1, w

    t += dt

    if t >= tn:
        write_data(w, t, dx, dy, m, nx, ny, out_dir)

        max_diff = max_abs_diff(wp, w)

        print(f'Write data in file t={t:.3f}, convergence of w={max_diff:.5f}')

        stats.append([n, t, max_diff])

        wp = w.copy()

        m += 1

        tn += out_time_step

    n += 1

write_statistics(stats, out_dir)
