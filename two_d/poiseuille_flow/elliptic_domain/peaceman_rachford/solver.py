from math import ceil, cos, pi

import numpy as np

from slae.direct.tridiagonal import solve_uniform

from utils.vec import max_abs_diff

from writers.one_d.data import write_statistics
from writers.two_d.data import write_data

a = 2.
b = 1.

dx = 0.01
dy = 0.01

def generate_domain():
    return (-a, -b, 2*a, 2*b)

def point_in_domain(x, y):
    if x*x/a/a + y*y/b/b <= 1:
        xl = x - dx
        xr = x + dx

        yd = y - dy
        yu = y + dy

        if (xl*xl/a/a + y*y/b/b > 1 or
            xr*xr/a/a + y*y/b/b > 1 or
            x*x/a/a + yd*yd/b/b > 1 or
            x*x/a/a + yu*yu/b/b > 1):

            return (True, True)

        return (True, False)

    return (False, False)

nu = 1.

r = 0.15

end_time = 0.3
out_time_step = 0.003

out_dir = f'.out/2d/poiseuille_flow/elliptic_domain/peaceman_rachford/nu={nu}/dx={dx}, dy={dy}/r={r}'

dt = min(r*dx*dx/nu, r*dy*dy/nu)

rx = nu*dt/dx/dx
ry = nu*dt/dy/dy

(x0, y0, l, h) = generate_domain()

nx = ceil(l / dx) + 1
ny = ceil(h / dy) + 1

w = np.zeros(nx*ny)
w1 = np.zeros(nx*ny)

for i in range(nx):
    for j in range(ny):
        x = x0 + i*dx
        y = y0 + j*dy

        p = j + i*ny

        (is_domain, _) = point_in_domain(x, y)

        if is_domain:
            w[p] = cos(pi*x/4) * cos(pi*y/2)
        else:
            w[p] = 0.

wp = w.copy()

t = 0
tn = out_time_step

n = 1
m = 0

stats = []

write_data(w, t, dx, dy, m, nx, ny, out_dir)

m += 1

vx = np.zeros(nx)
vy = np.zeros(ny)

wx = np.zeros(nx)
wy = np.zeros(ny)

while t <= end_time:
    for j in range(1, ny-1):
        vx[0] = w[j]

        for i in range(1, nx-1):
            p = j + i*ny

            pu = j + 1 + i*ny
            pd = j - 1 + i*ny

            vx[i] = 0.5 * ry * w[pu] + (1 - ry)*w[p] + 0.5 * ry * w[pd]

        vx[nx-1] = w[j + (nx - 1)*ny]

        solve_uniform(-0.5*rx, 0, 1 + rx, 1, 1, -0.5*rx, 0, vx, wx)

        for i in range(nx):
            p = j + i*ny

            w1[p] = wx[i]

    for i in range(1, nx - 1):
        vy[0] = w1[i*ny]

        for j in range(1, ny - 1):
            vy[j] = 0.5*rx*w1[j+(i+1)*ny] + (1 - rx)*w1[j+i*ny] + 0.5*rx*w1[j+(i-1)*ny]

        vy[ny-1] = w1[ny-1+i*ny]

        solve_uniform(-0.5*ry, 0, 1 + ry, 1, 1, -0.5*ry, 0, vy, wy)

        for j in range(ny):
            p = j + i*ny

            w[p] = wy[j]

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
