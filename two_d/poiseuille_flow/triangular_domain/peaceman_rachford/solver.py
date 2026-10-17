from math import ceil, sqrt, exp, cos, sin, pi

import numpy as np

from slae.direct.tridiagonal import solve_uniform

from utils.vec import max_abs_diff

from writers.one_d.data import write_statistics
from writers.two_d.data import write_data

a = 1.

dx = 0.01
dy = 0.01

def generate_domain():
    return (-a/2, -a, a, a)

def is_point_in_domain(x, y):
    return y <= 0 and y >= -a*sqrt(3)/2 + sqrt(3)*x and y >= -a*sqrt(3)/2 - sqrt(3)*x

def point_in_domain(x, y):
    if is_point_in_domain(x, y):
        xl = x - dx
        xr = x + dx

        yd = y - dy
        yu = y + dy

        if (not is_point_in_domain(xl, y) or
            not is_point_in_domain(xr, y) or
            not is_point_in_domain(x, yd) or
            not is_point_in_domain(x, yu)):

            return (True, True)

        return (True, False)

    return (False, False)

nu = 1.

r = 0.15

end_time = 0.3
out_time_step = 0.003

out_dir = f'.out/2d/poiseuille_flow/triangular_domain/peaceman_rachford/nu={nu}/dx={dx}, dy={dy}/r={r}'

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
            w[p] = cos(pi*x) * sin(pi*y)
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

while t <= end_time:
    for i in range(nx):
        (is_domain, _) = point_in_domain(x0 + i*dx, y0 + (ny - 1)*dy)

        if is_domain:
            x = x0 + i*dx
            y = y0 + (ny - 1)*dy

            wb = cos(pi*x)*sin(pi*y)*exp(-2*pi*pi*nu*t)
            w1b = cos(pi*x)*sin(pi*y)*exp(-2*pi*pi*nu*(t + dt))

            p = ny - 1 + i*ny

            w1[p] = (wb + w1b)/2

    for j in range(ny - 1):
        vx = []

        for i in range(nx):
            (is_domain, is_boundary) = point_in_domain(x0 + i*dx, y0 + j*dy)

            if is_domain:
                if is_boundary:
                    x = x0 + i*dx
                    y = y0 + j*dy

                    wb = cos(pi*x)*sin(pi*y)*exp(-2*pi*pi*nu*t)
                    w1b = cos(pi*x)*sin(pi*y)*exp(-2*pi*pi*nu*(t + dt))

                    vx.append((wb + w1b)/2)
                else:
                    p = j + i*ny

                    pu = j + 1 + i*ny
                    pd = j - 1 + i*ny

                    vx.append(0.5*ry*w[pu] + (1 - ry)*w[p] + 0.5*ry*w[pd])

        wx = np.zeros(len(vx))

        solve_uniform(-0.5*rx, 0, 1 + rx, 1, 1, -0.5*rx, 0, vx, wx)

        k = 0

        for i in range(nx):
            (is_domain, _) = point_in_domain(x0 + i*dx, y0 + j*dy)

            if is_domain:
                p = j + i*ny

                w1[p] = wx[k]

                k += 1

    for i in range(nx):
        vy = []

        for j in range(ny):
            (is_domain, is_boundary) = point_in_domain(x0 + i*dx, y0 + j*dy)

            if is_domain:
                if is_boundary:
                    x = x0 + i*dx
                    y = y0 + j*dy

                    vy.append(cos(pi*x)*sin(pi*y)*exp(-2*pi*pi*nu*(t + dt)))
                else:
                    p = j + i*ny

                    pr = j + (i + 1)*ny
                    pl = j + (i - 1)*ny

                    vy.append(0.5*rx*w1[pr] + (1 - rx)*w1[p] + 0.5*rx*w1[pl])

        wy = np.zeros(len(vy))

        solve_uniform(-0.5*ry, 0, 1 + ry, 1, 1, -0.5*ry, 0, vy, wy)

        k = 0

        for j in range(ny):
            (is_domain, _) = point_in_domain(x0 + i*dx, y0 + j*dy)

            if is_domain:
                p = j + i*ny

                w[p] = wy[k]

                k += 1

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
