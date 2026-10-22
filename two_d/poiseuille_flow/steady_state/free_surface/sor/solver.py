from math import ceil, sqrt, cos, sin, pi

import numpy as np

from writers.two_d.data import write_data

def set_boundary_condition(w, nx, ny):
    for i in range(nx):
        w[i*ny] = 0.

    for j in range(ny):
        w[j] = 0
        w[j + (nx - 1)*ny] = 0

def calc_max_abs_res(nx, ny, w):
    max_res = 0.

    for i in range(1, nx-1):
        for j in range(1, ny-1):
            p = j + i*ny

            pl = j + (i - 1)*ny
            pr = j + (i + 1)*ny

            pd = j - 1 + i*ny
            pu = j + 1 + i*ny

            max_res = max(
                max_res,
                abs(
                    (w[pr] - 2*w[p] + w[pl])/dx/dx +
                    (w[pu] - 2*w[p] + w[pd])/dy/dy +
                    5*pi*pi*sin(pi*i*dx)*sin(pi*j*dy/2)/4
                )
            )

    return max_res

dx = 0.01
dy = 0.01

b = dx/dy

epsilon = 1e-6

out_dir = f'.out/2d/poiseuille_flow/steady_state/free_surface/sor/dx={dx}, dy={dy}'

nx = ceil(1 / dx) + 1
ny = ceil(1 / dy) + 1

l = (cos(pi/(nx-1)) + b*b*cos(pi/(ny-1)))/(1 + b*b)
k = l*l
omega = 2*(1 - sqrt(1 - k))/k

w = np.zeros(nx*ny)

for i in range(nx):
    for j in range(ny):
        x = i*dx
        y = j*dy

        p = j + i*ny

        w[p] = 0

set_boundary_condition(w, nx, ny)

max_res = calc_max_abs_res(nx, ny, w)

m = 1

while max_res > epsilon:
    for i in range(1, nx - 1):
        for j in range(1, ny - 1):
            p = j + i*ny

            pr = j + (i + 1)*ny
            pl = j + (i - 1)*ny

            pu = j + 1 + i*ny
            pd = j - 1 + i*ny

            w[p] += omega*(
                w[pr] + w[pl] +
                b*b*(w[pu] + w[pd]) -
                2*(1 + b*b)*w[p] +
                5*dx*dx*pi*pi*sin(pi*i*dx)*sin(pi*j*dy/2)/4
            )/2/(b*b + 1)

        w[ny-1+i*ny] += omega*(
            w[ny-1+(i+1)*ny] + w[ny-1+(i-1)*ny] +
            2*b*b*w[ny-2+i*ny] -
            2*(1 + b*b)*w[ny-1+i*ny] +
            5*dx*dx*pi*pi*sin(pi*i*dx)/4
        )/2/(b*b + 1)

    max_res = calc_max_abs_res(nx, ny, w)

    print(f'Iteration={m}, convergence of w={max_res:.6f}')

    m += 1

write_data(w, 0, dx, dy, m, nx, ny, out_dir)
