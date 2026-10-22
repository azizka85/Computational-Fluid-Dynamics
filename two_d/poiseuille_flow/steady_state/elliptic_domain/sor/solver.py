from math import ceil, sqrt, cos, pi

import numpy as np

from slae.direct.tridiagonal import solve_uniform

from writers.two_d.data import write_data

g = 1.
nu = 1.

a = 2.
b = 1.

dx = 0.01
dy = 0.01

r = dx/dy

epsilon = 1e-6

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

def calc_max_abs_res(
    nx, ny,
    x0, y0,
    w
):
    max_res = 0.

    for i in range(1, nx-1):
        for j in range(1, ny-1):
            x = x0 + i*dx
            y = y0 + j*dy

            is_domain, is_boundary = point_in_domain(x, y)

            if is_domain and not is_boundary:
                p = j + i*ny

                pl = j + (i - 1)*ny
                pr = j + (i + 1)*ny

                pd = j - 1 + i*ny
                pu = j + 1 + i*ny

                max_res = max(
                    max_res,
                    abs(
                        (w[pr] - 2*w[p] + w[pl])/dx/dx +
                        (w[pu] - 2*w[p] + w[pd])/dy/dy -
                        g/nu
                    )
                )

    return max_res

out_dir = f'.out/2d/poiseuille_flow/steady_state/elliptic_domain/sor/dx={dx}, dy={dy}'

(x0, y0, l, h) = generate_domain()

nx = ceil(l / dx) + 1
ny = ceil(h / dy) + 1

f = (cos(pi/(nx-1)) + r*r*cos(pi/(ny-1)))/(1 + r*r)
k = f*f
omega = 2*(1 - sqrt(1 - k))/k

w = np.zeros(nx*ny)

max_res = calc_max_abs_res(nx, ny, x0, y0, w)

m = 1

while max_res > epsilon:
    for i in range(1, nx - 1):
        for j in range(1, ny - 1):
            x = x0 + i*dx
            y = y0 + j*dy

            is_domain, is_boundary = point_in_domain(x, y)

            if is_domain and not is_boundary:
                p = j + i*ny

                pl = j + (i - 1)*ny
                pr = j + (i + 1)*ny

                pd = j - 1 + i*ny
                pu = j + 1 + i*ny

                w[p] += omega*(
                    w[pr] + w[pl] +
                    r*r*(w[pu] + w[pd]) -
                    2*(1 + r*r)*w[p] -
                    dx*dx*g/nu
                )/2/(r*r + 1)

    max_res = calc_max_abs_res(nx, ny, x0, y0, w)

    print(f'Iteration={m}, convergence of w={max_res:.6f}')

    m += 1

write_data(w, 0, dx, dy, 0, nx, ny, out_dir)
