"""
2D single-phase Darcy flow: steady-state verification of well models (Peaceman).

    (k h/mu) (d2p/dx2 + d2p/dy2) + q delta(x) delta(y) = 0

5-point finite volume scheme, well at the centre node (cell K), solved by SOR.
Boundary nodes are held at the exact infinite-domain line-source solution

    p(r) = q mu/(2 pi k h) * ln(R/r)          (p = 0 at the reference radius R)

so the computed p_K can be compared with the exact pressure. Well models:

    p_w = p_K + q mu/(2 pi k h) * ln(r_o/r_w)

    naive         r_o = r_w                       (p_w = p_K)
    equal area    r_o = sqrt(dx dy/pi)
    Peaceman      r_o = 0.14 sqrt(dx^2 + dy^2)    (= 0.198 dx for dx = dy)

Equivalent radius of the computed cell pressure:  r_eq = R exp(-p_K/c0),  c0 = q mu/(2 pi k h)
"""
from math import ceil, sqrt, cos, log, exp, pi

import numpy as np

from writers.two_d.data import write_data


def set_boundary_condition(p, nx, ny, pb):
    for i in range(nx):
        p[i*ny] = pb[i*ny]
        p[ny - 1 + i*ny] = pb[ny - 1 + i*ny]

    for j in range(ny):
        p[j] = pb[j]
        p[j + (nx - 1)*ny] = pb[j + (nx - 1)*ny]


def calc_max_abs_res(nx, ny, p, f):
    max_res = 0.

    for i in range(1, nx-1):
        for j in range(1, ny-1):
            pc = j + i*ny

            pl = j + (i - 1)*ny
            pr = j + (i + 1)*ny

            pd = j - 1 + i*ny
            pu = j + 1 + i*ny

            max_res = max(
                max_res,
                abs(p[pr] + p[pl] + b*b*(p[pu] + p[pd]) - 2*(1 + b*b)*p[pc] + f[pc]*dx*dx)
            )

    return max_res


md = 9.869233e-16                  # 1 mD in m^2

k = 100*md                         # permeability, m^2
h = 10.                            # reservoir thickness, m
mu = 1e-3                          # viscosity, Pa s
q = 200./86400                     # injection rate, m^3/s (q < 0 for a producer)
rw = 0.1                           # well radius, m

dx = 20.
dy = 20.

lx = 2000.                         # domain size, m (the well is in the centre)
ly = 2000.

r_ref = 1000.                      # radius where p = 0 in the exact solution, m

epsilon = 1e-8                     # tolerance on the residual relative to q mu/(k h)

out_dir = f'.out/2d/darcy_flow/well_model/steady_state/sor/dx={dx}, dy={dy}'

b = dx/dy

nx = ceil(lx / dx) + 1
ny = ceil(ly / dy) + 1

if nx % 2 == 0: nx += 1            # odd number of nodes -> a centre node exists
if ny % 2 == 0: ny += 1

iw = nx//2
jw = ny//2

pw_node = jw + iw*ny

c0 = q*mu/(2*pi*k*h)

# exact solution at the nodes (used on the boundary and for comparison)
p_ex = np.zeros(nx*ny)

for i in range(nx):
    for j in range(ny):
        r = sqrt(((i - iw)*dx)**2 + ((j - jw)*dy)**2)

        p_ex[j + i*ny] = c0*log(r_ref/r) if r > 0 else 0.

# source: discrete delta = 1/(dx dy) at the well node
f = np.zeros(nx*ny)
f[pw_node] = q*mu/(k*h)/(dx*dy)

l = (cos(pi/(nx-1)) + b*b*cos(pi/(ny-1)))/(1 + b*b)
kk = l*l
omega = 2*(1 - sqrt(1 - kk))/kk

p = np.zeros(nx*ny)

set_boundary_condition(p, nx, ny, p_ex)

scale = q*mu/(k*h)

max_res = calc_max_abs_res(nx, ny, p, f)/scale

m = 1

while max_res > epsilon:
    for i in range(1, nx - 1):
        for j in range(1, ny - 1):
            pc = j + i*ny

            pr = j + (i + 1)*ny
            pl = j + (i - 1)*ny

            pu = j + 1 + i*ny
            pd = j - 1 + i*ny

            p[pc] += omega*(
                p[pr] + p[pl] +
                b*b*(p[pu] + p[pd]) -
                2*(1 + b*b)*p[pc] +
                f[pc]*dx*dx
            )/2/(b*b + 1)

    max_res = calc_max_abs_res(nx, ny, p, f)/scale

    if m % 50 == 0:
        print(f'Iteration={m}, relative residual={max_res:.3e}')

    m += 1

print(f'Converged in {m - 1} iterations, omega={omega:.4f}')

# ----------------------------------------------------------------- well models
pk = p[pw_node]
pw_exact = c0*log(r_ref/rw)

r_eq = r_ref*exp(-pk/c0)

print(f'\nCell pressure p_K = {pk/1e5:.4f} bar, exact well pressure p_w = {pw_exact/1e5:.4f} bar')
print(f'Equivalent radius r_eq = {r_eq:.4f} m = {r_eq/dx:.5f} dx '
      f'(Peaceman: {0.14*sqrt(dx*dx + dy*dy)/dx:.5f} dx)\n')

models = {
    'naive, p_w = p_K': rw,
    'equal area, r_o = sqrt(dx dy/pi)': sqrt(dx*dy/pi),
    'Peaceman, r_o = 0.14 sqrt(dx^2 + dy^2)': 0.14*sqrt(dx*dx + dy*dy),
}

for name, ro in models.items():
    pw = pk + c0*log(ro/rw)

    print(f'{name:42s} r_o = {ro:8.3f} m   p_w = {pw/1e5:8.4f} bar   '
          f'error = {100*(pw - pw_exact)/pw_exact:+7.2f}%')

# pressure in bar
write_data(p/1e5, 0, dx, dy, 0, nx, ny, out_dir)
