"""
2D single-phase Darcy flow: transient test of well models against the Theis solution.

    phi c_t h dp/dt = (k h/mu) (d2p/dx2 + d2p/dy2) + q delta(x) delta(y)

5-point finite volumes, explicit Euler (r = eta dt/dx^2, stable for rx + ry <= 1/2),
well at the centre node. Theis solution (infinite domain, line source):

    p(r, t) - p0 = q mu/(4 pi k h) * E1(r^2/(4 eta t)),   eta = k/(phi mu c_t)

bc = 'theis'   : boundary nodes follow the Theis solution (infinite domain)
bc = 'initial' : boundary nodes stay at p0 (constant-pressure boundary)

Well models:  p_w = p_K + q mu/(2 pi k h) * ln(r_o/r_w)
    naive r_o = r_w,  equal area r_o = sqrt(dx dy/pi),  Peaceman r_o = 0.14 sqrt(dx^2 + dy^2)

Output: pressure p - p0 in bar, time in hours.
"""
from math import ceil, sqrt, log, exp, pi

import numpy as np

from utils.fs import ensure_path
from utils.vec import max_abs_diff
from writers.one_d.data import write_statistics
from writers.two_d.data import write_data


def e1(x):
    """Exponential integral E1(x), x > 0: series for x <= 1, continued fraction for x > 1."""
    if x <= 1.:
        s = 0.
        term = 1.
        kk = 1

        while True:
            term *= -x/kk
            add = -term/kk
            s += add

            if abs(add) < 1e-16*abs(s):
                break

            kk += 1

        return -0.5772156649015329 - log(x) + s

    # modified Lentz algorithm for E1(x) = exp(-x) / (x + 1 - 1/(x + 3 - 4/(x + 5 - ...)))
    tiny = 1e-300
    bb = x + 1.
    cc = 1./tiny
    dd = 1./bb
    hh = dd

    for i in range(1, 200):
        an = -i*i
        bb += 2.
        dd = 1./(an*dd + bb)
        cc = bb + an/cc
        de = cc*dd
        hh *= de

        if abs(de - 1.) < 1e-16:
            break

    return hh*exp(-x)


def theis(r, t):
    return c1*e1(r*r/(4*eta*t))


md = 9.869233e-16                  # 1 mD in m^2

k = 100*md                         # permeability, m^2
h = 10.                            # thickness, m
mu = 1e-3                          # viscosity, Pa s
phi = 0.2                          # porosity
ct = 1e-9                          # total compressibility, 1/Pa
q = 200./86400                     # injection rate, m^3/s (q < 0 for a producer)
rw = 0.1                           # well radius, m

dx = 20.
dy = 20.

lx = 2000.
ly = 2000.

bc = 'theis'                       # 'theis' | 'initial'

r = 0.2                            # eta dt / dx^2

end_time = 10*86400.               # s
out_time_step = 6*3600.            # s

out_dir = f'.out/2d/darcy_flow/well_model/transient/explicit_euler/{bc}/dx={dx}, dy={dy}, r={r}'

eta = k/(phi*mu*ct)

dt = min(r*dx*dx/eta, r*dy*dy/eta)

rx = eta*dt/dx/dx
ry = eta*dt/dy/dy

nx = ceil(lx / dx) + 1
ny = ceil(ly / dy) + 1

if nx % 2 == 0: nx += 1
if ny % 2 == 0: ny += 1

iw = nx//2
jw = ny//2

pw_node = jw + iw*ny

c0 = q*mu/(2*pi*k*h)               # steady radial factor
c1 = q*mu/(4*pi*k*h)               # Theis factor

# distances of the boundary nodes from the well
boundary = []

for i in range(nx):
    for j in range(ny):
        if i == 0 or j == 0 or i == nx - 1 or j == ny - 1:
            boundary.append((j + i*ny, sqrt(((i - iw)*dx)**2 + ((j - jw)*dy)**2)))


def set_boundary_condition(p, t):
    for pc, rb in boundary:
        p[pc] = theis(rb, t) if (bc == 'theis' and t > 0) else 0.


models = {
    'naive': rw,
    'equal_area': sqrt(dx*dy/pi),
    'peaceman': 0.14*sqrt(dx*dx + dy*dy),
}

p = np.zeros(nx*ny)                # p - p0, Pa
P = p.reshape(nx, ny)              # view: P[i, j] = p[j + i*ny]

pp = p.copy()

src = dt*q/(phi*ct*h*dx*dy)        # pressure increment of the well node per step

t = 0.
tn = out_time_step

n = 1
m = 0

stats = []
well = []

write_data(p/1e5, t/3600, dx, dy, m, nx, ny, out_dir)

m += 1

print(f'dt = {dt:.1f} s, rx = {rx:.3f}, ry = {ry:.3f}, steps = {ceil(end_time/dt)}')

while t < end_time:
    pn = P.copy()

    P[1:-1, 1:-1] = (pn[1:-1, 1:-1]
                     + rx*(pn[2:, 1:-1] - 2*pn[1:-1, 1:-1] + pn[:-2, 1:-1])
                     + ry*(pn[1:-1, 2:] - 2*pn[1:-1, 1:-1] + pn[1:-1, :-2]))

    P[iw, jw] += src

    t += dt

    set_boundary_condition(p, t)

    if t >= tn:
        write_data(p/1e5, t/3600, dx, dy, m, nx, ny, out_dir)

        max_diff = max_abs_diff(pp, p)/1e5

        pk = p[pw_node]
        pw_theis = theis(rw, t)

        pw = {name: pk + c0*log(ro/rw) for name, ro in models.items()}

        print(f't={t/3600:7.2f} h   Theis={pw_theis/1e5:7.3f} bar   '
              + '   '.join(f'{name}={pw[name]/1e5:7.3f} ({100*(pw[name] - pw_theis)/pw_theis:+6.2f}%)' for name in models))

        stats.append([n, t/3600, max_diff])
        well.append([t/3600, pw_theis/1e5] + [pw[name]/1e5 for name in models])

        pp = p.copy()

        m += 1

        tn += out_time_step

    n += 1

write_statistics(stats, out_dir)

ensure_path(f'{out_dir}/well_pressure.csv')
np.savetxt(f'{out_dir}/well_pressure.csv', np.array(well), delimiter=', ', fmt='%.5f',
           header='t_h, p_w_theis_bar, ' + ', '.join(f'p_w_{name}_bar' for name in models), comments='')
