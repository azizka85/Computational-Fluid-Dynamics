from math import ceil

import numpy as np

from slae.direct.tridiagonal import solve_uniform

from writers.one_d.data import write_data

u = 1.

d = 0.2
l = 1.

dx = 0.01

a = 0.3

end_time = 0.3
out_time_step = 0.003

out_dir = f'.out/1d/linear_advection/implicit_euler/u={u}, l={l}, d={d}/dx={dx}, a={a}'

dt = a*dx/u

nx = ceil(l / dx) + 1

C = np.zeros(nx)
C1 = np.zeros(nx)

cutoff = int(d // dx) + 1

C[:cutoff] = 1.

t = 0.
tn = out_time_step

n = 1
m = 0

write_data(C, t, dx, m, out_dir)

m += 1

al = -0.5*a
al1 = -0.5*a

ac = 1
ac0 = 1
ac1 = 1 + 0.5*a

ar = 0.5*a
ar0 = 0

while t <= end_time:
    solve_uniform(al, al1, ac, ac0, ac1, ar, ar0, C, C1)

    C1[0] = 1.

    C, C1 = C1, C

    t += dt

    if t >= tn:
        write_data(C, t, dx, m, out_dir)

        print(f'Write data in file t={t:.3f}')

        m += 1

        tn += out_time_step

    n += 1
