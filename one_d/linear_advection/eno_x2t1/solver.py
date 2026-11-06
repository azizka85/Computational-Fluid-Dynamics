from math import ceil

import numpy as np

from writers.one_d.data import write_data

u = 1.

d = 0.2
l = 1.

dx = 0.01

a = 0.3

end_time = 0.3
out_time_step = 0.003

out_dir = f'.out/1d/linear_advection/eno_x2t1/u={u}, l={l}, d={d}/dx={dx}, a={a}'

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

while t <= end_time:
    for i in range(nx):
        Cl = C[1]
        Cr = C[nx-2]

        if i > 0:
            Cl = C[i-1]

        if i < nx-1:
            Cr = C[i+1]

        vl = (C[i] - Cl)/2
        vr = (Cr - C[i])/2

        v = vl if abs(vl) < abs(vr) else vr

        Cp = C[i] + v
        Cll = C[2]

        if i == 1:
            Cll = C[1]
        elif i > 1:
            Cll = C[i-2]

        vl = (Cl - Cll)/2
        vr = (C[i] - Cl)/2

        v = vl if abs(vl) < abs(vr) else vr

        Cm = Cl + v

        C1[i] = C[i] - a * (Cp - Cm)

    C, C1 = C1, C

    t += dt

    if t >= tn:
        write_data(C, t, dx, m, out_dir)

        print(f'Write data in file t={t:.3f}')

        m += 1

        tn += out_time_step

    n += 1
