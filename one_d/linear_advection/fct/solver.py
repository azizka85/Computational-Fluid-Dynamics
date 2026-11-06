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

out_dir = f'.out/1d/linear_advection/fct/u={u}, l={l}, d={d}/dx={dx}, a={a}'

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

        C1[i] = C[i] - a*(C[i] - Cl)

    for i in range(nx):
        Crr = C1[nx-3]

        if i == nx-2:
            Crr = C1[nx-2]
        elif i < nx-2:
            Crr = C1[i+2]

        Cr = C1[nx-2]

        if i < nx-1:
            Cr = C1[i+1]

        dp1 = Cr - C1[i]
        dp2 = Crr - Cr

        Cll = C1[2]

        if i == 1:
            Cll = C1[1]
        elif i > 1:
            Cll = C1[i-2]

        Cl = C1[1]

        if i > 0:
            Cl = C1[i-1]

        dm1 = C1[i] - Cl
        dm2 = Cl - Cll

        ap = a
        am = a

        nup = ap*(1 - ap)/2
        num = am*(1 - am)/2

        s = -1 if dp1 < 0 else 1

        gp = s*max(0., min(s*dm1, min(nup*abs(dp1), s*dp2)))

        s = -1 if dm1 < 0 else 1

        gm = s*max(0., min(s*dm2, min(num*abs(dm1), s*dp1)))

        C[i] = C1[i] + (gm - gp)

    t += dt

    if t >= tn:
        write_data(C, t, dx, m, out_dir)

        print(f'Write data in file t={t:.3f}')

        m += 1

        tn += out_time_step

    n += 1
