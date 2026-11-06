from math import ceil

import numpy as np

from writers.one_d.data import write_data

l = 1.

dx = 0.01

a = 0.3
r = 0.3

end_time = 0.3
out_time_step = 0.003

out_dir = f'.out/1d/ade/eno_x2t2/l={l}/dx={dx}, a={a}, r={r}'

dt = dx

nx = ceil(l / dx) + 1

C = np.zeros(nx)
C1 = np.zeros(nx)

C[nx-1] = 1.

t = 0.
tn = out_time_step

n = 1
m = 0

write_data(C, t, dx, m, out_dir)

m += 1

while t <= end_time:
    for i in range(1, nx-1):
        vl = (C[i] - C[i-1])/2
        vr = (C[i+1] - C[i])/2

        v = vl if abs(vl) < abs(vr) else vr

        Cp = C[i] + v
        Cll = C[0]

        if i > 1:
            Cll = C[i-2]

        vl = (C[i-1] - Cll)/2
        vr = (C[i] - C[i-1])/2

        v = vl if abs(vl) < abs(vr) else vr

        Cm = C[i-1] + v

        C1[i] = C[i] - a*(Cp - Cm) + r*(C[i+1] - 2*C[i] + C[i-1])

    C1[0] = 0.
    C1[nx-1] = 1.

    for i in range(1, nx-1):
        vl = (C1[i] - C1[i-1])/2
        vr = (C1[i+1] - C1[i])/2

        v = vl if abs(vl) < abs(vr) else vr

        Cp = C1[i] + v
        Cll = C1[0]

        if i > 1:
            Cll = C1[i-2]

        vl = (C1[i-1] - Cll)/2
        vr = (C1[i] - C1[i-1])/2

        v = vl if abs(vl) < abs(vr) else vr

        Cm = C1[i-1] + v

        C[i] = 0.5*C[i] + 0.5*C1[i] - a*(Cp - Cm)/2 + r*(C1[i+1] - 2*C1[i] + C1[i-1])/2

    C[0] = 0.
    C[nx-1] = 1.

    t += dt

    if t >= tn:
        write_data(C, t, dx, m, out_dir)

        print(f'Write data in file t={t:.3f}')

        m += 1

        tn += out_time_step

    n += 1
