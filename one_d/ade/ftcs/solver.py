from math import ceil

import numpy as np

from writers.one_d.data import write_data

l = 1.

dx = 0.01

a = 0.3
r = 0.3

end_time = 0.3
out_time_step = 0.003

out_dir = f'.out/1d/ade/ftcs/l={l}/dx={dx}, a={a}, r={r}'

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
        C1[i] = C[i] - a*(C[i+1] - C[i-1])/2 + r*(C[i+1] - 2*C[i] + C[i-1])

    C1[0] = 0.
    C1[nx-1] = 1.

    C, C1 = C1, C

    t += dt

    if t >= tn:
        write_data(C, t, dx, m, out_dir)

        print(f'Write data in file t={t:.3f}')

        m += 1

        tn += out_time_step

    n += 1
