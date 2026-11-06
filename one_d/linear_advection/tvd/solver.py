from math import ceil

import numpy as np

from writers.one_d.data import write_data

u = 1.

d = 0.2
l = 1.

dx = 0.01

a = 0.3

epsilon = 0.1

end_time = 0.3
out_time_step = 0.003

out_dir = f'.out/1d/linear_advection/tvd/u={u}, l={l}, d={d}/dx={dx}, a={a}, epsilon={epsilon}'

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

def q(x):
    if abs(x) < 2*epsilon:
        return x*x/4/epsilon + epsilon

    return abs(x)

while t <= end_time:
    for i in range(nx):
        Cl = C[1]
        Cr = C[nx-2]

        if i > 0:
            Cl = C[i-1]

        if i < nx-1:
            Cr = C[i+1]

        Cll = C1[2]
        Crr = C1[nx-3]

        if i == 1:
            Cll = C1[1]
        elif i > 1:
            Cll = C1[i-2]

        if i == nx-2:
            Crr = C1[nx-2]
        elif i < nx-2:
            Crr = C1[i+2]

        f = u*C[i]
        fr = u*Cr
        fl = u*Cl

        fp = (f + fr)/2
        fm = (fl + f)/2

        ar = a
        al = a

        arr = a
        all = a

        dr = Cr - C[i]
        dl = C[i] - Cl

        drr = Crr - Cr
        dll = Cl - Cll

        gr = (q(ar) - ar*ar)*dr/2
        gl = (q(al) - al*al)*dl/2

        grr = (q(arr) - arr*arr)*drr/2
        gll = (q(all) - all*all)*dll/2

        sr = 1 if gr > 0 else -1
        sl = 1 if gl > 0 else -1

        srr = 1 if grr > 0 else -1

        g = sr*max(0., min(abs(gr), gl*sr))
        gp = srr*max(0., min(abs(grr), gr*srr))
        gm = sl*max(0., min(abs(gl), gll*sl))

        br = (gp - g)/dr if dr > 0 else 0
        bl = (g - gm)/dl if dl > 0 else 0

        hp = (g + gp - q(ar + br)*dr)/2
        hm = (g + gm - q(al + bl)*dl)/2

        C1[i] = C[i] - (fp - fm)*dt/dx - hp + hm

    C, C1 = C1, C

    t += dt

    if t >= tn:
        write_data(C, t, dx, m, out_dir)

        print(f'Write data in file t={t:.3f}')

        m += 1

        tn += out_time_step

    n += 1
