from math import ceil

import numpy as np

from slae.direct.tridiagonal import solve_uniform

from writers.one_d.data import write_data

u_top = 1.

h = 1.
dy = 0.01

out_dir = f'.out/1d/couette_flow/steady_state/U={u_top}, H={h}/dy={dy}'

ny = ceil(h / dy) + 1

u = np.zeros(ny)
u1 = np.zeros(ny)

u[0] = 0
u[ny-1] = u_top

al = 1
al1 = 0

ac = -2
ac0 = 1
ac1 = 1

ar = 1
ar0 = 0

solve_uniform(al, al1, ac, ac0, ac1, ar, ar0, u, u1)

write_data(u1, 0, dy, 0, out_dir)
