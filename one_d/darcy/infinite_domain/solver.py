"""
1D single-phase Darcy flow: injecting (producing) well in an infinite domain.

Dimensionless pressure equation (Lecture 7):
    dp/dt = d2p/dx2 + delta(x),   p(x, 0) = 0,   p -> 0 as |x| -> inf
Scales:  x = L*x~,  t = L^2/eta * t~,  p - p0 = q*mu*L/k * p~,  eta = k/(phi*mu*c_t)

Analytical solution:
    p(x, t) = sqrt(t/pi)*exp(-x^2/(4t)) - |x|/2 * erfc(|x|/(2*sqrt(t)))

Vertex-centred finite volumes + implicit Euler:
    V_j (p_j^{n+1} - p_j^n)/dt = T_{j+1/2}(p_{j+1} - p_j) - T_{j-1/2}(p_j - p_{j-1}) + Q_j

bc = 'dirichlet' : uniform grid on [-L, L],  p = 0 at x = +-L
bc = 'neumann'   : uniform grid on [-L, L],  dp/dn = 0 at x = +-L
bc = 'mapped'    : x = Lc*tan(pi*xi/2), xi in [-1, 1],  p = 0 at xi = +-1 (x = +-inf)
                   V_j = J(xi_j)*dxi,  T_{j+1/2} = 1/(J(xi_{j+1/2})*dxi),  J = dx/dxi

Output: write_data writes node j at position j*h, so the files show
    x + L        for 'dirichlet' and 'neumann'  (h = dx)
    xi + 1       for 'mapped'                   (h = dxi, computational coordinate)
"""
from math import sqrt, exp, erfc, pi, tan, cos

import numpy as np

from utils.fs import ensure_path
from writers.one_d.data import write_data, write_statistics


def solve_tridiagonal(a, b, c, d):
    """Thomas algorithm: a[j] p[j-1] + b[j] p[j] + c[j] p[j+1] = d[j]."""
    n = len(d)

    cp = np.zeros(n)
    dp = np.zeros(n)

    cp[0] = c[0]/b[0]
    dp[0] = d[0]/b[0]

    for j in range(1, n):
        den = b[j] - a[j]*cp[j-1]

        cp[j] = c[j]/den
        dp[j] = (d[j] - a[j]*dp[j-1])/den

    p = np.zeros(n)
    p[n-1] = dp[n-1]

    for j in range(n-2, -1, -1):
        p[j] = dp[j] - cp[j]*p[j+1]

    return p


def p_exact(x, t):
    x = abs(x)

    return sqrt(t/pi)*exp(-x*x/(4*t)) - x/2*erfc(x/(2*sqrt(t)))


bc = 'mapped'          # 'dirichlet' | 'neumann' | 'mapped'

l = 1.                 # half-length of the truncated domain [-L, L]
lc = 1.                # mapping scale, x = Lc*tan(pi*xi/2)

n = 100                # number of intervals (even: the well is node n/2)

dt = 0.001
end_time = 100.
out_time_step = 1.

out_dir = f'.out/1d/darcy/well/{bc}/n={n}, L={l}, Lc={lc}, dt={dt}'

# ----------------------------------------------------------------- grid
if bc == 'mapped':
    h = 2./n                                            # dxi

    xi = np.array([-1. + j*h for j in range(n + 1)])

    x = np.array([lc*tan(pi*s/2) for s in xi])
    x[0], x[n] = -np.inf, np.inf

    jac = lambda s: pi*lc/2/cos(pi*s/2)**2              # J = dx/dxi

    v = np.array([jac(s)*h for s in xi])                # end values are not used
    tr = np.array([1./(jac(xi[j] + h/2)*h) for j in range(n)])
else:
    h = 2.*l/n                                          # dx

    x = np.array([-l + j*h for j in range(n + 1)])

    v = np.full(n + 1, h)
    v[0] = v[n] = h/2                                   # half control volumes

    tr = np.full(n, 1./h)

jw = n//2                                               # well node, x = 0

q = np.zeros(n + 1)
q[jw] = 1.

# ----------------------------------------------------------------- matrix
# row j: a p_{j-1} + b p_j + c p_{j+1} = V_j/dt p_j^n + Q_j
a = np.zeros(n + 1)
b = np.zeros(n + 1)
c = np.zeros(n + 1)

for j in range(n + 1):
    tl = tr[j-1] if j > 0 else 0.                       # T_{j-1/2}; 0 outside -> no flux
    tg = tr[j] if j < n else 0.                         # T_{j+1/2}

    a[j] = -tl
    c[j] = -tg
    b[j] = v[j]/dt + tl + tg

if bc in ('dirichlet', 'mapped'):                       # p = 0 at the end nodes
    a[0] = c[0] = 0.
    b[0] = 1.

    a[n] = c[n] = 0.
    b[n] = 1.

# ----------------------------------------------------------------- time loop
p = np.zeros(n + 1)
d = np.zeros(n + 1)

t = 0.
tn = out_time_step

it = 1
m = 0

stats = []
well = []

write_data(p, t, h, m, out_dir)

m += 1

while t < end_time - 1e-12:
    for j in range(n + 1):
        d[j] = v[j]/dt*p[j] + q[j]

    if bc in ('dirichlet', 'mapped'):
        d[0] = d[n] = 0.

    p = solve_tridiagonal(a, b, c, d)

    t += dt

    if t >= tn - 1e-9:
        write_data(p, t, h, m, out_dir)

        # max |p - p_exact| over the finite nodes
        max_diff = max(abs(p[j] - p_exact(x[j], t)) for j in range(n + 1) if np.isfinite(x[j]))

        pw_exact = p_exact(0., t)

        print(f'Write data in file t={t:.3f}, p_w={p[jw]:.5f}, exact={pw_exact:.5f}, '
              f'error={100*(p[jw] - pw_exact)/pw_exact:+.2f}%')

        stats.append([it, t, max_diff])
        well.append([t, p[jw], pw_exact])

        m += 1

        tn += out_time_step

    it += 1

write_statistics(stats, out_dir)

ensure_path(f'{out_dir}/well_pressure.csv')
np.savetxt(f'{out_dir}/well_pressure.csv', np.array(well), delimiter=', ', fmt='%.6f',
           header='t, p_w, p_w_exact', comments='')
