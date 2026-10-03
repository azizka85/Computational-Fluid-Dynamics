import numpy as np

L = 1.

dx = 0.5

x = np.arange(0, L + dx, dx)

n = len(x)

def f(x):
    return 1 + np.cos(np.pi*x/2) + np.sin(np.pi*x)

def linear_interpolation(x, xm1, xc, xp1, ym1, yc, yp1):
    if x < xc:
        fm1 = (xc - x) / (xc - xm1)
        fc = (x - xm1) / (xc - xm1)

        return ym1 * fm1 + yc * fc

    fp1 = (x - xc) / (xp1 - xc)
    fc = (xp1 - x) / (xp1 - xc)

    return yc * fc + yp1 * fp1

def quadratic_interpolation(x, xm2, xm1, xc, xp1, xp2, ym2, ym1, yc, yp1, yp2):
    if x < xc:
        fm2 = (x - xm1)*(x - xc)/(xm2 - xm1)/(xm2 - xc)
        fm1 = (x - xm2)*(x - xc)/(xm1 - xm2)/(xm1 - xc)
        fc = (x - xm2)*(x - xm1)/(xc - xm2)/(xc - xm1)

        return ym2 * fm2 + ym1 * fm1 + yc * fc

    fc = (x - xp1)*(x - xp2)/(xc - xp1)/(xc - xp2)
    fp1 = (x - xc)*(x - xp2)/(xp1 - xc)/(xp1 - xp2)
    fp2 = (x - xc)*(x - xp1)/(xp2 - xc)/(xp2 - xp1)

    return yc * fc + yp1 * fp1 + yp2 * fp2

xt = np.zeros(0)

yl = np.zeros(0)
yq = np.zeros(0)

dxp = 0.01

for i in range(1, n-1):
    xmm = (x[i] + x[i-1]) / 2
    xpm = (x[i] + x[i+1]) / 2

    xp = np.arange(x[i-1], x[i+1] + dxp, dxp)

    xt = np.concatenate((xt, xp))
    yl = np.concatenate((
        yl,
        np.array([
            linear_interpolation(
                xc,
                x[i-1], x[i], x[i+1],
                f(x[i-1]), f(x[i]), f(x[i+1])
            )
            for xc in xp
        ])
    ))

    yq = np.concatenate((
        yq,
        np.array([
            quadratic_interpolation(
                xc,
                x[i-1], xmm, x[i], xpm, x[i+1],
                f(x[i-1]), f(xmm), f(x[i]), f(xpm), f(x[i+1])
            )
            for xc in xp
        ])
    ))

output.RowData.append(xt, 'x')
output.RowData.append(yl, 'yl')
output.RowData.append(yq, 'yq')
output.RowData.append(f(xt), 'ye')
