import numpy as np

r = np.array([0.1, 0.5, 1, 10, 100])

output.RowData.append(r, 'r')

b = [pi/180, pi/30, pi/10, pi/3, pi]

for i in range(len(b)):
    g_bi = np.abs(1/(1 + 2*r*(1 - cos(b[i]))))
    g_cn = np.abs((1 - r*(1 - cos(b[i])))/(1 + r*(1 - cos(b[i]))))

    output.RowData.append(g_bi, f'bi{i}')
    output.RowData.append(g_cn, f'cn{i}')
