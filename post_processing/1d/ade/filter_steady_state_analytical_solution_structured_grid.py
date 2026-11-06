import numpy as np

L = 1.

Pe = 100.

x = inputs[0].Points[:,0]

if Pe > 40:
    C = np.zeros_like(x)
    C[x == L] = 1
else:
    C = (exp(Pe*x/L) - 1)/(exp(Pe) - 1)

output.PointData.append(C, f'C_anal')
output.PointData.append(inputs[0].PointData['C'], f'C')
