import numpy as np

U = 1.
d = 0.2

x = inputs[0].Points[:,0]
t = inputs[0].FieldData['Time']

C = np.zeros_like(x)

C[x <= d + U*t] = 1.

output.PointData.append(C, f'Ca')
output.PointData.append(inputs[0].PointData['C'], f'C')
