import numpy as np

U = 1.

H = 1.

y = inputs[0].Points[:,1]
u = U*y/H

output.PointData.append(u, f'u_anal')
output.PointData.append(inputs[0].PointData['u'], f'u')
