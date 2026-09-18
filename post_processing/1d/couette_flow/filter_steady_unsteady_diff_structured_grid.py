import numpy as np

u1 = inputs[0].PointData['u']
u2 = inputs[1].PointData['u']

diff = np.abs(u1 - u2)

output.PointData.append(diff, f'diff')
