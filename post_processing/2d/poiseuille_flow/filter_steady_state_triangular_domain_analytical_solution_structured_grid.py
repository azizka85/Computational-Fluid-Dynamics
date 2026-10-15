import numpy as np

g = 1.
nu = 1.

a = 1.

x = inputs[0].Points[:,0]
y = inputs[0].Points[:,1]

A = g/2/sqrt(3)/nu/a

w = A*y*(y + a*sqrt(3)/2 - sqrt(3)*x)*(y + a*sqrt(3)/2 + sqrt(3)*x)
w[(y > 0) | (y < -a*np.sqrt(3)/2 + np.sqrt(3)*x) | (y < -a*np.sqrt(3)/2 - np.sqrt(3)*x)] = 0

output.PointData.append(w, f'wa')
output.PointData.append(inputs[0].PointData['w'], f'w')
