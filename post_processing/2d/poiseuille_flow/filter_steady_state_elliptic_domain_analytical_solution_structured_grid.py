import numpy as np

g = 1.
nu = 1.

a = 2.
b = 1.

x = inputs[0].Points[:,0]
y = inputs[0].Points[:,1]

A = -g*a**2*b**2/(2*nu*(a**2 + b**2))

w = A*(1 - (x/a)**2 - (y/b)**2)
w[x**2/a**2 + y**2/b**2 > 1] = 0

output.PointData.append(w, f'wa')
output.PointData.append(inputs[0].PointData['w'], f'w')
