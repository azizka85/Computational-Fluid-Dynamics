import numpy as np

nu = 1.

a = 2.
b = 1.

x = inputs[0].Points[:,0]
y = inputs[0].Points[:,1]

t = inputs[0].FieldData['Time']

w = inputs[0].PointData['w']
wa = np.cos(np.pi*x/4) *  np.cos(np.pi*y/2) * np.exp(-5*np.pi**2*nu*t/16)

wa[x**2/a**2 + y**2/b**2 > 1] = 0

output.PointData.append(wa, f'wa')
output.PointData.append(w, f'w')
