import numpy as np

nu = 1.

a = 1.

x = inputs[0].Points[:,0]
y = inputs[0].Points[:,1]

t = inputs[0].FieldData['Time']

w = inputs[0].PointData['w']
wa = np.cos(np.pi*x) *  np.cos(np.pi*y) * np.exp(-2*np.pi**2*nu*t)

wa[(y > 0) | (y < -a*np.sqrt(3)/2 + np.sqrt(3)*x) | (y < -a*np.sqrt(3)/2 - np.sqrt(3)*x)] = 0

output.PointData.append(wa, f'wa')
output.PointData.append(w, f'w')
