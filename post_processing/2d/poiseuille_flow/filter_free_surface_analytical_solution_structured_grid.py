import numpy as np

nu = 1.

x = inputs[0].Points[:,0]
y = inputs[0].Points[:,1]

t = inputs[0].FieldData['Time']

w = np.sin(np.pi*x) *  np.sin(np.pi*y/2) * np.exp(-5*np.pi**2*nu*t/4)

output.PointData.append(w, f'wa')
output.PointData.append(inputs[0].PointData['w'], f'w')
