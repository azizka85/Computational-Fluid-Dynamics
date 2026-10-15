import numpy as np

a = 1.

x = inputs[0].Points[:,0]
y = inputs[0].Points[:,1]

w = np.sin(np.pi*x) *  np.cos(np.pi*y)
w[(y > 0) | (y < -a*np.sqrt(3)/2 + np.sqrt(3)*x) | (y < -a*np.sqrt(3)/2 - np.sqrt(3)*x)] = 0

output.PointData.append(w, f'wa')
output.PointData.append(inputs[0].PointData['w'], f'w')
