import numpy as np

x = inputs[0].Points[:,0]
y = inputs[0].Points[:,1]

w = np.sin(np.pi*x) *  np.sin(np.pi*y)

output.PointData.append(w, f'wa')
output.PointData.append(inputs[0].PointData['w'], f'w')
