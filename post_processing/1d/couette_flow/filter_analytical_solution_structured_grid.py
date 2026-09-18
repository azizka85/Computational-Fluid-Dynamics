import numpy as np

U = 1.
H = 1.

nu = 1.

y = inputs[0].Points[:,1]
t = inputs[0].FieldData['Time']

u = U*y/H

for n in range(1,11):
    u += 2*U*exp(-n**2*np.pi**2*nu*t/H**2)*cos(np.pi*n)*sin(np.pi*n*y/H)/np.pi/n

output.PointData.append(u, f'u10')
output.PointData.append(inputs[0].PointData['u'], f'u')
