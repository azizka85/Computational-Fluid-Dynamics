import numpy as np

U = 1.
H = 1.

nu = 1.

y = inputs[0].Points[:,1]
t = inputs[0].FieldData['Time']

u = U

for n in range(0,11):
    u -= 4*U*exp(-((2*n + 1.)/2)**2*np.pi**2*nu*t/H**2)*sin(np.pi*(2*n + 1.)*y/2/H)/np.pi/(2*n + 1)

output.PointData.append(u, f'u10')
output.PointData.append(inputs[0].PointData['u'], f'u')
