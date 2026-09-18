import numpy as np

t = inputs[0].FieldData['t'][0]
H = inputs[0].FieldData['H'][0]
U = inputs[0].FieldData['U'][0]
nu = inputs[0].FieldData['nu'][0]

dy = 0.1

y = np.arange(0, H + dy, dy) 

output.RowData.append(y, 'y')

def get_solution(m):
	u = U*y/H

	for n in range(1, m + 1):
		u += 2*U*exp(-n**2*np.pi**2*nu*t/H**2)*cos(np.pi*n)*sin(np.pi*n*y/H)/np.pi/n
		
	return u

terms = [1, 3, 5, 10, 100]

for term in terms:
	u = get_solution(term)
	
	output.RowData.append(u, f'u{term:03}')
