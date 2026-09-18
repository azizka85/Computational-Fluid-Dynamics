import numpy as np

nu = 1.
U = 1.

H = 1.

dy = 0.01

t = 0.

output.FieldData.append(t, 't')
output.FieldData.append(H, 'H')
output.FieldData.append(U, 'U')
output.FieldData.append(nu, 'nu')

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
