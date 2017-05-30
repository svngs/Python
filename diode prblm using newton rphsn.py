import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

N = 10
v = np.zeros(N+1)
I = np.zeros(N)
f = np.zeros(N)
fd = np.zeros(N)
err = np.zeros(N+1)
v0 = 2
I0 = 1e-9
vt = 25e-3
R = 1e3
v[0] = 0.35

for j in range(N) :
	f[j]   =v0-R*(I0*(np.exp(v[j]/vt)-1))- v[j]
	fd[j]  =-R*(I0/vt)*(np.exp(v[j]/vt))- 1
	v[j+1] = v[j]-(f[j]/fd[j])

if((np.abs(v[j+1]-v[j])/np.abs(v[j]))<1e-2):
	
	print('V = ',v[j+1])
	I = I0*(np.exp(v[j]/vt)- 1)
	print('I = ',I)

for i in range (N):
	err[i] = np.abs(v[i+1]/v[i] - 1)
plt.plot(range(N+1), err)
plt.title('Newton Raphson Method')
plt.xlabel('Values of N')
plt.ylabel('Error')
plt.show()
