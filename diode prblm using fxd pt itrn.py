import scipy as sp
import numpy as np
import matplotlib.pyplot as plt

N = 10
V0 = 2
I0 = 1e-9
Vt = 25e-3
R = 1e3
I = np.zeros(N+1)
Vd = np.zeros(N+1)
err = np.zeros(N+1)

I[0] = 2e-3
Vd[0] = -2
for i in range (N):
	I[i+1] = V0/R - (Vt/R)*(np.log(1+I[i]/I0))
	Vd[i+1]= V0 - I[i+1]*R
	
for i in range (N):
	err[i] = np.abs(Vd[i+1]/Vd[i] - 1)
plt.plot(range(N+1), err)
plt.title('Fixed Point Iteration')
plt.xlabel('Values of N')
plt.ylabel('Error')
plt.show()
	
