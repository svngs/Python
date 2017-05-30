import numpy as np
import scipy as sp
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

N = 10
v = np.zeros(N)
v0 = 2.0
vt = 25e-3
I0 = 1e-9
R = 1e3

def f(i) :
    return i - I0*(np.exp((v0 - i*R)/vt)-1)
    
Id = fsolve(f, 1.6e-3)
Vd = v0 - Id[0]*R
print(Id[0])
print(Vd)
