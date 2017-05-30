import numpy as np
import matplotlib.pyplot as plt
import PhysUnit as Unit

xmin = 0
xmax = 1

ymin = 0
ymax = 1

N = 40

A = np.zeros([N*N,N*N])
B = np.zeros([N*N])
V = np.zeros([N*N])
V0 = np.zeros([N*N])

x = np.linspace( xmin, xmax, N)
y = np.linspace( ymin, ymax, N)

Delta = x[1] - x[0]

def  do_assemble( ):
    global A, B, V, N, x, Delta
    for j in range(N*N) :
        A[j,j] = -4.0/ Delta**2
        if ( j == 0 ) :    
            A[0, 1] = 1.0/ Delta**2
            A[j, j+N ] = 1.0/ Delta**2
        elif ( j == N*N-1 ) :
            A[N*N-1, N*N-2] = 1.0/ Delta**2
            A[j, j-N ] = 1.0/ Delta**2
        elif (j < N)   :
            A[j, j-1 ] = 1.0/ Delta**2
            A[j, j+1 ] = 1.0/ Delta**2
            A[j, j+N] = 1.0/ Delta**2
        elif (j > (N*N-1-N))   :
            A[j, j+1 ] = 1.0/ Delta**2
            A[j, j-1 ] = 1.0/ Delta**2
            A[j, j-N] = 1.0/ Delta**2
        else :
            A[j, j+1 ] = 1.0/ Delta**2
            A[j, j-1 ] = 1.0/ Delta**2
            A[j, j-N ] = 1.0/ Delta**2
            A[j, j+N ] = 1.0/ Delta**2

            B[j] = 0.0
            
def do_Dirichlet( ):
    global A, B, V, N, x, Delta
    eps = 7./3 - 4./3 -1
    L = 1.0/eps
    for k in range(N)   :
        A[k, k] = L  / Delta**2
        B[k]= 0.0*L /Delta**2
    for k in range(N-1, N*N, N) :
        A[k, k] = L   / Delta**2
        B[k]    = 0.0* L / Delta**2
    for k in range(N*N-1, N*N)  :
        A[k, k] = L   / Delta**2
        B[k]    = 0.0* L / Delta**2
    for k in range (0, N*N-N+1,N)   :
        A[k,k] = L/Delta**2
        B[k]    =1.0*L/Delta**2

def do_solve( ):
    global A, B, V, N,  x, Delta,V0
    V = np.linalg.solve( A, B )
    V0= V.reshape(N,N)

def do_plot( ):
    global A, B, V, N,  x, Delta,V0
    plt.contourf( x,y, V0 )
    plt.xlim(xmax=max(x))
    plt.show()

    
if __name__ == '__main__':
    do_assemble() ;
    do_Dirichlet( ) ;
    do_solve( );
    do_plot();  
