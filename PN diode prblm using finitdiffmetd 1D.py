import numpy as np
import matplotlib.pyplot as plt
import PhysUnit as Unit

xmin = -3*Unit.um
xmax = 3*Unit.um

N = 600
A = np.zeros([N+1,N+1])
B = np.zeros(N+1)
Phi = np.zeros(N+1)

x = np.linspace( xmin, xmax, N+1)
Delta = x[1] - x[0]

def  do_assemble( ):
    global A, B, Phi, N, x, Delta
    for j in range(N+1) :
        A[j,j] = -2.0 / Delta**2
        if ( j == 0 ) :
            A[0, 1] = 1.0 / Delta**2
        elif ( j == N ) :
            A[N, N-1] = 1.0 / Delta**2
        else :
            A[j, j-1 ] = 1.0 / Delta**2
            A[j, j+1 ] = 1.0 / Delta**2
            
    for j in range(N+1) :
        
        if (j > N/3 and j< N/2) :
            B[j] = 1e15*((1.6e-19*Unit.C*pow(Unit.cm,-3))/(12*Unit.eps0))

        elif ( j > N/2 and j < 2*N/3) :
            B[j] = -1e15*((1.6e-19*Unit.C*pow(Unit.cm,-3))/(12*Unit.eps0))

def do_Dirichlet( ):
    global A, B, Phi, N, x, Delta
    
    eps = 7./3 - 4./3 -1
    L = 1.0/eps
    
    A[0,0] = L / Delta**2
    B[0] = 0.0 * L / Delta**2
    
def do_Neumann( ):
    global A, B, Phi, N, x, Delta 
    
    eps = 7./3 - 4./3 -1
    L = 1.0/eps
    
    A[N,N] = L / Delta**2
    A[N, N-1] =  -L / Delta**2   
    B[N] = 0.0 * L / Delta
    
def do_solve( ):
    global A, B, Phi, N,  x, Delta
    Phi = np.linalg.solve( A, B )

def do_plot( ):
    global A, B, Phi, N,  x, Delta
    
    plt.plot( x, Phi )
    plt.title('P-N Diode Characteristics')
    plt.xlim(xmax=max(x))
    plt.show()

if __name__ == '__main__':
    do_assemble() ;
    do_Dirichlet( ) ;
    do_Neumann( ) ;
    do_solve( );
    do_plot();    
