import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh



Part 1

#Building the numerical Hamiltonian: H = T (kinetic energy) +  V (potential energy)


#parameters
L = 8.0
N = 300

x = np.linspace(-L, L, N)
dx = x[1] - x[0]   #mathematicallly equivalent to 2L/N - 1


diag = -2.0 * np.ones(N)
off_diag = 1.0 * np.ones(N - 1)

D2 = (
    np.diag(diag)
    + np.diag(off_diag, k=1)
    + np.diag(off_diag, k=-1)
) / dx**2


assert D2.shape == (N, N) 

#check
f = x**2
d2f_numeric = D2 @ f

#kinetic energy operator T

T = -0.5 * D2

#potential energy operator V(x) -  1/2*x^2

V = 0.5 * x**2
V_mat = np.diag(V)


#total energy hamiltonian
H = T + V_mat

assert np.allclose(H, H.T)







