import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh



#Part 1

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



#Part 2 -  Eigenvalues problem


eigvals, eigvecs = np.linalg.eigh(H)

assert np.all(np.diff(eigvals) >= 0)

n_states = 10
energies = eigvals[:n_states]
states = eigvecs[:, :n_states]   # shape: (N, 10)

for n in range(n_states):
    psi = states[:, n]
    norm = np.sqrt(np.sum(np.abs(psi)**2) * dx)
    states[:, n] = psi / norm



print("First 10 energy eigenvalues:")
for n in range(10):
    print(f"n = {n:2d}   E_num = {eigvals[n]:.6f}")


print("\nComparison with analytic energies:")
print(" n    E_num      E_ana      rel_error")
print("---------------------------------------")

for n in range(10):
    E_num = eigvals[n]
    E_ana = n + 0.5
    rel_err = abs(E_num - E_ana) / E_ana
    print(f"{n:2d}  {E_num:9.6f}  {E_ana:9.6f}  {rel_err:9.2e}")

print("\nGrid parameters:")
print(f"L = {L}")
print(f"N = {N}")
print(f"dx = {dx:.6f}")


psi0 = states[:, 0]
psi1 = states[:, 1]

even_check = np.allclose(psi0, psi0[::-1], atol=1e-2)
odd_check  = np.allclose(psi1, -psi1[::-1], atol=1e-2)

print("\nParity checks:")
print("Ground state even:", even_check)
print("First excited odd:", odd_check)



