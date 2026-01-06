import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh


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


nstates_plot = 5    
nstates_table = 10   

for n in range(min(states.shape[1], nstates_table)):
    psi = states[:, n]
    norm = np.sqrt(np.sum(np.abs(psi)**2) * dx)
    states[:, n] = psi / norm


hbar = 1.0

print("\nExpectation values + uncertainties (ħ=1):")
print("n   <x>       Δx       <p>       Δp      ΔxΔp")
print("------------------------------------------------")

for n in range(min(states.shape[1], nstates_table)):
    psi = states[:, n]

    # <x>, <x^2>
    x_mean  = np.sum(np.conjugate(psi) * (x * psi)) * dx
    x2_mean = np.sum(np.conjugate(psi) * ((x**2) * psi)) * dx

    # p = -i ħ d/dx (momentum operator)
    dpsi = np.gradient(psi, dx)
    ppsi = -1j * hbar * dpsi
    p_mean = np.sum(np.conjugate(psi) * ppsi) * dx

    # p^2 = -ħ^2 d^2/dx^2
    d2psi = np.gradient(np.gradient(psi, dx), dx)
    p2psi = -(hbar**2) * d2psi
    p2_mean = np.sum(np.conjugate(psi) * p2psi) * dx

    # Variances
    var_x = np.real(x2_mean - x_mean**2)
    var_p = np.real(p2_mean - p_mean**2)

    Dx = np.sqrt(max(var_x, 0.0))
    Dp = np.sqrt(max(var_p, 0.0))
    hup = Dx * Dp

    print(f"{n:<2d}  {np.real(x_mean):>8.4f}  {Dx:>8.4f}  {np.real(p_mean):>8.4f}  {Dp:>8.4f}  {hup:>8.4f}")


#Potential well - V(x) = 1/2mω^2x^2 - parabolic shape
plt.figure()
plt.plot(x, V)
plt.xlabel("x")
plt.ylabel("V(x)")
plt.title("Harmonic oscillator potential")
plt.show()


#Eigenfunctions offset by energies

plt.figure()
plt.plot(x, V, label="V(x)")

scale = 0.6
for n in range(min(nstates_plot, states.shape[1])):
    plt.plot(x, scale * states[:, n] + eigvals[n], label=f"ψ{n} + E{n}")
    plt.hlines(eigvals[n], x[0], x[-1], linestyles="dotted", linewidth=1)

plt.xlabel("x")
plt.ylabel("Energy / offset wavefunctions")
plt.title("First eigenfunctions (offset by energies) + potential")
plt.legend(fontsize=8)
plt.show()


#Probability densities

plt.figure()
for n in range(min(nstates_plot, states.shape[1])):
    plt.plot(x, np.abs(states[:, n])**2, label=f"|ψ{n}|²")
plt.xlabel("x")
plt.ylabel("Probability density")
plt.title("Probability densities (n=0..4)")
plt.legend()
plt.show()

# Energy spectrum comparison: numerical vs analytical

n_compare = min(10, len(eigvals))
E_ana = np.array([n + 0.5 for n in range(n_compare)])

plt.figure()
plt.scatter(range(n_compare), eigvals[:n_compare], label="Numerical")
plt.plot(range(n_compare), E_ana, "--", label="Analytic (n+1/2)")
plt.xlabel("n")
plt.ylabel("Energy")
plt.title("Energy spectrum: numerical vs analytic")
plt.legend()
plt.show()



