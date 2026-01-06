# quantum-harmonic-oscillator-
Numerical solution of the 1D QHO using finite differences approach


This project implements a finite-difference numerical solution of the
time-independent Schrödinger equation for the one-dimensional quantum
harmonic oscillator. Numerical eigenvalues and eigenfunctions are
computed and validated against analytic solutions.

We consider a particle of mass m in a harmonic potential
V(x) = ½ m ω² x² and solve the time-independent Schrödinger equation

Ĥψ = Eψ

where Ĥ is the Hamiltonian operator. 

The spatial domain is discretized using a uniform grid and the second
derivative is approximated using a central finite-difference approach.
This results in a tridiagonal Hamiltonian matrix, which is diagonalized
numerically to obtain energy eigenvalues and eigenfunctions.

- Grid size (N = 300): chosen to balance accuracy and computational cost
- Domain size (x ∈ [−8, 8]): large enough to capture wavefunction decay but remains within reasonable range of computational cost
- Units: ħ = m = ω = 1 (natural units)

## Repository Structure
- `qho.py` – main solver
- `plots/` – eigenfunction and probability density plots
- `report.pdf` – short project report 
