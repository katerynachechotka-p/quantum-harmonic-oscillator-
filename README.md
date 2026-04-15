# Quantum Harmonic Oscillator — Numerical Solution

Numerical solution of the 1D time-independent Schrödinger equation using 
finite differences, with an extension to the anharmonic case to investigate 
where perturbative approximations break down.

## Method

Discretizes the spatial domain into N=300 points on x ∈ [−8, 8], approximates 
the second derivative with a central finite-difference stencil, and assembles 
a tridiagonal Hamiltonian matrix. Diagonalization via `numpy.linalg.eigh` 
yields energy eigenvalues and eigenfunctions. Units: ħ = m = ω = 1.

## Results — Harmonic Oscillator

Energy eigenvalues validated against analytic formula Eₙ = n + 0.5:

| n | E_numerical | E_analytic | Relative error |
|---|---|---|---|
| 0 | 0.500000 | 0.5 | ~10⁻⁶ |
| 1 | 1.500001 | 1.5 | ~10⁻⁶ |
| 2 | 2.500003 | 2.5 | ~10⁻⁶ |

Relative error remains below 2% for all 10 lowest states. Ground state 
saturates the Heisenberg uncertainty bound: ΔxΔp ≈ 0.5. Parity checks 
confirm ground state is even, first excited state is odd.

## Results — Anharmonic Oscillator V(x) = ½x² + λx⁴

First-order perturbation theory prediction: E₀ ≈ E₀⁽⁰⁾ + ⟨ψ₀|λx⁴|ψ₀⟩

| λ | E_numerical | E_perturbation | Relative error |
|---|---|---|---|
| 0.001 | 0.500657 | 0.500660 | 5.20×10⁻⁶ |
| 0.010 | 0.507162 | 0.507406 | 4.80×10⁻⁴ |
| 0.050 | 0.532533 | 0.537388 | 9.12×10⁻³ |
| 0.100 | 0.559019 | 0.574866 | 2.83×10⁻² |

**Key finding:** Perturbation theory breaks down at λ ≥ 0.05 (error exceeds 
1%). At λ = 0.1, relative error reaches 2.83%, demonstrating that numerical 
diagonalization becomes essential in the strongly anharmonic regime.

![Anharmonic Analysis](plots/anharmonic_analysis.png)

**Physical significance:** This perturbative breakdown is directly analogous 
to failure of perturbative QCD at low energies (requiring Lattice QCD) and 
breakdown of weak-coupling approximations in quantum many-body systems.

## Repository Structure
qho.py - #main QHO solver

anharmonic_extension.py - #anharmonic extension + perturbation theory

plots/ - #all output figures



