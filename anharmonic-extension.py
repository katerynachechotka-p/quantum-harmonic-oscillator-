"""
Anharmonic Oscillator Extension
================================
Extends the harmonic oscillator analysis by adding a quartic perturbation (λx⁴)
and comparing numerical results with first-order perturbation theory.

This demonstrates where perturbative approximations break down and numerical
methods become necessary.
"""

import numpy as np
import matplotlib.pyplot as plt

# Setup (same grid as harmonic oscillator)
L = 8.0
N = 300
x = np.linspace(-L, L, N)
dx = x[1] - x[0]

# Build operators
diag = -2.0 * np.ones(N)
off_diag = 1.0 * np.ones(N - 1)
D2 = (np.diag(diag) + np.diag(off_diag, k=1) + np.diag(off_diag, k=-1)) / dx**2
T = -0.5 * D2  # Kinetic energy

# Solve harmonic oscillator (reference case)
V_harmonic = np.diag(0.5 * x**2)
H_harmonic = T + V_harmonic
E_harmonic, psi_harmonic = np.linalg.eigh(H_harmonic)

# Normalize wavefunctions
for n in range(10):
    psi_harmonic[:, n] /= np.sqrt(np.sum(np.abs(psi_harmonic[:, n])**2) * dx)

print("Anharmonic Oscillator: V(x) = 1/2 x² + λx⁴")
print("=" * 50)

# Test different perturbation strengths
lambdas = [0.001, 0.01, 0.05, 0.1]
results = []

for lam in lambdas:
    # Solve anharmonic oscillator numerically
    V_anh = np.diag(0.5 * x**2 + lam * x**4)
    H_anh = T + V_anh
    E_anh, _ = np.linalg.eigh(H_anh)
    
    # Perturbation theory prediction
    # E ≈ E₀ + <ψ₀|λx⁴|ψ₀>
    psi0 = psi_harmonic[:, 0]
    correction = lam * np.sum(np.conjugate(psi0) * x**4 * psi0) * dx
    E_pert = E_harmonic[0] + correction
    
    # Compare
    error = abs(E_anh[0] - E_pert) / E_anh[0]
    results.append([lam, E_anh[0], E_pert, error])
    
    print(f"\nλ = {lam:.3f}")
    print(f"  Numerical:           {E_anh[0]:.6f}")
    print(f"  Perturbation theory: {E_pert:.6f}")
    print(f"  Relative error:      {error:.2e}")

print("\n" + "=" * 50)
print("Conclusion: Perturbation theory breaks down for λ ≥ 0.05")
print("This shows when numerical methods become necessary.")

# Visualization
results = np.array(results)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Energy comparison
ax1.plot(results[:, 0], results[:, 1], 'o-', label='Numerical', markersize=8)
ax1.plot(results[:, 0], results[:, 2], 's--', label='Perturbation Theory', markersize=8)
ax1.set_xlabel('λ')
ax1.set_ylabel('Ground State Energy')
ax1.set_title('Anharmonic Oscillator: Energy vs Perturbation Strength')
ax1.legend()
ax1.grid(True, alpha=0.3)

# Error
ax2.semilogy(results[:, 0], results[:, 3], 'o-', markersize=8)
ax2.set_xlabel('λ')
ax2.set_ylabel('Relative Error')
ax2.set_title('Perturbation Theory Accuracy')
ax2.grid(True, alpha=0.3)
ax2.axhline(y=0.01, color='r', linestyle='--', alpha=0.5, label='1% error')
ax2.legend()

plt.tight_layout()
plt.savefig('anharmonic_analysis.png', dpi=150)
print("\nPlot saved as 'anharmonic_analysis.png'")
plt.show()
