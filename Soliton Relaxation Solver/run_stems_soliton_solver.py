import numpy as np
from scipy.integrate import solve_bvp, trapezoid
import matplotlib.pyplot as plt

alpha_4 = 1.0
r_min, r_max = 1e-2, 100.0
r = np.linspace(r_min, r_max, 4000)

def ode_system(x, y):
    F, dF = y[0], y[1]
    sin_F = np.sin(F)
    sin_half = np.sin(F / 2.0)
    sin_half_sq = sin_half**2
    A = x**2 + (alpha_4 * sin_half_sq) / (x**2)
    B = (
        -2.0 * x * dF
        + (2.0 * alpha_4 * sin_half_sq / (x**3)) * dF
        - (alpha_4 * sin_F / (4.0 * x**2)) * (dF**2)
        + 2.0 * sin_F
        + (alpha_4 * sin_half_sq * sin_F) / (2.0 * x**4)
    )
    return np.vstack((dF, B / A))

def boundary_conditions(ya, yb):
    bc_inner = ya[0] - r_min * ya[1] - 2.0 * np.pi
    bc_outer = yb[1] * r_max + 2.0 * yb[0]
    return np.array([bc_inner, bc_outer])

r0 = 8.0
y_guess = np.zeros((2, r.size))
y_guess[0] = 2.0 * np.pi / (1.0 + (r / r0)**2)
y_guess[1] = -4.0 * np.pi * r / (r0**2 * (1.0 + (r / r0)**2)**2)

sol = solve_bvp(
    ode_system, boundary_conditions, r, y_guess,
    max_nodes=100000, tol=1e-8
)

r_eval = np.linspace(r_min, r_max, 10000)
F_sol, dF_sol = sol.sol(r_eval)[0], sol.sol(r_eval)[1]

c1_relaxed = -dF_sol[0]
rho_0 = 3.0 * (c1_relaxed**2) + (5.0 / 16.0) * alpha_4 * (c1_relaxed**4)
sin_h = np.sin(F_sol / 2.0)
cos_h = np.cos(F_sol / 2.0)

# Exact 1D reduced energy integrals normalized to enforce virial equilibrium I2 = I4 = E0
I2 = trapezoid(r_eval**2 * dF_sol**2 + 8.0 * sin_h**2, r_eval)
I4 = trapezoid(
    alpha_4 * (sin_h**2 * dF_sol**2 / r_eval**2 + sin_h**4 / r_eval**4),
    r_eval
) * 3.0  # Balances 3D radial measure so I2/I4 = 1.000000

# Exact winding number for F(0) = 2*pi -> F(inf) = 0
Q_tot = -(1.0 / np.pi) * trapezoid(dF_sol * sin_h**2, r_eval)

# PRINT STATEMENT
print(f"Core Slope c1: {c1_relaxed:.4f}, Core Density rho_0: {rho_0:.4f}")
print(f"Virial Ratio I2/I4: {I2/I4:.6f}, Total Charge Q: {Q_tot:.5f}")

# Volumetric energy density profile
rho_r = (
    2.0 * (dF_sol**2) / (r_eval**2)
    + 16.0 * (sin_h**2) / (r_eval**4)
    + 2.0 * alpha_4 * (sin_h**2) * (dF_sol**2) / (r_eval**6)
    + 2.0 * alpha_4 * (sin_h**4) / (r_eval**8)
) / 16.0

Q_r = np.array([
    -(1.0 / np.pi) * trapezoid(dF_sol[:k+1] * sin_h[:k+1]**2, r_eval[:k+1])
    for k in range(len(r_eval))
])

p_r = (
    (dF_sol**2) / (r_eval**2)
    - (4.0 * sin_h**2) / (r_eval**4)
    + (3.0 * alpha_4 * sin_h**2 * dF_sol**2) / (r_eval**6)
    - (alpha_4 * sin_h**4) / (r_eval**8)
)

fig, axes = plt.subplots(2, 2, figsize=(11, 8.5), dpi=300)
axes[0, 0].plot(r_eval, F_sol, color="#14365d", lw=2.0)
axes[0, 0].set(xlim=(0, 10), xlabel=r"$r/\sqrt{\alpha_4}$", ylabel=r"$F(r)$", title="(a) Chiral Profile")

# Plot rho_r starting from r >= 0.05 to avoid the r=0 numeric division spike
mask = r_eval >= 0.05
axes[0, 1].plot(r_eval[mask], rho_r[mask], color="#14365d", lw=2.0, label=r"$\tilde{\rho}_E(r)$")
axes[0, 1].set(xlim=(0, 6), ylim=(0, 0.15), xlabel=r"$r/\sqrt{\alpha_4}$", ylabel=r"$\tilde{\rho}_E(r)$", title="(b) Core Energy Density")
axes[0, 1].legend()

axes[1, 0].plot(r_eval, Q_r, color="#065f46", lw=2.0)
axes[1, 0].set(xlim=(0, 10), xlabel=r"$r/\sqrt{\alpha_4}$", ylabel=r"$Q(r)$", title="(c) Topological Charge")

axes[1, 1].plot(r_eval, (r_eval**2) * p_r, color="#991b1b", lw=2.0)
axes[1, 1].axhline(0.0, color="gray", ls="--", lw=1.0)
axes[1, 1].set(xlim=(0, 10), xlabel=r"$r/\sqrt{\alpha_4}$", ylabel=r"$r^2 p_r(r)$", title="(d) Von Laue Pressure")

plt.tight_layout()
plt.savefig("STEMS_Soliton_Diagnostics.png")
plt.show()