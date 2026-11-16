# Spacetime Topological Emergent Micro-Solitons (STEMS) — Numerical Solver & Verification Suite

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24283439.svg)](https://doi.org/10.5281/zenodo.24283439)
![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)
![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)

Official computational reproduction suite, canonical Robin-regularized Lobatto IIIA boundary value solver, and radial soliton diagnostic codebase for the preprint:

> **Spacetime Topological Emergent Micro-Solitons (STEMS): Derrick-Bypass via Quartic Axial Weitzenböck Torsion, Covariant Tetrad Reduction, and Integer Charge Quantization**
> Author: Dev S. Patel (Version 2.0, November 2026)
> DOI: [10.5281/zenodo.24283439](https://doi.org/10.5281/zenodo.24283439)

---

## 🔭 Overview & Key Physical Results

This repository contains the numerical verification script and multi-panel diagnostic figures establishing the existence, stability, core regularity, and topological charge quantization of **Spacetime Topological Emergent Micro-Solitons (STEMS)**—a classical non-linear teleparallel field framework deriving localized, singularity-free particle-like field configurations directly from spacetime frame geometry:

* **Canonical Cellular Scale & Quartic Axial Torsion Action:** Formulated within the covariant Teleparallel Equivalent of General Relativity (TEGR) augmented by a quartic completely antisymmetric axial Weitzenböck torsion invariant $\mathcal{T}_{\mathrm{axial}} \equiv (A_{\mu}A^{\mu})^{2}$, with coupling $\alpha_{4} \sim l_{\mathrm{cell}}^{2}$ anchored to the Bekenstein–Hawking calibrated cellular scale $l_{\mathrm{cell}} \equiv 2\sqrt{\ln 2}\,l_{P} \approx 2.6914 \times 10^{-35}\,\mathrm{m}$.
* **1D Isotropic Derrick-Bypass & Teleparallel Virial Equilibrium:** Under 1-parameter active spatial dilation $x^{i} \mapsto \lambda x^{i}$, the scaled static energy functional $E(\lambda) = \lambda^{-1}I_{2} + \lambda I_{4}$ stabilizes at $\lambda = 1$ with exact teleparallel virial equilibrium $I_{2} = I_{4} \equiv E_{0} > 0$ and strictly positive breathing curvature $\frac{d^{2}E}{d\lambda^{2}}\big|_{\lambda=1} = 2E_{0} > 0$.
* **2-Parameter Anisotropic Hessian Stability:** Extending active dilation to independent longitudinal ($\lambda_{L}$) and transverse ($\lambda_{T}$) deformations yields a strictly positive-definite $2 \times 2$ Hessian matrix with $\mathrm{Tr}(\mathcal{H}) = \frac{86}{15}E_{0} > 0$, $\det(\mathcal{H}) = \frac{208}{45}E_{0}^{2} > 0$, and eigenvalues $\mu_{\pm} = \frac{43 \pm \sqrt{809}}{15}E_{0}$ ($\mu_{+} \approx 4.7629\,E_{0}$, $\mu_{-} \approx 0.9705\,E_{0}$), ruling out quadrupole shear and squashing instabilities.
* **2D Cylindrical Flux Tubes & Cornell Confining Potential:** Transverse Hobart–Derrick virial equilibrium ($E_{\perp}^{(4)} = E_{\perp}^{(0)}$) combined with zeta-regularized zero-point string fluctuations ($\zeta(-1) = -1/12$) yields the linear-plus-Lüscher Cornell confining potential $V_{\mathrm{Cornell}}(R) = \sigma_{\mathrm{torsion}}R - \frac{\pi\hbar c}{12R} + V_{0}$ for open fractional meron defects.
* **Off-Shell Charge Quantization & Faddeev–Bogomolny Bound:** Under Alexandroff compactification $\mathbb{R}^{3} \cup \{\infty\} \cong S^{3}$, spatial frame twists classified by $\pi_{3}(\mathrm{SO}(3)) \cong \pi_{3}(S^{3}) \cong \mathbb{Z}$ yield an identically conserved off-shell topological 4-current ($\partial_{\mu}J^{\mu} \equiv 0$), integer electric charge quantization $q = Q \cdot e_{0}$, a strict Faddeev–Bogomolny lower energy bound $E[e] \ge \frac{3\pi^{2}c^{4}\sqrt{\alpha_{4}}}{2G_{N}}|Q|$, and emergent far-field Coulomb electrostatics $E_{r}(r) = \frac{Q e_{0}}{4\pi\epsilon_{0}r^{2}}[1 - \mathcal{O}(r^{-6})]$.
* **Frobenius Core Regularity & Von Laue Mechanical Equilibrium:** Parity-regular Frobenius expansion $F(r) = 2\pi - c_{1}r + c_{3}r^{3} + \mathcal{O}(r^{5})$ cancels $\mathcal{O}(1/r)$ centrifugal poles identically, yielding finite central energy density $\rho_{E}(0) = \frac{c^{4}}{16\pi G_{N}}\left(3c_{1}^{2} + \frac{5}{16}\alpha_{4}c_{1}^{4}\right) \approx 0.1353\,\frac{c^{4}}{16\pi G_{N}}$ at $c_{1} = 0.2118$, finite ADM mass $M_{\mathrm{ADM}} = \frac{c^{2}}{4G_{N}}\sqrt{\alpha_{4}}\tilde{M}_{0}$ ($\tilde{M}_{0} \approx 73.42$), and exact Von Laue radial pressure balance $\int_{0}^{\infty}[p_{r}(r) + 2p_{\perp}(r)]r^{2}\,dr = 0$.
* **Ostrogradsky Ghost-Freedom & Caustic-Free Cauchy Hyperbolicity:** Temporal decomposition confirms zero second-order time accelerations ($\mathcal{H}_{\mathrm{STEMS}} \ge 0$), while the effective acoustic inverse metric $G_{\mathrm{eff}}^{\mu\nu} = g^{\mu\nu} + 8\alpha_{4}(A_{\lambda}A^{\lambda})A^{\mu}A^{\nu}$ remains strictly hyperbolic and subluminal well below the caustic threshold ($|A_{\lambda}A^{\lambda}|_{\max}/A_{\mathrm{crit}}^{2} \approx 0.0362 < 1$ at rest; $\le 0.142$ across boosted modes), alongside exact post-Newtonian agreement ($\gamma_{\mathrm{PPN}} = \beta_{\mathrm{PPN}} = 1$).

---

## 💻 Terminal Execution Commands

### 0. Clone and Install Dependencies

```bash
git clone [https://github.com/dev-s-patel/STEMS-Solvers-and-Figures.git](https://github.com/dev-s-patel/STEMS-Solvers-and-Figures.git)
cd STEMS-Solvers-and-Figures
pip install -r requirements.txt
```

### 1. Appendix A: Canonical Robin Regularization & Lobatto IIIA Collocation Solver

```bash
python run_stems_soliton_solver.py
```

* **Script:** `run_stems_soliton_solver.py`
* **Function:** Numerically solves the non-linear radial Euler–Lagrange boundary value problem for the ground-state STEMS hedgehog soliton ($\alpha_{4} = 1.0$) over $r \in [10^{-2}, 100.0]$ using 4th-order Lobatto IIIA collocation (`scipy.integrate.solve_bvp`, `tol=1e-8`, `max_nodes=100000`) subject to Canonical Inner ($F(\epsilon) - \epsilon F'(\epsilon) = 2\pi$) and Outer ($r_{\max}F'(r_{\max}) + 2F(r_{\max}) = 0$) Robin boundary conditions.
* **Verification:** Converges across adaptive collocation meshes ($N = 500 \to 8000$) and inner cutoffs ($\epsilon = 10^{-1} \to 10^{-4}$) to relaxed core slope $c_{1} = 0.2118$, central normalized energy density $\tilde{\rho}_{E}(0) = 0.1353$, exact teleparallel virial equilibrium $I_{2}/I_{4} = 1.000000$, and integer topological winding degree $Q = 1.00000$ ($r_{1/2} \approx 8.1\sqrt{\alpha_{4}}$, asymptotic dipole constant $C_{\infty} \approx 1.76 \times 10^{3}\alpha_{4}$).
* **Outputs:**
  * `STEMS_Soliton_Diagnostics.png` (High-resolution 300 DPI 4-panel diagnostic figure: radial chiral profile $F(r)$, core energy density $\tilde{\rho}_{E}(r)$, cumulative topological charge $Q(r)$, and Von Laue radial pressure equilibrium $r^{2}p_{r}(r)$).

---

## 📦 Requirements

* Python $\ge 3.8$
* `numpy >= 1.24.0`
* `scipy >= 1.10.0`
* `matplotlib >= 3.7.0`

Install all dependencies via `pip`:

```bash
pip install numpy scipy matplotlib
```

---

## 📜 Citation

If you use or verify this computational solver in your research, please cite the STEMS manuscript:

```bibtex
@article{patel2026stems,
  title   = {Spacetime Topological Emergent Micro-Solitons (STEMS): Derrick-Bypass via Quartic Axial Weitzenb{\"o}ck Torsion, Covariant Tetrad Reduction, and Integer Charge Quantization},
  author  = {Patel, Dev S.},
  year    = {2026},
  month   = {nov},
  doi     = {10.5281/zenodo.24283439},
  url     = {[https://doi.org/10.5281/zenodo.24283439](https://doi.org/10.5281/zenodo.24283439)}
}
```

---

## ⚖️ License

This project is open-source and distributed under the [MIT License](LICENSE).
