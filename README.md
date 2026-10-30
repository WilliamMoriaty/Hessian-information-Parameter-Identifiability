# Hessian Information for Parameter Identifiability

Code and data accompanying the paper on **Hessian-corrected Fisher information for parameter identifiability** in data-driven models. The repository reproduces the fitting, identifiability analysis, and uncertainty quantification for three case studies: an HIV viral-dynamics model, the JAK-STAT signaling benchmark, and a network model of amyloid-β (Aβ) spreading.

<img width="13813" height="7440" alt="Fig1" src="https://github.com/user-attachments/assets/a6491860-85c0-4827-acff-4d6089382e27" />

## Overview

Classical local identifiability analysis relies on the Fisher information matrix (FIM) $F = J^\top J$, which uses only first-order sensitivities. When $F$ is singular or ill-conditioned, or when the fit leaves a non-zero residual, first-order information alone can misclassify parameters and misstate their uncertainty. This work augments the FIM with second-order (Hessian) information of the observation map.

The workflow, shown in the figure above, has five steps.

**1. Parameter fitting.** Given a model $\varphi(t;\theta)$, an observation $h(\varphi(t;\theta))$ and data $\{t_n, \hat h_n\}_{n=1}^N$, parameters are estimated by least squares in log coordinates $\eta = \ln\theta$:

$$
\eta^* = \arg\min_{\eta\in\Theta} \ell(\eta;\hat h),\qquad
\ell(\eta;\hat h) = \tfrac12\sum_{n=1}^N \lVert h(t_n,\eta)-\hat h_n\rVert_2^2,\qquad
r^* = h(t_n,\eta^*)-\hat h_n .
$$

**2. FIM and Hessian tensor.** The FIM $F$ and the Hessian tensor $H$ (second directional response of $h$) are computed from forward first- and second-order sensitivity equations. An eigendecomposition $F = Q\Lambda Q^\top$, $Q = [U\ V]$, separates range $F$ (columns $U$) from ker $F$ (columns $V$).

**3. Parameter identifiability: Lipschitz vs. Hölder ($r^*=0$).** A parameter is identifiable if $\sup_{\eta_\delta\in M_\delta}\lVert\eta_\delta-\eta^*\rVert_2 \le C\,\omega(\lVert\delta\rVert_2)$ for data perturbations $\delta$.
- *Theorem 1:* if $F \succ 0$, then $\omega(\varepsilon)=\varepsilon$ (Lipschitz stability).
- If $F$ is singular, the error splits as $\Delta\eta = \Delta\eta_R + \Delta\eta_N$. The range component satisfies $\lVert\Delta\eta_R\rVert_2 \le \varepsilon/\sigma_+$, so first-order information suffices. The kernel component is unconstrained by $F$ alone but controlled by $H$: $\lVert\Delta\eta_N\rVert_2 \le 4\varepsilon^{1/2}/\sqrt{c_2}$ with $c_2 = \min_{\beta\in\ker F,\,\lVert\beta\rVert_2=1}\lVert\Pi^\perp H[\beta,\beta]\rVert_2 > 0$ and $\Pi^\perp = I - JF^\dagger J^\top$ (Hölder-1/2 stability).

**4. Hessian correction of the Fisher information.** For non-zero residuals, the corrected curvature is

$$K = F + r^*\otimes H .$$

The residual term lifts small eigenvalues of $F$ along ker $F$. *Theorem 3* gives a crossover at $\varepsilon_* = \gamma\nu/c_2$, where $\nu = \lVert r^*\rVert_2$ and $\gamma = \min_{\xi\in\ker F,\,\lVert\xi\rVert_2=1}\langle\hat r, \Pi^\perp H[\xi,\xi]\rangle$ is the residual curvature on ker $F$: stability is Lipschitz for $\varepsilon\le\varepsilon_*$ and Hölder-1/2 for $\varepsilon\ge\varepsilon_*$.

**5. Coordinate identifiability and uncertainty quantification.** Partitioning $K$ with respect to coordinate $p$ gives the profiled curvature $\kappa_p = K_{pp} - k_p^\top\tilde K^{-1}k_p$ (*Theorem 4*); coordinate $\eta_p$ is identifiable iff $\kappa_p > 0$, with Lipschitz coefficient $L_p$ for the profiled loss. Uncertainty is then propagated with the quadratic posterior $\eta\mid\hat h \sim \mathcal N(\eta^*, \Sigma_\eta)$, comparing $K$-based and $F$-based intervals for parameters and model predictions (and against MCMC).

## Repository structure

```
.
├── HIV/                         HIV viral dynamics (patients 3 and 8)
│   ├── HIV.ipynb                fitting + identifiability, patient 3
│   ├── HIV8.ipynb               fitting + identifiability, patient 8
│   ├── HIV_patient3.txt         viral load data
│   ├── HIV_patient8.txt
│   ├── HIV_patient3_results.mat results exported by the notebooks
│   ├── HIV_patient8_results.mat
│   ├── identifiability.m        figures, patient 3
│   └── identifiability_8.m      figures, patient 8
├── JAK-STAT-Benchmark/          JAK-STAT pathway (Raue et al. 2009 benchmark)
│   ├── Pathway.ipynb            fitting, profile likelihood, identifiability, crossover
│   ├── data.mat                 Swameye et al. (2003) dataset 1
│   ├── fit_results.mat
│   ├── profile_likelihood.mat
│   ├── identifiability_analysis.mat
│   ├── crossover_analysis.mat
│   ├── identifiability.m        identifiability figures
│   └── crossover_regime.m       Lipschitz vs. Hölder-1/2 crossover figure
└── Abeta/                       Aβ spreading on a 68-region brain network
    ├── Abeta_CN_interior_recon.ipynb
    ├── Lp68.mat                 graph Laplacian (68 ROIs)
    ├── Patient_751_5.mat        patient Aβ data
    ├── Abeta_CN_interior_recon_analysis_results.mat
    └── identifiability.m        figures
```

## Case studies

**HIV viral dynamics.** The target-cell-limited model
$\dot T=\lambda-dT-kVT,\ \dot I=kVT-\delta I,\ \dot V=pI-cV$
with observation $h=\log_{10}(2V\times10^3)$. With $c=3$ fixed, the five free parameters $\eta=\log(\lambda,d,k,\delta,p)$ are fitted to patient viral-load data. Outputs include spectra and eigenvectors of $F_\eta$ and $K_\eta$, eigenvector perturbation maps, residual-Hessian indices, and $K$/$F$-based credible intervals compared with MCMC.

**JAK-STAT signaling (benchmark).** The four-state STAT model of Raue et al. (2009) with a 10-step linear-chain delay, driven by spline-interpolated pEpoR input and fitted to Swameye et al. (2003) data. The benchmark contains known structural non-identifiabilities ($p_2, x_1(0), s_1, s_2$) and a practical non-identifiability ($p_3$), which are used to validate the Hessian-based diagnostics against profile likelihoods. `crossover_regime.m` visualizes the Lipschitz-to-Hölder crossover of the weak $\ln p_2$ direction.

**Aβ spreading.** A network Fisher–KPP model $\dot x = -Lx + \lambda\odot x\odot(K-x)$ on a 68-region brain graph, with $\theta=(\lambda_1,\dots,\lambda_{68},K_1,\dots,K_{68})$ (136 parameters) fitted to longitudinal regional Aβ PET data. Rates and carrying capacities are optimized alternately.

## Requirements

**Python** (3.9+): `numpy`, `scipy`, `matplotlib`, plus `h5py` and `fides` for the Aβ notebook.

```bash
pip install numpy scipy matplotlib h5py fides
```

**MATLAB** (R2020b or later recommended, for `tiledlayout`) for the figure scripts.

## Usage

Each case study follows the same two-step pattern: run the notebook to fit the model and export results to `.mat`, then run the MATLAB script to produce figures.

```bash
git clone https://github.com/WilliamMoriaty/Hessian-information-Parameter-Identifiability.git
cd Hessian-information-Parameter-Identifiability/HIV
jupyter notebook HIV.ipynb        # run all cells in order
```

```matlab
% in MATLAB, from the same folder
identifiability          % patient 3
identifiability_8        % patient 8 (after running HIV8.ipynb)
```

The exported `.mat` result files are already included, so the MATLAB figures can be reproduced without rerunning the notebooks. Note that the MCMC sections (especially for the 136-parameter Aβ model) can take substantial time; check convergence diagnostics before interpreting the resulting intervals.


## License

Please add a license (e.g. MIT) before reuse.
