# Relapse-Driven Bistability in a Model of Gang Dynamics

Reproducibility code for:

> T. Kwofie, L. Rodriguez Rodriguez, X. Wang, J. M. Hyman and Y. Kang,
> *Relapse-Driven Bistability in a Nonlinear Model of Gang Dynamics with
> Deterrence and Reformed-Led Advocacy.* Mathematical Biosciences and
> Engineering (in revision, MBE-8616).

**[Read the companion site →](https://machyman.github.io/kwofie2026relapse/)**

Everything the paper reports — every figure, every equilibrium, every sensitivity
value — is produced by the code here, in **Python (Colab notebooks)** and in
**MATLAB**. No empirical individual-level data are used; the numerical results
follow from the model and an illustrative parameter set.

---

## The result in one paragraph

A susceptible–gang–reformed (SGR) model in which people are recruited into a gang,
reform out of it, and can **relapse** back. Recruitment is damped by perceived-risk
deterrence and by reformed-led advocacy. When relapse is strong enough, the
bifurcation at the invasion threshold turns **backward**: a stable gang-free state
and a stable high-gang state coexist *below* the threshold, separated by a saddle.
The practical consequence is that **pushing the invasion threshold below one is not
sufficient to eliminate an established gang** — the system has to be moved across
the saddle. The sensitivity analysis says where the leverage is: relapse and
reformation are co-dominant and opposite in sign, while recruitment is a
comparatively small driver of the *established* gang, because susceptibles are
depleted there.

$$\frac{dg}{dt} = \frac{\beta_{sg}\, s\, g}{1 + \eta\, r\, s + \epsilon\, g}
  + \beta_{rg}\, r\, g - (\gamma + \mu)\, g$$

---

## Quick start

```bash
git clone https://github.com/machyman/kwofie2026relapse.git
cd kwofie2026relapse
pip install -r requirements.txt

make verify      # reproduce the published results  (~1 min)
make test        # run the Python test suite
make figures     # regenerate every figure by executing the notebooks
make matlab      # run the MATLAB / Octave validation suite
```

### Run it in the browser (no install)

| Notebook | What it does | |
|---|---|---|
| `01_model_and_equilibria.ipynb` | The model, the equilibria, and the verification suite | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/machyman/kwofie2026relapse/blob/main/notebooks/01_model_and_equilibria.ipynb) |
| `02_bifurcation_diagrams.ipynb` | Figures 3–7: the backward bifurcation and the intervention sweeps | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/machyman/kwofie2026relapse/blob/main/notebooks/02_bifurcation_diagrams.ipynb) |
| `03_prcc_sensitivity.ipynb` | Figure 8 and Table 4: the PRCC sensitivity analysis | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/machyman/kwofie2026relapse/blob/main/notebooks/03_prcc_sensitivity.ipynb) |

---

## Repository structure

| Path | Purpose |
|---|---|
| `notebooks/01_model_and_equilibria.ipynb` | **Source of truth for the Python implementation.** Defines the model and writes `python/sgr_engine.py`. Runs a 10-check verification suite. |
| `notebooks/02_bifurcation_diagrams.ipynb` | Reproduces Figures 3–7. Verifies the transcritical point, the saddle-node locations, and the width of the bistable interval. |
| `notebooks/03_prcc_sensitivity.ipynb` | Reproduces Figure 8 and Table 4. Writes the shared Latin Hypercube design used for cross-language parity. |
| `python/sgr_engine.py` | **Generated** by notebook 01 — do not edit by hand. Model, equilibria, stability. Imported by notebooks 02–03, the tests, and `verify.py`. |
| `matlab/sgrBaseline.m` | Illustrative baseline parameters. |
| `matlab/sgrRHS.m` | Right-hand side of the reduced model (ODE-solver ready). |
| `matlab/sgrJacobian.m` | Analytic Jacobian of the reduced system. |
| `matlab/sgrEquilibria.m` | All feasible equilibria, by polynomial elimination, with stability classification. |
| `matlab/sgrBifurcation.m` | Parameter sweep returning equilibrium branches and saddle-node locations. |
| `matlab/sgrPRCC.m` | PRCC sensitivity analysis, written toolbox-free (see below). |
| `matlab/sgrValidate.m` | 10-check MATLAB validation suite; nonzero exit code on failure. |
| `tests/test_engine.py` | Equilibria, bistability, residuals, and the advocacy-decoupling invariant. |
| `tests/test_prcc.py` | Regression against the published PRCC values and the shared design. |
| `data/prcc_lhs_samples.csv` | The Latin Hypercube design (1000 × 6, seed 42). Shared so MATLAB and Python evaluate identical samples. |
| `data/prcc_results.json` | Published PRCC values with bootstrap confidence intervals. |
| `verify.py` | Reproducibility entry point. Re-derives the published results and runs the tests. |
| `Makefile` | `verify`, `test`, `figures`, `matlab` targets. |
| `requirements.txt` | Pinned Python dependencies. |
| `docs/` | Source of the companion GitHub Pages site. |
| `.github/workflows/pages.yml` | Builds and deploys the site. |
| `CITATION.cff`, `LICENSE` | Citation metadata; MIT license. |

---

## How the implementations stay consistent

Three implementations of one model is a standing invitation to drift. Two
conventions prevent it.

1. **One definition of the model.** Notebook 01 is the authored source; its engine
   cell writes `python/sgr_engine.py`. Notebooks 02 and 03, the tests, and
   `verify.py` all import that file. Nothing redefines the model.

2. **Parity is measured, not assumed.** MATLAB and Python are compared on the same
   Latin Hypercube design (`data/prcc_lhs_samples.csv`), because seeded random
   streams differ between the two languages and independently generated samples
   could never agree exactly. On the shared design the two implementations agree
   to **3 × 10⁻¹¹** across all six PRCC values, and the equilibria agree to
   **4 × 10⁻¹³**.

Equilibria are found by **polynomial elimination**, not by iteration. This is not a
stylistic choice: the model is bistable, and a Newton solver seeded from a single
initial guess silently misses the **saddle** branch — which is precisely the branch
the paper is about.

## MATLAB without toolboxes

The MATLAB code uses only base MATLAB (`roots`, `eig`, `ode45`) and runs unchanged
under **GNU Octave**. Every Statistics-Toolbox function the analysis would normally
need is implemented from scratch: `lhsdesign`, `tiedrank`, `partialcorr`, `bootstrp`
and `corr` are replaced by local equivalents in `matlab/sgrPRCC.m`.

> Note for anyone verifying this: `corr` *exists* in Octave but is a Statistics
> Toolbox function in MATLAB. Running cleanly under Octave therefore does **not**
> prove toolbox independence; it has to be checked against MATLAB's toolbox list.

---

## Reproducibility notes

- Deterministic: the Latin Hypercube design and the bootstrap use a fixed seed (42),
  so `make verify` reproduces the published numbers exactly.
- The notebooks carry a `FULL` switch. `FULL = False` (the default) completes in
  well under five minutes; `FULL = True` produces the publication-resolution figures.
- The figures in the article itself are the authors' original EPS files. Re-running
  the notebooks produces equivalent diagrams — identical branch structure and
  saddle-node values — but not byte-identical images across matplotlib versions.

## License and citation

Released under the MIT License (`LICENSE`). If you use this code, please cite the
paper and the software; citation metadata is in `CITATION.cff`.
