---
layout: home
title: Overview
nav_order: 1
---

# Relapse-Driven Bistability in a Model of Gang Dynamics
{: .no_toc }

Companion site for the paper and its reproducibility code.
{: .fs-6 .fw-300 }

[Run it in Colab](#run-it-in-your-browser){: .btn .btn-primary .mr-2 }
[View on GitHub](https://github.com/machyman/kwofie2026relapse){: .btn }

---

> T. Kwofie, L. Rodriguez Rodriguez, X. Wang, J. M. Hyman and Y. Kang,
> *Relapse-Driven Bistability in a Nonlinear Model of Gang Dynamics with Deterrence
> and Reformed-Led Advocacy.*
> Mathematical Biosciences and Engineering (in revision, MBE-8616).

---

## The question

Suppose a community drives gang recruitment down far enough that, on paper, a gang
can no longer invade — each active member replaces less than one of themselves.
Should the gang disappear?

**Not necessarily.** This paper gives a mechanism for why not, and the mechanism is
**relapse**.

## The model

People are susceptible ($$s$$), active in a gang ($$g$$), or reformed ($$r$$), with
$$s + g + r = 1$$. They are recruited, they reform, and — crucially — they can
**relapse** from the reformed class back into the gang. Recruitment is damped by
perceived-risk deterrence ($$\epsilon$$) and by reformed-led advocacy ($$\eta$$):

$$
\frac{ds}{dt} = \mu - \frac{\beta_{sg}\, s\, g}{D} + \beta_{rs}\, r\, s - \mu s,
\qquad
\frac{dg}{dt} = \frac{\beta_{sg}\, s\, g}{D} + \beta_{rg}\, r\, g - (\gamma + \mu)\, g,
$$

$$
\frac{dr}{dt} = \gamma g - \beta_{rg}\, r\, g - \beta_{rs}\, r\, s - \mu r,
\qquad
D = 1 + \eta\, r\, s + \epsilon\, g .
$$

The invasion threshold is $$\mathcal{R}_0^{\rm simp} = \beta_{sg}/(\mu + \gamma)$$.

## The result

When relapse is strong enough, the bifurcation at $$\mathcal{R}_0^{\rm simp} = 1$$
turns **backward**. A stable gang-free state and a stable high-gang state then
coexist *below* the threshold, separated by a saddle whose stable manifold is the
basin boundary.

At the paper's illustrative parameters the bistable window spans an order of
magnitude — roughly $$0.09 < \mathcal{R}_0^{\rm simp} < 1$$. Inside it, **whether the
gang persists depends on where the community starts, not only on the parameters.**

> **Reducing recruitment below the invasion threshold is not sufficient to eliminate
> an established gang.** The system must be pushed back across the saddle.

## Where the leverage is

A global sensitivity analysis (PRCC, evaluated at the established-gang equilibrium so
the answer does not depend on the initial condition) gives a clear ranking:

| Parameter | PRCC | Reading |
|:--|--:|:--|
| Reformation $$\gamma$$ | $$-0.98$$ | Dominant, suppressive |
| Relapse $$\beta_{rg}$$ | $$+0.98$$ | Dominant, driving |
| Recruitment $$\beta_{sg}$$ | $$+0.17$$ | Smaller positive driver |
| Deterrence $$\epsilon$$, reintegration $$\beta_{rs}$$, advocacy $$\eta$$ | $$\approx 0$$ | Not resolved at the equilibrium |

Relapse and reformation are **co-dominant and opposite in sign** — the established
gang is held in place by the balance between them, which is the same mechanism that
bends the bifurcation backward. Recruitment matters for *invasion*, but is a weak
lever on an *established* gang, because susceptibles are nearly exhausted there.

Advocacy is not unimportant: it acts at first order on the **threshold** itself, which
is a different quantity from equilibrium prevalence. The paper keeps those two
questions separate, and so should any policy reading of this figure.

---

## Reproduce the results

```bash
git clone https://github.com/machyman/kwofie2026relapse.git
cd kwofie2026relapse
pip install -r requirements.txt
make verify        # re-derives every published number
```

Everything in the paper is reproducible from this repository. The Latin Hypercube
design and the bootstrap use a fixed seed, so the numbers come out the same every
time.

## Run it in your browser

No installation required.

| Notebook | | |
|:--|:--|:--|
| **The model and its equilibria** — including the saddle that a Newton solver would miss | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/machyman/kwofie2026relapse/blob/main/notebooks/01_model_and_equilibria.ipynb) |
| **Bifurcation diagrams** — Figures 3–7 | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/machyman/kwofie2026relapse/blob/main/notebooks/02_bifurcation_diagrams.ipynb) |
| **PRCC sensitivity** — Figure 8 and Table 4 | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/machyman/kwofie2026relapse/blob/main/notebooks/03_prcc_sensitivity.ipynb) |

## MATLAB

A complete MATLAB implementation lives in `matlab/`. It uses **only base MATLAB** —
no toolboxes — and runs unchanged under GNU Octave, so it is usable without a
Statistics Toolbox license:

```matlab
addpath('matlab');
sgrValidate();          % 10/10 checks against the published values
[E, kinds] = sgrEquilibria(sgrBaseline('beta_sg', 0.4325));
```

The two implementations are checked against each other on a shared sampling design:
they agree to $$3 \times 10^{-11}$$ on the PRCC values and $$4 \times 10^{-13}$$ on the
equilibria.

---

## How to cite

```bibtex
@article{kwofie2026relapse,
  title   = {Relapse-Driven Bistability in a Nonlinear Model of Gang Dynamics
             with Deterrence and Reformed-Led Advocacy},
  author  = {Kwofie, Theophilus and Rodriguez Rodriguez, Lucero and Wang, Xia
             and Hyman, James M. and Kang, Yun},
  journal = {Mathematical Biosciences and Engineering},
  year    = {2026},
  note    = {In revision (MBE-8616)}
}
```

Machine-readable citation metadata is in
[`CITATION.cff`](https://github.com/machyman/kwofie2026relapse/blob/main/CITATION.cff).
Released under the MIT License.

---

<div class="fs-3 fw-300 text-grey-dk-000" markdown="1">
Site built from commit
[`{{ site.data.version.commit }}`](https://github.com/machyman/kwofie2026relapse/commit/{{ site.data.version.commit }})
on {{ site.data.version.built }} · manuscript {{ site.data.version.manuscript_version }}
</div>
