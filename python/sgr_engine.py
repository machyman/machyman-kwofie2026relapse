"""
SGR gang-dynamics engine: model, equilibria, stability.

GENERATED FILE — do not edit by hand.
Source of truth: notebooks/01_model_and_equilibria.ipynb (re-run it to regenerate).

Companion paper: Kwofie, Rodriguez Rodriguez, Wang, Hyman and Kang,
"Relapse-Driven Bistability in a Nonlinear Model of Gang Dynamics with
Deterrence and Reformed-Led Advocacy," MBE-8616.
"""
from dataclasses import dataclass, replace

import numpy as np


@dataclass(frozen=True)
class Params:
    """Model parameters. Advocacy eta is independent of the reintegration rate."""
    beta_sg: float   # recruitment rate (week^-1)
    beta_rg: float   # relapse rate (week^-1)
    beta_rs: float   # reintegration, reformed -> susceptible (week^-1)
    mu: float        # demographic turnover (week^-1)
    gamma: float     # reformation rate (week^-1)
    eps: float       # perceived-risk deterrence (dimensionless)
    eta: float       # advocacy effectiveness (dimensionless)


#: Illustrative baseline (paper, Table 2). Not calibrated to a community.
ILLUSTRATIVE = Params(beta_sg=0.855, beta_rg=1.63, beta_rs=0.0001,
                      mu=0.005, gamma=0.86, eps=0.95, eta=0.0)


def r0_simp(p: Params) -> float:
    """Simplified reproduction number beta_sg / (mu + gamma)."""
    return p.beta_sg / (p.mu + p.gamma)


def denom(g: float, r: float, p: Params) -> float:
    """Recruitment-reduction denominator D = 1 + eta*r*s + eps*g."""
    s = 1.0 - g - r
    return 1.0 + p.eta * r * s + p.eps * g


def rhs(y, p: Params):
    """Reduced (g, r) right-hand side; s = 1 - g - r. ODE-solver ready."""
    g, r = y
    s = 1.0 - g - r
    D = denom(g, r, p)
    dg = p.beta_sg * s * g / D + p.beta_rg * r * g - (p.gamma + p.mu) * g
    dr = p.gamma * g - p.beta_rg * r * g - p.beta_rs * r * s - p.mu * r
    return np.array([dg, dr])


def jacobian(y, p: Params) -> np.ndarray:
    """Analytic Jacobian of the reduced system at y = (g, r)."""
    g, r = y
    s = 1.0 - g - r
    D = denom(g, r, p)
    dDg = -p.eta * r + p.eps
    dDr = p.eta * (s - r)
    dgg = p.beta_sg * ((-g + s) * D - s * g * dDg) / D**2 + p.beta_rg * r - (p.gamma + p.mu)
    dgr = p.beta_sg * ((-g) * D - s * g * dDr) / D**2 + p.beta_rg * g
    drg = p.gamma - p.beta_rg * r + p.beta_rs * r
    drr = -p.beta_rg * g - p.beta_rs * (s - r) - p.mu
    return np.array([[dgg, dgr], [drg, drr]])


def classify(y, p: Params) -> str:
    """'stable' | 'saddle' | 'unstable' from the Jacobian eigenvalues."""
    re = np.linalg.eigvals(jacobian(y, p)).real
    if np.all(re < -1e-9):
        return "stable"
    if np.any(re > 1e-9) and np.any(re < -1e-9):
        return "saddle"
    return "unstable"


def equilibria(p: Params, tol: float = 1e-9):
    """All feasible equilibria [(s, g, r), ...], by polynomial elimination.

    Deterministic: recovers every branch, including the saddle, with no initial
    guess. See the notebook for the derivation of the polynomial in r.
    """
    a, b, c = p.beta_sg, p.beta_rg, p.beta_rs
    m, gm, ep, et = p.mu, p.gamma, p.eps, p.eta
    k = gm + m

    # inner factor Dn + eta*r*S + eps*N, ascending powers of r
    i0 = gm
    i1 = (c - b) + et * gm + ep * (c + m)
    i2 = -et * (b + k) - ep * c
    i3 = et * b

    # P(r) = a*S(r) - inner(r) * (k - b*r)
    coeffs = [i3 * b,                      # r^4
              -(i3 * k - i2 * b),          # r^3
              a * b - (i2 * k - i1 * b),   # r^2
              -a * (b + k) - (i1 * k - i0 * b),
              a * gm - i0 * k]
    while len(coeffs) > 1 and abs(coeffs[0]) < 1e-14:
        coeffs = coeffs[1:]                # eta = 0 reduces the quartic to a cubic

    eqs = [(1.0, 0.0, 0.0)]                # the gang-free equilibrium always exists
    if len(coeffs) > 1:
        for root in np.roots(coeffs):
            if abs(root.imag) > 1e-8:
                continue
            r = float(root.real)
            if r < -tol or r > 1 + tol:
                continue
            Dn = gm + (c - b) * r
            if abs(Dn) < 1e-12:
                continue
            g = ((c + m) * r - c * r * r) / Dn
            s = 1.0 - g - r
            if g < tol or s < -tol:
                continue
            if not any(max(abs(np.array(e) - [s, g, r])) < 1e-6 for e in eqs):
                eqs.append((s, g, r))
    return sorted(eqs, key=lambda e: e[1])


def gstar(p: Params) -> float:
    """Gang fraction at the high-prevalence stable equilibrium E_+ (0 if none).

    Initial-condition independent: this is the sensitivity response used in the
    PRCC analysis, so the result does not depend on basin membership.
    """
    best = 0.0
    for (s, g, r) in equilibria(p):
        if g > 1e-6 and classify((g, r), p) == "stable" and g > best:
            best = g
    return best
