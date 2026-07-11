"""Engine tests: the model's published equilibria, stability, and invariants."""
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "python"))

from sgr_engine import ILLUSTRATIVE, replace, equilibria, classify, rhs, r0_simp  # noqa: E402


def test_gang_free_equilibrium_exists():
    E = equilibria(ILLUSTRATIVE)
    assert np.max(np.abs(np.array(E[0]) - [1.0, 0.0, 0.0])) < 1e-9


def test_worked_instance_three_equilibria():
    """At R0 = 0.5 the model has three equilibria, including the saddle."""
    p = replace(ILLUSTRATIVE, beta_sg=0.4325)
    E = equilibria(p)
    ref = np.array([[1.0, 0.0, 0.0],
                    [0.627830769838, 0.006983477402, 0.365185752760],
                    [0.035925517346, 0.440121503665, 0.523952978989]])
    assert len(E) == 3
    assert np.max(np.abs(np.array(E) - ref)) < 1e-6


def test_bistability_below_threshold():
    """The paper's central claim: a stable gang state persists for R0 < 1."""
    p = replace(ILLUSTRATIVE, beta_sg=0.4325)
    E = equilibria(p)
    kinds = [classify((g, r), p) for (s, g, r) in E]
    assert r0_simp(p) < 1.0
    assert kinds[0] == "stable" and kinds[-1] == "stable" and "saddle" in kinds


def test_residuals_at_machine_precision():
    for (s, g, r) in equilibria(ILLUSTRATIVE):
        assert np.max(np.abs(rhs((g, r), ILLUSTRATIVE))) < 1e-10


def test_advocacy_is_decoupled():
    """eta must act independently of beta_rs (post-reformulation model)."""
    p0 = replace(ILLUSTRATIVE, eta=0.0)
    p1 = replace(ILLUSTRATIVE, eta=1.0)
    assert p0.beta_rs == p1.beta_rs
    assert abs(rhs((0.3, 0.4), p0)[0] - rhs((0.3, 0.4), p1)[0]) > 1e-6
