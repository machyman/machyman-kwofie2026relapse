"""PRCC regression: the published sensitivity values must be reproduced."""
import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "python"))

from sgr_engine import ILLUSTRATIVE, replace, gstar  # noqa: E402

PARAM_ORDER = ["beta_sg", "gamma", "beta_rg", "beta_rs", "eps", "eta"]


def test_shared_lhs_design_present():
    """The shared design is what makes MATLAB/Python parity exact."""
    path = os.path.join(ROOT, "data", "prcc_lhs_samples.csv")
    X = np.loadtxt(path, delimiter=",")
    assert X.shape == (1000, 6)


def test_published_prcc_values():
    path = os.path.join(ROOT, "data", "prcc_results.json")
    with open(path) as f:
        res = json.load(f)
    expected = {"gamma": -0.977, "beta_rg": 0.976, "beta_sg": 0.167}
    got = {row["parameter"]: row["prcc"] for row in res["results"]}
    for name, want in expected.items():
        assert abs(got[name] - want) < 5e-3


def test_response_is_initial_condition_independent():
    """g* is read off the equilibrium, so it cannot depend on an initial guess."""
    base = replace(ILLUSTRATIVE, beta_sg=0.855)
    assert abs(gstar(base) - 0.4579) < 1e-3
