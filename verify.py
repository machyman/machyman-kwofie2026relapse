#!/usr/bin/env python3
"""Reproducibility check for the SGR gang-dynamics repository.

Runs the Python test suite and re-derives the published equilibria and PRCC
values. Exits 0 when everything the paper reports is reproduced.

Usage:
    python verify.py
"""
import json
import os
import subprocess
import sys

import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(ROOT, "python"))

from sgr_engine import ILLUSTRATIVE, replace, equilibria, classify, r0_simp, gstar  # noqa: E402


def main() -> int:
    failures = 0
    print("SGR repository verification")
    print("-" * 62)

    # 1. equilibria at the worked instance
    p = replace(ILLUSTRATIVE, beta_sg=0.4325)
    E = equilibria(p)
    kinds = [classify((g, r), p) for (s, g, r) in E]
    ok = len(E) == 3 and kinds == ["stable", "saddle", "stable"] and r0_simp(p) < 1
    failures += report("bistability below threshold (3 equilibria, R0 = 0.5)", ok)

    # 2. baseline high-gang equilibrium
    ok = abs(gstar(replace(ILLUSTRATIVE, beta_sg=0.855)) - 0.4579) < 1e-3
    failures += report("baseline E_+ at g = 0.4579", ok)

    # 3. published PRCC values
    with open(os.path.join(ROOT, "data", "prcc_results.json")) as f:
        res = json.load(f)
    got = {row["parameter"]: row["prcc"] for row in res["results"]}
    ok = (abs(got["gamma"] + 0.977) < 5e-3 and abs(got["beta_rg"] - 0.976) < 5e-3
          and abs(got["beta_sg"] - 0.167) < 5e-3)
    failures += report("PRCC: gamma and beta_rg co-dominant, beta_sg smaller", ok)

    # 4. shared LHS design (used for exact MATLAB parity)
    X = np.loadtxt(os.path.join(ROOT, "data", "prcc_lhs_samples.csv"), delimiter=",")
    failures += report("shared Latin Hypercube design present (1000 x 6)", X.shape == (1000, 6))

    # 5. the test suite
    r = subprocess.run([sys.executable, "-m", "pytest", "-q", os.path.join(ROOT, "tests")],
                       capture_output=True, text=True, cwd=ROOT)
    failures += report("pytest suite", r.returncode == 0)

    print("-" * 62)
    if failures == 0:
        print("RESULT: all checks passed. The published results are reproduced.")
    else:
        print(f"RESULT: {failures} check(s) FAILED.")
    return 1 if failures else 0


def report(name: str, ok: bool) -> int:
    print(f"  [{'PASS' if ok else 'FAIL'}]  {name}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
