# SGR gang-dynamics repository (MBE-8616)
#
#   make verify     reproduce the published results (Python)
#   make test       run the Python test suite
#   make figures    regenerate every figure by executing the notebooks
#   make matlab     run the MATLAB/Octave validation suite
#   make all        verify + figures

PYTHON  ?= python3
OCTAVE  ?= octave
NOTEBOOKS = notebooks/01_model_and_equilibria.ipynb \
            notebooks/02_bifurcation_diagrams.ipynb \
            notebooks/03_prcc_sensitivity.ipynb

.PHONY: all verify test figures matlab clean

all: verify figures

verify:
	$(PYTHON) verify.py

test:
	$(PYTHON) -m pytest -q tests

# Executing the notebooks IS how the figures are produced: notebook 01 also
# regenerates python/sgr_engine.py, so the engine can never drift from its source.
figures:
	$(PYTHON) -m jupyter nbconvert --to notebook --execute --inplace \
	    --ExecutePreprocessor.timeout=1200 $(NOTEBOOKS)

matlab:
	$(OCTAVE) --no-gui -q --eval "addpath('matlab'); n = sgrValidate(); exit(n > 0);"

clean:
	rm -rf python/__pycache__ tests/__pycache__ .pytest_cache
