#!/usr/bin/env bash
# Reproduces all tables and figures from the files in data/ (step 2, the stochastic design, is optional).
set -euo pipefail
cd "$(dirname "$0")"
uv run python scripts/01_score_sequences.py
uv run python scripts/04_figures_lambert_data.py
uv run python scripts/05_figures_experiment.py
uv run python scripts/06_tables.py
uv run python scripts/07_md_candidates.py
uv run python scripts/08_mutation_analysis.py
uv run python scripts/09_mnist_toy.py
# uv run python scripts/02_design_mutants.py --n-seqs 30000   # new designs (hours); results in results/design/
