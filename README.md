# Biophysics-informed generative models for thermostable ribozymes

![Graphical abstract](docs/graphical-abstract-ribozymes.png)

Code and data of *"Biophysics-informed generative models improve activity and thermostability of designed
ribozymes"*. Starting from Azoarcus group I intron variants that self-splice at 37 °C but not at 60 °C, we sample
variants with 5–10 additional substitutions from an effective Hamiltonian that combines

- a DCA (Potts) model of the intron family, `H_gen`, which keeps the sequences functionally plausible, and
- a structure-based potential from the Turner nearest-neighbour model (ViennaRNA),
  `S_sb = -Σ_(i,j)∈consensus p_ij`, the (negative) expected number of consensus base pairs at equilibrium,

`H_eff(x | x_start) = β H_gen(x) + λ_S S_sb(x) + λ_D d_H(x_start, x)`.

Both terms are needed to recover activity at 60 °C (hit rate 72.5% for the jointly optimised designs versus 7.8% and
0.4% for either term alone).

## Repository layout

```
data/        input data (see data/README.md)
thermodca/      small python package: DCA model, Turner-model scores, design sampler, plotting style
scripts/     numbered scripts that reproduce the pipeline (below)
results/     tables, figures and design outputs written by the scripts
run_all.sh   runs the whole analysis
```

## Installation

The environment is managed with [uv](https://docs.astral.sh/uv/):

```bash
uv sync                          # creates .venv from pyproject.toml / uv.lock
uv run python scripts/01_score_sequences.py
```

adabmDCA 2.0 is installed from the `develop` branch of <https://github.com/spqb/adabmDCApy> (see
`[tool.uv.sources]` in `pyproject.toml`; switch to `main` once it is released there). Figures use LaTeX
(`latex`, `dvipng` and the Computer Modern fonts must be installed on the system). A GPU is optional (it only
makes step 2 fast).

## Pipeline

| script | content | paper |
|---|---|---|
| `01_score_sequences.py` | `H_gen`, `S_sb` (37 °C, 60 °C), ΔG of the consensus structure for the Lambert et al. data and for all tested sequences | – |
| `02_design_mutants.py` | samples new designs from `H_eff` (multi-step Monte Carlo: Gibbs sweeps for the fast terms, Metropolis on `S_sb`) and selects them by region (Table I) | Sec. III.E, IV.C |
| `04_figures_lambert_data.py` | `S_sb` vs `H_gen` and activity at 60 °C; activity vs distance from Azoarcus | Fig. 2, SI Figs. 2–3, SI Fig. 7 |
| `05_figures_experiment.py` | phase space of start 5, activity violins, Δ activity, overlap/noise controls | Figs. 3–5, SI Figs. 4, 6 |
| `06_tables.py` | hit rates, activity quartiles, ΔΔG, ΔH-bonds, pool statistics | Table II, SI Tables II–IV |
| `07_md_candidates.py` | best/worst design per start (for MD), FORNA input files | Appendix C |
| `08_mutation_analysis.py` | per-mutation effects, epistasis, base-pair content (exploratory) | not in the paper |
| `09_mnist_toy.py` | MNIST 0→8 toy example | Fig. 1 |

Run everything with `./run_all.sh` (it calls `uv run`) or run a single script with `uv run python scripts/<name>.py`; each script has a docstring describing its inputs and outputs. Step 3 does not
exist: the high-throughput measurements are already in `data/experiment.csv`.

## Conventions

- **Start sequences** are the ten variants of Lambert et al. that we mutate (`start_id` 1–10). Only the Azoarcus
  sequence is a "wild type".
- **`S_sb` is negative**: `S_sb = -Σ p_ij`; low (more negative) is favourable.
- **Groups** follow the paper: `lowDCA_lowSB` = ↓H_gen ↓S_sb (joint), `lowDCA_highSB` = ↓H_gen ↑S_sb (evolutionary
  score only), `highDCA_lowSB` = ↑H_gen ↓S_sb (structure score only), plus the `noise` and `overlap` controls.
- **Activity** is `log10(f_sel / f_ref)`. It is normalised by subtracting, for each temperature, the mean over the
  detected overlap-group sequences; a *hit* has normalised activity ≥ −1. Sequences with no detected activity are
  empty (`NaN`) in `data/experiment.csv` and in the Lambert et al. table.
- Sequences have 197 nt: `GUG` + 193-nt aligned intron + `G`. The DCA model acts on the 193 central nucleotides and
  uses the adabmDCA 2.0 alphabet `-ACGU` (gap first).

The sampling parameters of the designs are in [docs/design_parameters.md](docs/design_parameters.md).

## Not reproduced here

- `02_design_mutants.py` is tested but not tuned for every group. With 600 chains it selects designs for the joint
  (`lowDCA_lowSB`, 150/600) and evolutionary-only (`lowDCA_highSB`, 27/600) groups, but none for the structure-only
  group (`highDCA_lowSB`): the sampling parameters of that group (λ_D, number of steps) used in the experiment are not
  recorded, and λ_D = 2, 4, 6 did not give designs within the budget (5–10 substitutions, `S_sb < -33`, `H_gen > 110`).

The following supplementary material has no code in this repository yet:

- SI Fig. 5 (MNIST classifier-only control);
- the structure drawings were made with FORNA / ChimeraX; step 7 only writes the FORNA
  input files;
- the sampling of the MNIST examples of Fig. 1 (only the sampled configurations are stored in
  `data/mnist_toy_samples.npz`, the RBM is not included);
- structure prediction (AlphaFold 3) and molecular dynamics of Appendix C (only the choice of the mutants is
  implemented, step 7).

## Reference

Lambert et al., *Exploring the space of self-reproducing ribozymes using generative models*,
Nat. Commun. 16, 7836 (2025).
