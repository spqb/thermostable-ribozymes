# Inconsistencies found while cleaning the repository

Things that do not match between the paper, the supplementary material, the original files and the new scripts.
Each item says what was found, why it matters, and what the repository currently does.

## Results that disagree with the data

### 1. SI Fig. 4: threshold line vs caption

- **Found:** the caption of SI Fig. 4 says the dashed line is the 10th percentile of the overlap-group activity. The
  line drawn (and used to define hits everywhere) is at −1, i.e. 10% of the mean overlap activity. The two agree at
  37 °C (10th percentile = −1.02) but not at 60 °C, where the 10th percentile of the detected overlap sequences is
  −1.53 (−1.74 counting undetected ones as the lowest values).
- **Also:** the 99.9th percentile of the noise group is −1.05 at 37 °C, barely below the −1 threshold, so the
  statement that the threshold "lies above the 99.9th percentile" holds only marginally at 37 °C.
- **Repository:** `scripts/05_figures_experiment.py` draws the line at −1, as the paper does.

### 2. SI Table III: activities of the starting sequences

- **Found:** the table reports Lambert activities at 37 °C of −1.57, −1.23, −1.87, … (from the old `WT_summary.csv`).
  The Lambert dataset in `data/lambert_dataset.csv` (the latest version, `Matching_Data_corrected.csv`) gives
  −1.66, −1.23, −1.87, −1.58, −1.94, −1.51, −2.13, −0.96, −1.72, −1.87. Only start 1 differs by more than rounding
  (−1.57 vs −1.66).
- **Repository:** `results/tables/table_starting_sequences.csv` is built from `data/lambert_dataset.csv`.

## Design procedure

### 3. Parameters of the structure-only group and of the other regions

- **Found:** only β, λ_S (from the file names: β = 1.5, λ_S = 2.0 for ↓H ↓S; β = 1.25, λ_S = −3.0 for ↓H ↑S;
  β = 0, λ_S = 3.5 for ↑H ↓S) and a single λ_D = 2.0 (from the old notebook) are recorded.
  With λ_D = 2, 4 and 6 the structure-only group (↑H ↓S) yields no design within the budget (5–10 substitutions
  from the start, `S_sb < −33`, `H_gen > 110`) in the test runs; the λ_D used for the experiment is not known.
- **Repository:** `thermodca/config.py` uses λ_D = 2.0 for all three regions; see the note in the main README.

### 4. Distance from Azoarcus

- **Found:** the paper (Methods) restricts the designs to 25–30 substitutions from Azoarcus, but SI Table IV shows
  designs up to 35 (the ↑H ↓S group of starts 6, 7, 8 and 10, and the ↓H ↑S group of start 6).
- **Repository:** `AZO_RANGE = (25, 30)` in `thermodca/config.py`.

### 5. Gap symbol

- **Found:** the paper says the local field of the gap symbol was set to a very negative value; the old design
  notebook did not do it (no gap appeared in the synthesised sequences).
- **Repository:** the new sampler suppresses the gap (`GAP_FIELD = -1e4`).

## Smaller points

- **Number of sequences in the alignment:** Sec. III.C of the paper says `M = 817`, the Methods say 816.
- **ΔG panels of SI Figs. 2–3:** the original code replaced `[`/`]` by `(`/`)` before evaluating the free energy,
  which turns the pseudoknot into nested stems and gives a different pairing pattern; the Table II values use the
  structure without the pseudoknot. The SI figures are reproduced with the original convention
  (`pseudoknot="merge"` in `thermodca/rna.py`).
- **Region names:** because `S_sb` changed sign (the original "SB score" was `−S_sb`), the original `lowDCA_highSB` is now
  `lowDCA_lowSB` (↓H ↓S), the original `lowDCA_lowSB` is now `lowDCA_highSB` (↓H ↑S) and `highDCA_highSB` is now
  `highDCA_lowSB` (↑H ↓S).
- **Figures not reproduced:** SI Fig. 5 and SI Fig. 7 (the latter needs the full set of Lambert et al. designs, not
  only the 3,429 sequences in `data/lambert_dataset.csv`).
