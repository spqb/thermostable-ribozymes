# Parameters used to generate the designs

The designs were sampled from

`H_eff(x | x_start) = β H_gen(x) + λ_S S_sb(x) + λ_D d_H(x_start, x)`

(Eq. 6 of the paper; `S_sb = -Σ p_ij`, so low `S_sb` is favourable). In the code, β rescales the whole DCA model,
λ_S is the weight of the Metropolis test on `S_sb` and λ_D is the weight that keeps the sequences close to the starting
sequence (an extra field `+λ_D` on the nucleotide of the starting sequence at each position).

## Values per region

| Region (paper) | Region (original file names) | β | λ_S | λ_D |
|---|---|---|---|---|
| ↓H_gen ↓S_sb (`lowDCA_lowSB`, joint) | `lowDCA_highSB_beta1.5_sb2.0` | 1.5 | +2.0 | not recorded |
| ↓H_gen ↑S_sb (`lowDCA_highSB`, evolutionary score only) | `lowDCA_lowSB_beta1.25_sb-3.00` | 1.25 | −3.0 | not recorded |
| ↑H_gen ↓S_sb (`highDCA_lowSB`, structure score only) | `highDCA_highSB_beta0_sb3.5` | 0 | +3.5 | not recorded |

Where the values come from:

- β and λ_S are read from the names of the original sequence files (`beta<β>_sb<λ_S>`). The "sb" value is the
  weight of the structure term in the energy `λ_S · S_sb`, so a positive value favours low `S_sb` (more consensus
  contacts) and a negative one favours high `S_sb`. This is consistent with the regions they populate (Table I).
- β = 0 for the structure-only group means that `H_gen` was not used in the sampling (only the structure term and the
  distance term); the sequences are then filtered on `H_gen > 110`.
- **λ_D is not recorded.** The only value found in the original notebook is λ_D = 2.0 (variable `is_power`), used for
  the ↓H_gen ↑S_sb region. `thermodca/config.py` assumes λ_D = 2.0 for the three regions. In test runs with λ_D = 2, 4
  and 6 the structure-only group did not yield designs within the budget (see below), so its original λ_D was
  probably different (or the sampling lengths were).

## Other settings

| | |
|---|---|
| Chains per starting sequence and region | 30,000 |
| "Slow" steps (Metropolis test on `S_sb`) | 20 |
| "Fast" sweeps between two slow steps | 10 (one sweep = L attempted Gibbs updates) |
| `S_sb` temperature | 37 °C |
| Substitutions from the starting sequence | 5–10 |
| Substitutions from Azoarcus | 25–30 (some designs reach 35, see SI Table IV) |
| Gap state | never proposed |

Selection of the designs (Table I of the paper), among the sampled sequences within the distance budgets:

| Region | `H_gen` | `S_sb` | Designs per starting sequence |
|---|---|---|---|
| ↓H_gen ↓S_sb | < `H_gen(x_start)` | < −33 | 150 |
| ↓H_gen ↑S_sb | < `H_gen(x_start)` | > −10 | 50 |
| ↑H_gen ↓S_sb | > 110 | < −33 | 50 |

## Reproducing the sampling

`scripts/02_design_mutants.py` uses these values (`DESIGN_PARAMS` in `thermodca/config.py`). Sampling is stochastic and no
random seed of the original runs was recorded, so the synthesised sequences themselves are in `data/experiment.csv`
and new runs give different designs.

A test with 600 chains for starting sequence 5 selected 150 designs for the joint region and 27 for the
evolutionary-only region, but none for the structure-only region. With 200 chains and λ_D = 4 the chains stayed
within 5–10 substitutions of the start on average but none satisfied `S_sb < −33` and `H_gen > 110` together; with
λ_D = 6 they barely moved.
