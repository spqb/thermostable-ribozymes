"""Step 2 - design variants of the starting sequences with the effective Hamiltonian (Eq. 6).

Three regions of the (H_gen, S_sb) plane are populated for each starting sequence (see thermodca/config.py):
  lowDCA_lowSB  (low H_gen, low S_sb)   jointly optimised designs
  lowDCA_highSB (low H_gen, high S_sb)  evolutionary score only
  highDCA_lowSB (high H_gen, low S_sb)  structure-based score only
Sampling is stochastic: running this script produces new designs, it does not regenerate the sequences that
were synthesised (those are in data/experiment.csv).

Reads   data/{lambert_dataset,azoarcus}.csv, data/dca_model_lambert.dat
Writes  results/design/selected_designs.csv        designs satisfying the region/budget constraints
        results/design/candidates_start<N>_<group>.csv.gz   all sampled chains (with --save-candidates)

Cost: every slow step scores n_seqs sequences with ViennaRNA, ~n_slow * n_seqs partition functions per run.
Use e.g. --n-seqs 500 for a quick test.
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import torch

from thermodca import data, dca, design
from thermodca.config import DESIGN, DESIGN_GROUPS, DESIGN_PARAMS, N_PER_REGION


def main(args):
    DESIGN.mkdir(parents=True, exist_ok=True)
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    azo = data.load_azoarcus()
    starts = data.load_starting_sequences()
    params = dca.load_model()
    selected = []
    for sid in args.starts:
        start = starts.loc[sid, "sequence"]
        start_energy = float(dca.dca_energy([start], params)[0])
        for group in args.groups:
            p = DESIGN_PARAMS[group]
            print(f"start {sid}, {group}: {p}")
            seqs = design.sample_designs(start, params, azo["dot_bracket"], p["beta"], p["lambda_S"], p["lambda_D"],
                                         args.n_seqs, args.n_slow, args.n_fast, args.n_workers)
            cand = design.score_candidates(seqs, params, azo["dot_bracket"], start, azo["sequence"], args.n_workers)
            if args.save_candidates:
                cand.to_csv(DESIGN / f"candidates_start{sid}_{group}.csv.gz", index=False)
            sel = design.select_designs(cand, group, start_energy, N_PER_REGION[group], args.seed)
            print(f"   {len(sel)} designs selected out of {len(cand)} chains")
            sel.insert(0, "group", group)
            sel.insert(1, "start_id", sid)
            selected.append(sel)
    out = pd.concat(selected, ignore_index=True)
    out.insert(0, "id", [f"start{s}_{g}_{i:03d}" for i, (s, g) in enumerate(zip(out.start_id, out.group), 1)])
    out.to_csv(DESIGN / "selected_designs.csv", index=False)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--starts", type=int, nargs="+", default=list(range(1, 11)))
    p.add_argument("--groups", nargs="+", default=DESIGN_GROUPS, choices=DESIGN_GROUPS)
    p.add_argument("--n-seqs", type=int, default=30000, help="independent chains per start and group")
    p.add_argument("--n-slow", type=int, default=20, help="number of Metropolis steps on S_sb")
    p.add_argument("--n-fast", type=int, default=10, help="Gibbs sweeps between two Metropolis steps")
    p.add_argument("--n-workers", type=int, default=None)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--save-candidates", action="store_true")
    main(p.parse_args())
