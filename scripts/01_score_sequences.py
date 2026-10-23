"""Step 1 - score all sequences with the DCA model and the Turner-model potentials.

Reads   data/lambert_dataset.csv, data/experiment.csv, data/azoarcus.csv, data/dca_model_lambert.dat
Writes  results/tables/lambert_scores.csv      H_gen, S_sb(37C, 60C), Delta G of the consensus structure
        results/tables/experiment_scores.csv   H_gen, S_sb(37C), distances from Azoarcus / starting sequence

H_gen: DCA energy (adabmDCA). S_sb = -sum of the equilibrium pairing probabilities of the consensus contacts
(ViennaRNA partition function). Delta G: Turner free energy of the sequence in the consensus structure.
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd

from thermodca import data, dca, rna
from thermodca.config import TABLES


def main(n_workers):
    TABLES.mkdir(parents=True, exist_ok=True)
    azo = data.load_azoarcus()
    structure = azo["dot_bracket"]
    params = dca.load_model()

    # ---- Lambert et al. dataset
    lam = data.load_lambert()
    seqs = lam["sequence"].tolist()
    out = lam[["name", "model", "start_id"]].copy()
    out["dca_energy"] = dca.dca_energy(seqs, params)
    for T in (37, 60):
        out[f"sb_{T}"] = rna.sb_potential(seqs, structure, T, n_workers)
        # "merge" reproduces the Delta G values shown in the supplementary figures (see rna.py)
        out[f"dg_{T}"] = rna.consensus_free_energy(seqs, structure, T, "merge", n_workers)
    out.to_csv(TABLES / "lambert_scores.csv", index=False)
    print("Lambert dataset:", out.shape)

    # ---- designed sequences and controls
    exp = data.load_experiment()
    seqs = exp["sequence"].tolist()
    out = exp[["id", "group", "start_id"]].copy()
    out["dca_energy"] = dca.dca_energy(seqs, params)
    out["sb_37"] = rna.sb_potential(seqs, structure, 37.0, n_workers)
    out["dist_azo"] = data.hamming(seqs, azo["sequence"])
    starts = data.load_starting_sequences()["sequence"]
    out["dist_start"] = pd.array([pd.NA] * len(out), dtype="Int64")
    for sid, ref in starts.items():
        mask = (exp["start_id"] == sid).to_numpy()
        out.loc[mask, "dist_start"] = data.hamming(exp.loc[mask, "sequence"].tolist(), ref)
    out.to_csv(TABLES / "experiment_scores.csv", index=False)
    print("Experiment:", out.shape)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--n-workers", type=int, default=None, help="processes for ViennaRNA (default: all cores)")
    main(p.parse_args().n_workers)
