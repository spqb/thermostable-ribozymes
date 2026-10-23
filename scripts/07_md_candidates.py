"""Step 7 - choose the designs used for the molecular-dynamics study and prepare FORNA inputs.

For each starting sequence we pick, among the designs with a hit at 37C and detected activity at 60C,
  good mutant: jointly optimised design (low H_gen, low S_sb) with the largest  delta = activity_60 - activity_37
  bad mutant : DCA-only design        (low H_gen, high S_sb) with the smallest delta
Rank 1 is the extreme design, rank 2 the second most extreme (used when the first is not suitable).

Reads   data/experiment.csv, data/lambert_dataset.csv, data/azoarcus.csv, results/tables/experiment_scores.csv
Writes  results/md_candidates/md_candidates.csv        chosen designs (ranks 1 and 2)
        results/md_candidates/best_designs.csv          design with the highest activity at 60C for each start
                                                        (jointly optimised group)
        results/md_candidates/forna/<name>.txt          FORNA input files (sequence, structure, colours)

FORNA (http://rna.tbi.univie.ac.at/forna/): use the lines under "# Sequence and structure" as the molecule and
the lines under "# Color codes" as the colour scheme.
Colours of the mutant files (relative to the starting sequence and to Azoarcus):
  forestgreen    position identical to Azoarcus that was mutated
  coral          mutation back to the Azoarcus nucleotide
  cornflowerblue mutation to a nucleotide that is neither the start nor the Azoarcus one
  lightgrey      position different from Azoarcus, left unchanged
The "<start>_mutation_frequency" files colour the starting sequence (white -> red) by the frequency of
mutation among the top-25% designs (activity at 60C) of the jointly optimised group.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd

from thermodca import data
from thermodca.analysis import INACTIVE_DELTA, delta_activity_table, mutation_frequency
from thermodca.config import MD

COLORS = dict(away_from_azo="forestgreen", back_to_azo="coral", new="cornflowerblue", unchanged="lightgrey")


def mutant_colors(azo, start, mutant):
    """{1-based position: colour} for the FORNA colour scheme of a mutant."""
    out = {}
    for pos, (a, s, m) in enumerate(zip(azo, start, mutant), start=1):
        if a == s:
            if a != m:
                out[pos] = COLORS["away_from_azo"]
        elif a == m:
            out[pos] = COLORS["back_to_azo"]
        elif s == m:
            out[pos] = COLORS["unchanged"]
        else:
            out[pos] = COLORS["new"]
    return out


def write_forna(path, name, sequence, structure, colors):
    lines = ["# Sequence and structure", f">{name}", sequence, structure, "", "# Color codes", f">{name}"]
    lines += colors if isinstance(colors, list) else [f"{pos}:{c}" for pos, c in colors.items()]
    Path(path).write_text("\n".join(lines) + "\n")


def main():
    MD.mkdir(parents=True, exist_ok=True)
    (MD / "forna").mkdir(exist_ok=True)
    azo = data.load_azoarcus()
    starts = data.load_starting_sequences()
    df = data.load_experiment_scored()
    hits = delta_activity_table(df)
    hits = hits[hits["delta_activity"] > INACTIVE_DELTA]  # detected activity at 60C

    rows = []
    for sid in starts.index:
        sub = hits[hits["start_id"] == sid]
        good = sub[sub["group"] == "lowDCA_lowSB"].sort_values("delta_activity", ascending=False)
        bad = sub[sub["group"] == "lowDCA_highSB"].sort_values("delta_activity", ascending=True)
        for kind, table in [("good", good), ("bad", bad)]:
            for rank, (_, row) in enumerate(table.head(2).iterrows(), start=1):
                rows.append(dict(start_id=sid, kind=kind, rank=rank, id=row["id"], group=row["group"],
                                 sequence=row["sequence"], norm_37=row["norm_37"], norm_60=row["norm_60"],
                                 delta_activity=row["delta_activity"], sb_37=row["sb_37"],
                                 dca_energy=row["dca_energy"], dist_azo=row["dist_azo"]))
                if rank == 1:
                    write_forna(MD / "forna" / f"start_{sid}_{kind}_mutant.txt", f"start_{sid}_{kind}_mutant",
                                row["sequence"], azo["dot_bracket"],
                                mutant_colors(azo["sequence"], starts.loc[sid, "sequence"], row["sequence"]))
        # mutation frequency among the top-25% (activity at 60C) designs of the joint group
        joint = df[(df["start_id"] == sid) & (df["group"] == "lowDCA_lowSB")]
        top = joint.nlargest(max(1, len(joint) // 4), "norm_60")
        freq = mutation_frequency(top["sequence"].tolist(), starts.loc[sid, "sequence"])
        write_forna(MD / "forna" / f"start_{sid}_mutation_frequency.txt", f"start_{sid}", starts.loc[sid, "sequence"],
                    azo["dot_bracket"], ["range=white:red", " ".join(f"{f:.4f}" for f in freq)])
    best = (df[df["group"] == "lowDCA_lowSB"].sort_values("norm_60", ascending=False)
            .groupby("start_id").head(1).sort_values("start_id"))
    best["start_id"] = best["start_id"].astype(int)
    best[["start_id", "id", "sequence", "norm_37", "norm_60", "dca_energy", "sb_37", "dist_start"]].to_csv(
        MD / "best_designs.csv", index=False)
    out = pd.DataFrame(rows)
    out.to_csv(MD / "md_candidates.csv", index=False)
    print(out[out["rank"] == 1][["start_id", "kind", "id", "delta_activity"]].to_string(index=False))


if __name__ == "__main__":
    main()
