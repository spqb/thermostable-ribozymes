"""Step 6 - summary tables of the experiment.

Reads   data/experiment.csv, data/lambert_dataset.csv, results/tables/*_scores.csv (step 1)
Writes  results/tables/table_pooled.csv(.tex)      Table II  : results pooled over the ten starting sequences
        results/tables/table_per_start.csv(.tex)   SI Table II: results for each starting sequence and group
        results/tables/table_starting_sequences.csv SI Table III: starting sequences (Lambert activities)
        results/tables/table_pools.csv             SI Table IV: distance statistics of the design pools

Definitions
  hit rate   fraction of designs with overlap-normalised activity >= -1
  top 25%    75th percentile of the normalised activity (undetected sequences count as -5)
  ddG        mean change of the Turner free energy (60C) of the design in the consensus structure
             without the pseudoknot, relative to its starting sequence
  dH-bonds   mean change of the number of hydrogen bonds (GC: 3, AU/GU: 2) of the consensus pairs
             without the pseudoknot formed by the design, relative to its starting sequence
"""
import itertools
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

from thermodca import data, rna
from thermodca.config import DESIGN_GROUPS, TABLES

DDG_TEMPERATURE = 60.0


UNDETECTED = -5.0  # value given to sequences without detected activity when computing percentiles


def top_quartile(x):
    return np.percentile(x.fillna(UNDETECTED).to_numpy(), 75)


def mean_pairwise_distance(seqs):
    arr = np.frombuffer("".join(seqs).encode(), dtype="S1").reshape(len(seqs), -1)
    n = len(seqs)
    total = sum((arr[i] != arr[i + 1 :]).sum() for i in range(n - 1))
    return total / (n * (n - 1) / 2)


def add_structure_columns(df, starts, structure):
    """ddG and dH-bonds of each design relative to its starting sequence."""
    ref_dg = {i: rna.consensus_free_energy([s], structure, DDG_TEMPERATURE, "remove")[0]
              for i, s in starts["sequence"].items()}
    ref_hb = {i: rna.hydrogen_bonds(s, structure, nested_only=True) for i, s in starts["sequence"].items()}
    df = df.copy()
    df["ddg"] = rna.consensus_free_energy(df["sequence"].tolist(), structure, DDG_TEMPERATURE, "remove") \
        - df["start_id"].map(ref_dg)
    df["dhb"] = [rna.hydrogen_bonds(s, structure, nested_only=True) - ref_hb[i]
                 for s, i in zip(df["sequence"], df["start_id"])]
    return df


def summarize(g):
    return pd.Series({
        "N": len(g),
        "top25_act_37": top_quartile(g["norm_37"]),
        "top25_act_60": top_quartile(g["norm_60"]),
        "hit_rate_37": 100 * g["hit_37"].mean(),
        "hit_rate_60": 100 * g["hit_60"].mean(),
        "ddG_kcal_mol": g["ddg"].mean(),
        "dH_bonds": g["dhb"].mean(),
    })


GROUP_TEX = {
    "lowDCA_lowSB": r"$\downarrow \mathcal{H}_\mathrm{gen}\ \downarrow \mathcal{S}_\mathrm{sb}$",
    "lowDCA_highSB": r"$\downarrow \mathcal{H}_\mathrm{gen}\ \uparrow \mathcal{S}_\mathrm{sb}$",
    "highDCA_lowSB": r"$\uparrow \mathcal{H}_\mathrm{gen}\ \downarrow \mathcal{S}_\mathrm{sb}$",
}


def tex_row(label, r, start=None):
    num = lambda x: f"{x + 0.0:.2f}".replace("-0.00", "0.00").replace("-", "$-$")
    cells = [label, str(int(r["N"])), num(r["top25_act_37"]), num(r["top25_act_60"]), f"{r['hit_rate_37']:.1f}\\%",
             f"{r['hit_rate_60']:.1f}\\%", num(r["ddG_kcal_mol"]), num(r["dH_bonds"])]
    if start is not None:
        cells = [start] + cells
    return " & ".join(cells) + r" \\"


def write_latex(table, per_start, path):
    """Booktabs tabular with the columns of Table II (``per_start``: one block per starting sequence)."""
    lead = "l" if per_start else ""
    sub = " & " if per_start else ""
    head = [r"\begin{tabular}{" + lead + "lrrrrrrr}", r"\toprule",
            ("Start & " if per_start else "") + r"Sampling condition & $N$ & \multicolumn{2}{c}{Top 25\% act.} & "
            r"\multicolumn{2}{c}{Hit rate} & $\Delta\Delta G$ & $\Delta$ H-bonds \\",
            sub + r" & & (37$^\circ$C) & (60$^\circ$C) & (37$^\circ$C) & (60$^\circ$C) & (kcal/mol) & (consensus) \\",
            r"\midrule"]
    body = []
    if per_start:
        last = table.index.get_level_values("start_id").max()
        for sid, block in table.groupby(level="start_id"):
            for k, ((_, g), r) in enumerate(block.iterrows()):
                body.append(tex_row(GROUP_TEX[g], r, start=f"st\\_{sid}" if k == 0 else ""))
            if sid != last:
                body.append(r"\midrule")
    else:
        body = [tex_row(GROUP_TEX[g], r) for g, r in table.iterrows()]
    Path(path).write_text("\n".join(head + body + [r"\bottomrule", r"\end{tabular}"]) + "\n")


def main():
    TABLES.mkdir(parents=True, exist_ok=True)
    structure = data.load_azoarcus()["dot_bracket"]
    starts = data.load_starting_sequences()
    df = data.load_experiment_scored()
    designs = df[df["group"].isin(DESIGN_GROUPS)].copy()
    designs["start_id"] = designs["start_id"].astype(int)
    designs = add_structure_columns(designs, starts, structure)

    pooled = designs.groupby("group").apply(summarize, include_groups=False).loc[DESIGN_GROUPS]
    pooled["N"] = pooled["N"].astype(int)
    pooled.round(2).to_csv(TABLES / "table_pooled.csv")
    per_start = designs.groupby(["start_id", "group"]).apply(summarize, include_groups=False)
    per_start = per_start.reindex(pd.MultiIndex.from_product([sorted(designs["start_id"].unique()), DESIGN_GROUPS],
                                                             names=["start_id", "group"]))
    per_start["N"] = per_start["N"].astype(int)
    per_start.round(2).to_csv(TABLES / "table_per_start.csv")
    write_latex(pooled.round(2), False, TABLES / "table_pooled.tex")
    write_latex(per_start.round(2), True, TABLES / "table_per_start.tex")
    print(pooled.round(2).to_string())

    # SI Table III: starting sequences
    lam_scores = pd.read_csv(TABLES / "lambert_scores.csv")
    t3 = starts[["name", "activity_37", "activity_60", "distance_azo"]].rename(columns={"name": "lambert_name"})
    t3["sb_37"] = lam_scores.set_index("name").loc[t3["lambert_name"], "sb_37"].to_numpy()
    t3["dca_energy"] = lam_scores.set_index("name").loc[t3["lambert_name"], "dca_energy"].to_numpy()
    t3.to_csv(TABLES / "table_starting_sequences.csv")

    # SI Table IV: pools
    rows = []
    for (sid, g), sub in designs.groupby(["start_id", "group"]):
        rows.append(dict(start_id=sid, group=g, N=len(sub),
                         mean_intra_distance=round(mean_pairwise_distance(sub["sequence"].tolist()), 1),
                         min_dist_azo=sub["dist_azo"].min(), max_dist_azo=sub["dist_azo"].max(),
                         min_dist_start=sub["dist_start"].min(), max_dist_start=sub["dist_start"].max()))
    pd.DataFrame(rows).to_csv(TABLES / "table_pools.csv", index=False)


if __name__ == "__main__":
    main()
