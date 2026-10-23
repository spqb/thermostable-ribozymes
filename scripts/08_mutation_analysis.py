"""Step 8 - analysis of the individual mutations (exploratory; not shown in the main text of the paper).

For the jointly optimised designs (low H_gen, low S_sb) of a starting sequence:
  (a) per-position mutation frequency among the top-25% designs (activity at 60C) of each group
  (b) effect of each mutation of the best design on H_gen and S_sb, in the context of the starting sequence
      and in the context of the best design (reverse mutations), including pairs of mutated base-paired or
      adjacent positions
  (c) epistasis of a target mutation: ddX(S) = [X(S) - X(S \\ target)] - [X(target) - X(start)] for every subset S
      of the mutations of the best design that contains the target (X = H_gen or S_sb)
For all starting sequences:
  (d) G-C / Watson-Crick content of the consensus pairs in the top designs of each group
  (e) Hamming distances between the starting sequences

Reads   data/*.csv, results/tables/experiment_scores.csv (step 1)
Writes  results/figures/extra/start<N>_mutation_frequency, _mutation_effects, _epistasis_<mutation>
        results/tables/mutation_base_pair_content.csv, results/tables/start_distance_matrix.csv

Mutations are written <ref><position><alt> with 1-based positions of the 197-nt sequence, e.g. C71U.
Negative changes of H_gen and S_sb are favourable.
"""
import argparse
import itertools
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from gpsopt import data, dca, rna
from gpsopt.analysis import mutation_frequency
from gpsopt.config import DESIGN_GROUPS, FIGURES, GROUP_COLORS, GROUP_LABELS, TABLES
from gpsopt.plotting import save_figure, use_paper_style

EXTRA = FIGURES / "extra"


def make_scorer(kind, params, structure):
    """Function mapping a list of sequences to H_gen ("dca"), S_sb at 37C ("sb") or Delta G at 60C ("fe")."""
    if kind == "dca":
        return lambda seqs: dca.dca_energy(seqs, params)
    if kind == "sb":
        return lambda seqs: rna.sb_potential(seqs, structure, 37.0, n_workers=1)
    if kind == "fe":
        return lambda seqs: rna.consensus_free_energy(seqs, structure, 60.0, "remove", n_workers=1)
    raise ValueError(kind)


def list_mutations(ref, mutant):
    """[(position (0-based), ref nt, alt nt)] differing between two sequences."""
    return [(i, a, b) for i, (a, b) in enumerate(zip(ref, mutant)) if a != b]


def label(m):
    return f"{m[1]}{m[0] + 1}{m[2]}"


def apply(seq, muts):
    s = list(seq)
    for pos, _, alt in muts:
        s[pos] = alt
    return "".join(s)


def reverse(muts):
    return [(p, b, a) for p, a, b in muts]


# ---------------------------------------------------------------------------------------------- (a)
def plot_mutation_frequency(df, start_id, start_seq, structure, threshold=0.5):
    use_paper_style(16)
    paired = np.zeros(len(structure), dtype=bool)
    for i, j in rna.structure_pairs(structure):
        paired[[i, j]] = True
    fig, axs = plt.subplots(3, 1, figsize=(12, 6), dpi=300, sharex=True)
    for ax, g in zip(axs, ["highDCA_lowSB", "lowDCA_highSB", "lowDCA_lowSB"]):
        sub = df[(df["start_id"] == start_id) & (df["group"] == g)]
        top = sub.nlargest(int(0.25 * len(sub)), "norm_60")
        freq = mutation_frequency(top["sequence"].tolist(), start_seq)
        for pos in np.flatnonzero(paired):
            ax.axvspan(pos - 0.5, pos + 0.5, color="grey", alpha=0.2, zorder=0, linewidth=0)
        ax.bar(range(len(freq)), freq, color=GROUP_COLORS[g])
        above = np.flatnonzero(freq > threshold)
        ax.bar(above, freq[above], color="red")
        for pos in above:
            ax.text(pos, freq[pos], str(pos), ha="center", va="bottom", fontsize=14, color="red")
        ax.set_ylabel("Frequency")
        ax.set_title(GROUP_LABELS[g])
        ax.set_ylim(0, 1.2)
    axs[-1].set_xlabel("Position in Sequence")
    fig.tight_layout()
    save_figure(fig, f"extra/start{start_id}_mutation_frequency")
    plt.close(fig)


# ---------------------------------------------------------------------------------------------- (b)
def context_effects(context, muts, scorers, pairs):
    """Change of each score upon each mutation (and adjacent / base-paired pairs) in ``context``."""
    partner = {}
    for i, j in pairs:
        partner[i], partner[j] = j, i
    by_pos = {m[0]: m for m in muts}
    items = {label(m): [m] for m in muts}
    for m in muts:
        if m[0] + 1 in by_pos:
            items[f"{label(m)}+{label(by_pos[m[0] + 1])}"] = [m, by_pos[m[0] + 1]]
        c = partner.get(m[0])
        if c in by_pos and c > m[0]:
            items[f"{label(m)}+{label(by_pos[c])}"] = [m, by_pos[c]]
    seqs = [apply(context, v) for v in items.values()]
    return {k: scorers[k](seqs) - scorers[k]([context])[0] for k in scorers}, list(items)


def plot_mutation_effects(start_seq, best_seq, start_id, scorers, structure):
    use_paper_style(16)
    muts = list_mutations(start_seq, best_seq)
    pairs = rna.structure_pairs(structure)
    fig, axs = plt.subplots(1, 2, figsize=(12, 6), dpi=300)
    for ax, (context, ms, title, color) in zip(axs, [(start_seq, muts, "Start context", "cornflowerblue"),
                                                    (best_seq, reverse(muts), "Best design context", "coral")]):
        eff, names = context_effects(context, ms, scorers, pairs)
        x, y = eff["dca"], eff["sb"]
        xl, yl = (min(x.min() - 1, -1), max(x.max() + 1, 1)), (min(y.min() - 1, -1), max(y.max() + 1, 1))
        ax.set_xlim(*xl), ax.set_ylim(*yl)
        ax.fill_between([xl[0], 0], yl[0], 0, color="green", alpha=0.2)   # both scores decrease: favourable
        ax.fill_between([0, xl[1]], 0, yl[1], color="red", alpha=0.2)     # both increase: unfavourable
        ax.fill_between([0, xl[1]], yl[0], 0, color="yellow", alpha=0.2)
        ax.fill_between([xl[0], 0], 0, yl[1], color="yellow", alpha=0.2)
        ax.scatter(x, y, color=color, edgecolor="black", s=100)
        for xi, yi, n in zip(x, y, names):
            ax.annotate(n, (xi + 0.3, yi), fontsize=12, alpha=0.7)
        ax.axhline(0, color="black", linewidth=0.5), ax.axvline(0, color="black", linewidth=0.5)
        ax.grid(lw=0.5, ls="--", alpha=0.7)
        ax.set_xlabel(r"$\Delta \mathcal{H}_{\mathrm{gen}}$"), ax.set_ylabel(r"$\Delta \mathcal{S}_{\mathrm{sb}}$")
        ax.set_title(title)
    fig.tight_layout()
    save_figure(fig, f"extra/start{start_id}_mutation_effects")
    plt.close(fig)


# ---------------------------------------------------------------------------------------------- (c)
def epistasis_scan(start_seq, muts, target, scorer):
    """ddX for every subset of ``muts`` containing ``target``; returns (n_other_mutations, ddX, subsets)."""
    others = [m for m in muts if m != target]
    subsets = [(target,) + c for k in range(len(others) + 1) for c in itertools.combinations(others, k)]
    full = [apply(start_seq, s) for s in subsets]
    without = [apply(start_seq, s[1:]) for s in subsets]
    base = scorer([start_seq, apply(start_seq, [target])])
    dd = scorer(full) - scorer(without) - (base[1] - base[0])
    return np.array([len(s) - 1 for s in subsets]), dd, subsets


def plot_epistasis(start_seq, best_seq, start_id, target_label, scorers):
    use_paper_style(16)
    muts = list_mutations(start_seq, best_seq)
    target = next((m for m in muts if label(m) == target_label), None)
    if target is None:
        print(f"mutation {target_label} not in the best design of start {start_id}: {[label(m) for m in muts]}")
        return
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.5), dpi=300)
    for ax, kind, ylabel, color in [(axs[0], "dca", r"$\Delta\Delta \mathcal{H}_{\mathrm{gen}}$", "cornflowerblue"),
                                    (axs[1], "sb", r"$\Delta\Delta \mathcal{S}_{\mathrm{sb}}$", "coral")]:
        n_other, dd, subsets = epistasis_scan(start_seq, muts, target, scorers[kind])
        ax.scatter(n_other, dd, alpha=0.7, c=color)
        best = [np.argmin(np.where(n_other == k, dd, np.inf)) for k in range(1, n_other.max() + 1)]
        ax.scatter(n_other[best], dd[best], color="red", s=80, zorder=5, edgecolor="black")
        previous = {target}
        for k, b in zip(range(1, n_other.max() + 1), best):
            added = set(subsets[b]) - previous
            previous = set(subsets[b])
            ax.annotate("+" + ",".join(label(m) for m in added), (k, dd[b]), xytext=(5, 5),
                        textcoords="offset points", fontsize=11, color="red")
        ax.set_xlabel("Number of other mutations")
        ax.set_ylabel(ylabel)
        ax.grid(lw=0.5, ls="--", alpha=0.7)
    fig.suptitle(f"Epistasis of {target_label} (start {start_id})")
    fig.tight_layout()
    save_figure(fig, f"extra/start{start_id}_epistasis_{target_label}")
    plt.close(fig)


# ---------------------------------------------------------------------------------------------- (d), (e)
def base_pair_content(df, starts, structure, n_top=50):
    rows = []
    for sid, start in starts["sequence"].items():
        wt = rna.base_pair_counts(start, structure)
        for g in DESIGN_GROUPS:
            sub = df[(df["start_id"] == sid) & (df["group"] == g)].nlargest(n_top, "norm_60")
            counts = pd.DataFrame([rna.base_pair_counts(s, structure) for s in sub["sequence"]])
            rows.append({"start_id": sid, "group": g, "n_top": len(sub), "start_gc": wt["gc"], "start_wc": wt["wc"],
                         "mean_gc": counts["gc"].mean(), "std_gc": counts["gc"].std(ddof=0),
                         "mean_wc": counts["wc"].mean(), "std_wc": counts["wc"].std(ddof=0)})
    return pd.DataFrame(rows).round(2)


def main(args):
    EXTRA.mkdir(parents=True, exist_ok=True)
    azo = data.load_azoarcus()
    structure = azo["dot_bracket"]
    starts = data.load_starting_sequences()
    df = data.load_experiment_scored()
    params = dca.load_model()
    scorers = {k: make_scorer(k, params, structure) for k in ("dca", "sb")}

    sid = args.start
    start_seq = starts.loc[sid, "sequence"]
    joint = df[(df["start_id"] == sid) & (df["group"] == "lowDCA_lowSB")]
    best = joint.nlargest(1, "norm_60").iloc[0]
    print(f"start {sid}: best design {best['id']} (normalised activity at 60C: {best['norm_60']:.2f}), mutations:",
          [label(m) for m in list_mutations(start_seq, best["sequence"])])

    plot_mutation_frequency(df, sid, start_seq, structure)
    plot_mutation_effects(start_seq, best["sequence"], sid, scorers, structure)
    for target in args.targets:
        plot_epistasis(start_seq, best["sequence"], sid, target, scorers)

    base_pair_content(df, starts, structure).to_csv(TABLES / "mutation_base_pair_content.csv", index=False)
    dist = np.array([[int(sum(a != b for a, b in zip(s, t))) for t in starts["sequence"]] for s in starts["sequence"]])
    pd.DataFrame(dist, index=starts.index, columns=starts.index).to_csv(TABLES / "start_distance_matrix.csv")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--start", type=int, default=9, help="starting sequence analysed in (a)-(c)")
    p.add_argument("--targets", nargs="*", default=["C71U"], help="target mutations for the epistasis scan (c)")
    main(p.parse_args())
