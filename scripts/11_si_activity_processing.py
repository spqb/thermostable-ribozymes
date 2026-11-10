"""Step 11 - SI figures and table of the section "Activity processing and controls".

Reads   data/experiment.csv (through thermodca.data.load_experiment_scored)
        data/lambert_dataset.csv (for the cross-experiment check of the reference group)
Writes  results/figures/si_fig06_noise_threshold   noise-control distributions and the hit threshold
        results/figures/si_fig07_reference_group   reference-group coverage, normalised activity and
                                                   agreement with the measurement of Lambert et al.
        results/tables/si_table04_read_statistics.csv   read counts, coverage and missing/undetected
                                                        fractions per group and temperature

Conventions (see the repository README): activity is the raw log10(f_sel/f_ref) of data/experiment.csv;
the reference-normalised activity subtracts, per temperature, the mean over the detected sequences of the
355-sequence overlap group; a hit has normalised activity >= -1. A sequence with no reference reads is
missing (no activity), one with reference but no selected reads has an undetected activity.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from thermodca import data
from thermodca.config import HIT_THRESHOLD, TABLES
from thermodca.plotting import save_figure, use_paper_style

TEMPS = (37, 60)
TCOLOR = {37: "cornflowerblue", 60: "firebrick"}


def panel_labels(axs, offset=(-0.14, 1.06)):
    """Label the axes a, b, c, ... in reading order, as the SI captions refer to them."""
    for letter, ax in zip("abcdefgh", np.ravel(axs)):
        ax.text(*offset, rf"\textbf{{{letter}}}", transform=ax.transAxes,
                fontsize=17, va="bottom", ha="left")


# ------------------------------------------------------------------ Fig. S6: noise and hit threshold
def si_fig06_noise(df):
    """Raw noise distributions against a matched Gaussian, and normalised noise against the threshold."""
    use_paper_style(14)
    noise = df[df["group"] == "noise"]
    fig, axs = plt.subplots(2, 2, figsize=(11, 7.5), dpi=300)

    for ax, T in zip(axs[0], TEMPS):
        raw = noise[f"activity_{T}"].dropna().to_numpy()
        ax.hist(raw, bins=15, color=TCOLOR[T], edgecolor="black", density=True, alpha=0.8)
        grid = np.linspace(raw.min() - 0.3, raw.max() + 0.6, 400)
        ax.plot(grid, stats.norm.pdf(grid, raw.mean(), raw.std(ddof=1)), "k--", lw=2,
                label="Gaussian fit")
        ax.set_xlabel(f"Raw activity at ${T}^\\circ$C")
        ax.set_ylabel("Density")
        ax.set_title(f"$n = {len(raw)}$ detected, skewness ${stats.skew(raw):.2f}$")
        ax.legend(fontsize=11)

    for ax, T in zip(axs[1], TEMPS):
        norm = noise[f"norm_{T}"].dropna().to_numpy()
        ax.hist(norm, bins=15, color=TCOLOR[T], edgecolor="black", alpha=0.8)
        ax.axvline(HIT_THRESHOLD, color="black", ls="-.", lw=2, label="Hit threshold")
        ax.axvline(np.quantile(norm, 0.999), color="0.35", ls=":", lw=2,
                   label=r"$99.9$th percentile")
        ax.set_xlabel(f"Normalised activity at ${T}^\\circ$C")
        ax.set_ylabel("Frequency")
        ax.set_xlim(right=max(HIT_THRESHOLD + 0.35, norm.max() + 0.35))
        ax.set_ylim(top=1.25 * ax.get_ylim()[1])
        ax.legend(fontsize=11, loc="upper left")

    panel_labels(axs)
    fig.tight_layout()
    save_figure(fig, "si_fig06_noise_threshold")
    plt.close(fig)


# ------------------------------------------------------------- Fig. S7: the 355-sequence reference set
def si_fig07_reference(df, lambert):
    """Coverage of the reference group, its normalised activities, and the comparison with Lambert et al."""
    use_paper_style(14)
    ref = df[df["group"] == "overlap"]
    fig, axs = plt.subplots(2, 3, figsize=(15, 8), dpi=300)

    for ax, col, label in zip(axs[0], ["dist_azo", "dca_energy", "sb_37"],
                              [r"$d_H(\mathrm{Azo})$", r"$\mathcal{H}_\mathrm{gen}$",
                               r"$\mathcal{S}_\mathrm{sb}$"]):
        ax.hist(ref[col].dropna(), bins=20, color="firebrick", edgecolor="black", alpha=0.8)
        ax.set_xlabel(label)
        ax.set_ylabel("Count")

    for ax, T in zip(axs[1][:2], TEMPS):
        act = ref[f"norm_{T}"].dropna().to_numpy()
        n_undet = int(ref[f"norm_{T}"].isna().sum())
        ax.hist(act, bins=25, color=TCOLOR[T], edgecolor="black", alpha=0.8)
        ax.axvline(HIT_THRESHOLD, color="black", ls="-.", lw=2, label="Hit threshold")
        ax.axvline(0.0, color="0.35", ls=":", lw=2, label="Reference mean")
        ax.set_xlabel(f"Normalised activity at ${T}^\\circ$C")
        ax.set_ylabel("Frequency")
        ax.set_title(f"${len(act)}$ detected, ${n_undet}$ undetected")
        ax.legend(fontsize=11, loc="upper left")

    ax = axs[1][2]
    both = ref.merge(lambert, on="sequence", suffixes=("", "_lam"))
    keep = both["activity_37_lam"].notna() & both["norm_37"].notna()
    x = both.loc[keep, "activity_37_lam"].to_numpy()
    y = both.loc[keep, "norm_37"].to_numpy()
    pearson = np.corrcoef(x, y)[0, 1]
    spearman = stats.spearmanr(x, y).statistic
    ax.scatter(x, y, s=18, color="firebrick", alpha=0.7, edgecolor="none")
    ax.set_xlabel(r"Activity at $37^\circ$C, Lambert \textit{et al.}")
    ax.set_ylabel(r"Normalised activity at $37^\circ$C, this work")
    ax.set_title(f"$n = {int(keep.sum())}$, $r = {pearson:.2f}$, $\\rho = {spearman:.2f}$")

    panel_labels(axs)
    fig.tight_layout()
    save_figure(fig, "si_fig07_reference_group")
    plt.close(fig)
    return dict(n=int(keep.sum()), pearson=pearson, spearman=spearman)


# ----------------------------------------------------------------------------- Table S4: read statistics
def si_table04(df):
    rows = []
    for group, sub in df.groupby("group"):
        for T in TEMPS:
            ref, sel, norm = sub[f"ref_count_{T}"], sub[f"sel_count_{T}"], sub[f"norm_{T}"]
            rows.append({
                "group": group, "T": T, "N": len(sub),
                "ref_reads": ref.sum(), "sel_reads": sel.sum(),
                "ref_median": ref.median(), "ref_q1": ref.quantile(0.25), "ref_q3": ref.quantile(0.75),
                "missing": int((ref == 0).sum()),
                "undetected": int(((sel == 0) & (ref > 0)).sum()),
                "detected": int(norm.notna().sum()),
                "hits": int((norm >= HIT_THRESHOLD).sum()),
            })
    table = pd.DataFrame(rows)
    table["undetected_pct"] = 100 * table["undetected"] / table["N"]
    table["hit_pct"] = 100 * table["hits"] / table["N"]
    TABLES.mkdir(parents=True, exist_ok=True)
    table.to_csv(TABLES / "si_table04_read_statistics.csv", index=False)
    return table


def main():
    df = data.load_experiment_scored()
    lambert = data.load_lambert()

    table = si_table04(df)
    print(table.to_string(index=False, float_format=lambda v: f"{v:.1f}"))

    for T in TEMPS:
        raw = df.loc[df.group == "noise", f"activity_{T}"].dropna()
        norm = df.loc[df.group == "noise", f"norm_{T}"].dropna()
        print(f"{T}C noise: raw mean {raw.mean():.2f} sd {raw.std():.2f} skew {stats.skew(raw):.2f} "
              f"Shapiro p {stats.shapiro(raw).pvalue:.1e}; normalised max {norm.max():.2f}, "
              f"99.9th percentile {np.quantile(norm, 0.999):.2f}, hits {(norm >= HIT_THRESHOLD).sum()}")
        ref = df.loc[df.group == "overlap", f"activity_{T}"]
        print(f"{T}C reference mean subtracted: {ref.mean():.4f} (over {ref.notna().sum()} detected)")

    si_fig06_noise(df)
    print("reference-group agreement with Lambert et al.:", si_fig07_reference(df, lambert))


if __name__ == "__main__":
    main()
