"""Step 4 - figures based on the data of Lambert et al.

Reads   results/tables/lambert_scores.csv (step 1), data/lambert_dataset.csv
Writes  results/figures/fig2_sb_vs_hgen_37C        main text Fig. 2 (S_sb at 37C)
        results/figures/si_fig2_active_fraction_37C  SI Fig. 2: S_sb and Delta G at 37C
        results/figures/si_fig3_active_fraction_60C  SI Fig. 3: S_sb and Delta G at 60C
        results/tables/lambert_correlations.csv      Pearson correlations of the scores with the distance
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import pearsonr

from gpsopt import data
from gpsopt.config import TABLES
from gpsopt.plotting import save_figure, use_paper_style

# Activity thresholds on the original scale of Lambert et al. (SI: "illustrative cutoff")
ACTIVITY_THRESHOLD_37 = -2.76
ACTIVITY_THRESHOLD_60 = -2.0
AZO_RANGE = (25, 30)
MODELS = ("dca_f", "old_dca")
SB_THRESHOLD = {37: -33.0, 60: -30.0}  # S_sb thresholds shown in Fig. 2 / SI Fig. 3


def load_table():
    lam = data.load_lambert()
    scores = pd.read_csv(TABLES / "lambert_scores.csv")
    return lam.merge(scores.drop(columns=["model", "start_id"]), on="name")


def active_fraction_panel(ax, df, T):
    """Scatter of S_sb vs H_gen coloured by activity at 60C, with the active fraction above/below the cut."""
    thr = SB_THRESHOLD[T]
    sb = df[f"sb_{T}"]
    active = df["active_60"]
    below, above = sb <= thr, sb > thr
    pct_below = 100 * (active & below).sum() / max(below.sum(), 1)
    pct_above = 100 * (active & above).sum() / max(above.sum(), 1)

    kw = dict(alpha=1.0, edgecolors="k", linewidths=0.3)
    ax.scatter(df.loc[active, "dca_energy"], sb[active], c="forestgreen", label="Active at 60°C", **kw)
    ax.scatter(df.loc[~active, "dca_energy"], sb[~active], c="coral", label="Inactive at 60°C", **kw)
    ax.axhline(thr, color="black", linestyle="-.", linewidth=1.0, label=r"$\mathcal{S}_\mathrm{sb}$" + f" = {thr:g}")
    ax.text(0.975, 0.08, f"Active fraction: {pct_below:.1f}" + r"\%", transform=ax.transAxes, va="top", ha="right",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="#d9f2d9", edgecolor="#2f6f2f", alpha=0.9))
    ax.text(0.975, 0.92, f"Active fraction: {pct_above:.1f}" + r"\%", transform=ax.transAxes, va="bottom", ha="right",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="#ffe6cc", edgecolor="#b36200", alpha=0.9))
    ax.set_xlabel(r"$\mathcal{H}_\mathrm{gen}$")
    ax.set_ylabel(r"$\mathcal{S}_\mathrm{sb}$")
    ax.set_ylim(bottom=-45)
    ax.legend()
    return pct_above, pct_below


def delta_g_panel(ax, df, T):
    active = df["active_60"]
    kw = dict(alpha=1.0, edgecolors="k", linewidths=0.3)
    ax.scatter(df.loc[active, "dca_energy"], df.loc[active, f"dg_{T}"], c="forestgreen", label="Active at 60°C", **kw)
    ax.scatter(df.loc[~active, "dca_energy"], df.loc[~active, f"dg_{T}"], c="coral", label="Inactive at 60°C", **kw)
    ax.set_xlabel(r"$\mathcal{H}_\mathrm{gen}$")
    ax.set_ylabel(r"$\Delta G$")
    ax.legend()


def main():
    df = load_table()

    # ---- correlations with the distance from Azoarcus (all sequences and 25-30 subset)
    rows = []
    sel = (
        df["distance_azo"].between(*AZO_RANGE)
        & df["name"].str.contains("|".join(MODELS))
        & (df["activity_37"] > ACTIVITY_THRESHOLD_37)
    )
    sel_all = df["distance_azo"].between(0, 200) & df["name"].str.contains("|".join(MODELS)) & (
        df["activity_37"] > ACTIVITY_THRESHOLD_37
    )
    for label, sub in [("all", df), ("dca_models_active_37", df[sel_all]), ("dca_models_active_37_azo_25-30", df[sel])]:
        for col in ["sb_37", "sb_60", "dg_37", "dg_60", "dca_energy"]:
            r, p = pearsonr(sub[col], sub["distance_azo"])
            rows.append(dict(subset=label, score=col, n=len(sub), pearson_r=r, p_value=p))
    pd.DataFrame(rows).to_csv(TABLES / "lambert_correlations.csv", index=False)

    # ---- sequences active at 37C, 25-30 substitutions from Azoarcus, from the DCA models
    sub = df[sel].copy().reset_index(drop=True)
    sub["active_60"] = sub["activity_60"] > ACTIVITY_THRESHOLD_60
    print(f"{len(sub)} sequences; {sub['active_60'].sum()} active at 60C")

    # Fig. 2 (main text): S_sb at 37C
    use_paper_style(18)
    fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
    active_fraction_panel(ax, sub, 37)
    fig.tight_layout()
    save_figure(fig, "fig2_sb_vs_hgen_37C")
    plt.close(fig)

    # SI Figs. 2 and 3: S_sb and Delta G computed at 37C and 60C
    for T, name in [(37, "si_fig2_active_fraction_37C"), (60, "si_fig3_active_fraction_60C")]:
        fig, axs = plt.subplots(1, 2, figsize=(16, 6), dpi=300, sharex=True)
        pa, pb = active_fraction_panel(axs[0], sub, T)
        delta_g_panel(axs[1], sub, T)
        fig.tight_layout()
        save_figure(fig, name)
        plt.close(fig)
        print(f"T={T}C: active fraction above/below S_sb threshold: {pa:.1f}% / {pb:.1f}%")


if __name__ == "__main__":
    main()
