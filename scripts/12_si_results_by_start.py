"""Step 12 - SI table and figure of the section "Results by starting sequence".

Reads   data/experiment.csv (through thermodca.data.load_experiment_scored)
        results/tables/table_per_start.csv (step 6: Delta Delta G and Delta H-bonds per start and group)
        results/tables/experiment_scores.csv (S_sb of each design, through load_experiment_scored)
Writes  results/tables/si_table05_results_by_start.csv  activities, hit rates and activity changes for
                                                         each starting sequence and sampling region
        results/figures/si_fig08_results_by_start       hit rates, activity distributions and activity
                                                         changes for each starting sequence

Conventions (as in Sec. S6 of the SI): the activity a_T is the reference-normalised log10 activity; a hit
has a_T >= -1; undetected activities (control reads but no reacted read) are -inf, kept in every count,
hit rate and quantile. The activity change is Delta = a_60 - a_37, computed for the designs that are hits
at 37 C; Delta = -inf when the 60 C activity is undetected. Undetected values are drawn in a separate
band at the bottom of the plots, for display only.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from thermodca import data
from thermodca.config import DESIGN_GROUPS, GROUP_COLORS, GROUP_LABELS, HIT_THRESHOLD, TABLES
from thermodca.plotting import save_figure, use_paper_style

TEMPS = (37, 60)
STARTS = range(1, 11)
UNDET_Y = {"act": -4.3, "delta": -3.4}  # display position of undetected values


def panel_labels(axs, offset=(-0.07, 1.04)):
    for letter, ax in zip("abcdefgh", axs):
        ax.text(*offset, rf"\textbf{{{letter}}}", transform=ax.transAxes, fontsize=17, va="bottom", ha="left")


def with_inf(df):
    """Normalised activities with undetected = -inf (missing stays NaN; there is none among the designs)."""
    df = df.copy()
    for T in TEMPS:
        undetected = (df[f"ref_count_{T}"] > 0) & (df[f"sel_count_{T}"] == 0)
        df[f"a_{T}"] = df[f"norm_{T}"].where(~undetected, -np.inf)
        df[f"undet_{T}"] = undetected
        df[f"missing_{T}"] = ~(df[f"ref_count_{T}"] > 0)
    df["delta"] = df["a_60"] - df["a_37"]
    return df


def q75(x):
    return np.quantile(np.sort(x.to_numpy()), 0.75, method="linear")


def summarize(g):
    h37 = g[g["hit_37"]]
    return pd.Series({
        "N": len(g),
        "undetected_37": int(g["undet_37"].sum()),
        "undetected_60": int(g["undet_60"].sum()),
        "missing": int(g["missing_37"].sum() + g["missing_60"].sum()),
        "median_act_37": g["a_37"].median(),
        "median_act_60": g["a_60"].median(),
        "top25_act_37": q75(g["a_37"]),
        "top25_act_60": q75(g["a_60"]),
        "hits_37": int(g["hit_37"].sum()),
        "hits_60": int(g["hit_60"].sum()),
        "hits_both": int((g["hit_37"] & g["hit_60"]).sum()),
        "hit_rate_37": 100 * g["hit_37"].mean(),
        "hit_rate_60": 100 * g["hit_60"].mean(),
        # retention: fraction of the 37 C hits that are also hits at 60 C
        "retention": 100 * h37["hit_60"].mean() if len(h37) else np.nan,
        "median_delta_hits37": h37["delta"].median() if len(h37) else np.nan,
        "sb_max": g["sb_37"].max() if "sb_37" in g else np.nan,
    })


def results_table(designs):
    t = designs.groupby(["start_id", "group"]).apply(summarize, include_groups=False)
    t = t.reindex(pd.MultiIndex.from_product([list(STARTS), DESIGN_GROUPS], names=["start_id", "group"]))
    struct = pd.read_csv(TABLES / "table_per_start.csv").set_index(["start_id", "group"])
    t = t.join(struct[["ddG_kcal_mol", "dH_bonds"]])
    pooled = designs.groupby("group").apply(summarize, include_groups=False).loc[DESIGN_GROUPS]
    return t, pooled


def strip(ax, designs, col, groups, ylabel, undet_y, width=0.26, show_hit=True):
    rng = np.random.default_rng(0)
    off = np.linspace(-width, width, len(groups))
    for k, g in enumerate(groups):
        for s in STARTS:
            v = designs.loc[(designs.group == g) & (designs.start_id == s), col].to_numpy()
            if col == "delta":
                v = v[designs.loc[(designs.group == g) & (designs.start_id == s), "hit_37"].to_numpy()]
            fin = v[np.isfinite(v)]
            x0 = s + off[k]
            ax.scatter(x0 + rng.uniform(-0.05, 0.05, len(fin)), fin, s=4, color=GROUP_COLORS[g],
                       alpha=0.45, lw=0, rasterized=True)
            nu = int((~np.isfinite(v)).sum())
            if nu:
                ax.scatter(x0 + rng.uniform(-0.05, 0.05, nu), undet_y + rng.uniform(-0.15, 0.15, nu), s=4,
                           color=GROUP_COLORS[g], alpha=0.45, lw=0, rasterized=True)
            if len(v):
                med = np.median(v)
                ax.plot([x0 - 0.07, x0 + 0.07], [max(med, undet_y)] * 2, color="black", lw=1.6, zorder=5)
    ax.axhspan(undet_y - 0.35, undet_y + 0.35, color="0.92", zorder=0)
    lo = ax.get_ylim()[1]
    ax.set_xticks(list(STARTS))
    ax.set_xticklabels([f"st\\_{s}" for s in STARTS])
    ax.set_xlim(0.5, 10.5)
    for s in list(STARTS)[:-1]:
        ax.axvline(s + 0.5, color="0.85", lw=0.8, zorder=0)
    if show_hit:
        ax.axhline(HIT_THRESHOLD, color="black", ls="-.", lw=1.2)
    ax.set_ylabel(ylabel)
    ticks = [t for t in ax.get_yticks() if undet_y + 0.6 < t < lo]
    ax.set_yticks([undet_y] + ticks)
    ax.set_yticklabels(["undet."] + [f"${t:g}$" for t in ticks])
    ax.set_ylim(undet_y - 0.4, lo)


GROUP_TEX = {
    "lowDCA_lowSB": r"$\downarrow\mathcal{H}_\mathrm{gen}\downarrow\mathcal{S}_\mathrm{sb}$",
    "lowDCA_highSB": r"$\downarrow\mathcal{H}_\mathrm{gen}\uparrow\mathcal{S}_\mathrm{sb}$",
    "highDCA_lowSB": r"$\uparrow\mathcal{H}_\mathrm{gen}\downarrow\mathcal{S}_\mathrm{sb}$",
}


def num(x, nd=2):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "---"
    if np.isinf(x):
        return r"$-\infty$"
    return f"{x + 0.0:.{nd}f}".replace("-0.00", "0.00").replace("-", "$-$")


def write_latex(table, pooled, path):
    """Body rows of SI Table S5 (one block per starting sequence, plus the pooled block)."""
    rows = []
    for sid in STARTS:
        for k, g in enumerate(DESIGN_GROUPS):
            r = table.loc[(sid, g)]
            cells = [f"{sid}" if k == 0 else "", GROUP_TEX[g], f"{int(r['N'])}",
                     f"{int(r['undetected_37'])}/{int(r['undetected_60'])}",
                     num(r["median_act_37"]), num(r["median_act_60"]),
                     f"{r['hit_rate_37']:.1f}", f"{r['hit_rate_60']:.1f}",
                     "---" if np.isnan(r["retention"]) else f"{r['retention']:.1f}",
                     num(r["median_delta_hits37"])]
            rows.append(" & ".join(cells) + r" \\")
        rows.append(r"\hline" if sid != list(STARTS)[-1] else r"\hline")
    for g in DESIGN_GROUPS:
        r = pooled.loc[g]
        cells = ["All" if g == DESIGN_GROUPS[0] else "", GROUP_TEX[g], f"{int(r['N'])}",
                 f"{int(r['undetected_37'])}/{int(r['undetected_60'])}",
                 num(r["median_act_37"]), num(r["median_act_60"]),
                 f"{r['hit_rate_37']:.1f}", f"{r['hit_rate_60']:.1f}",
                 f"{r['retention']:.1f}", num(r["median_delta_hits37"])]
        rows.append(" & ".join(cells) + r" \\")
    Path(path).write_text("\n".join(rows) + "\n")


def si_fig08(designs, table):
    use_paper_style(13)
    fig = plt.figure(figsize=(11, 14), dpi=300)
    gs = fig.add_gridspec(4, 2, height_ratios=[1, 1.15, 1.15, 1.0], hspace=0.42, wspace=0.18,
                          top=0.90)
    axa, axb = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    axc, axd, axe = (fig.add_subplot(gs[i, :]) for i in (1, 2, 3))

    w = 0.27
    for ax, T in ((axa, 37), (axb, 60)):
        for k, g in enumerate(DESIGN_GROUPS):
            r = table.xs(g, level="group")[f"hit_rate_{T}"]
            ax.bar(np.array(list(STARTS)) + (k - 1) * w, r.to_numpy(), width=w, color=GROUP_COLORS[g],
                   edgecolor="black", lw=0.4, label=GROUP_LABELS[g])
        ax.set_xticks(list(STARTS))
        ax.set_xticklabels([str(s) for s in STARTS])
        ax.set_xlabel("Starting sequence")
        ax.set_ylim(0, 105)
        ax.set_title(f"${T}\\,^\\circ$C", pad=8)
    axa.set_ylabel(r"Hit rate (\%)")
    handles, labels = axa.get_legend_handles_labels()
    fig.legend(handles, labels, fontsize=12, ncol=3, loc="upper center",
               bbox_to_anchor=(0.5, 0.945), frameon=False)

    for ax, T in ((axc, 37), (axd, 60)):
        strip(ax, designs, f"a_{T}", DESIGN_GROUPS, f"Activity $a_{{{T}}}$", UNDET_Y["act"])
    strip(axe, designs, "delta", DESIGN_GROUPS[:2], r"$\Delta = a_{60} - a_{37}$", UNDET_Y["delta"],
          width=0.16, show_hit=False)
    axe.axhline(0, color="black", lw=0.8, alpha=0.5)
    for k, g in enumerate(DESIGN_GROUPS[:2]):
        axe.scatter([], [], color=GROUP_COLORS[g], s=20, label=GROUP_LABELS[g])
    axe.legend(fontsize=10, loc="upper right", framealpha=0.9)

    panel_labels([axa, axb], offset=(-0.15, 1.04))
    panel_labels([axc, axd, axe], offset=(-0.07, 1.03))
    for ax, l in zip((axc, axd, axe), "cde"):
        ax.texts[-1].set_text(rf"\textbf{{{l}}}")
    save_figure(fig, "si_fig08_results_by_start", bbox_inches="tight")
    plt.close(fig)


def main():
    df = with_inf(data.load_experiment_scored())
    designs = df[df["group"].isin(DESIGN_GROUPS)].copy()
    designs["start_id"] = designs["start_id"].astype(int)
    table, pooled = results_table(designs)
    out = table.round(2)
    out.to_csv(TABLES / "si_table05_results_by_start.csv")
    write_latex(table, pooled, TABLES / "si_table05_results_by_start.tex")
    pd.set_option("display.width", 250)
    print(out.to_string())
    print(pooled.round(2).to_string())

    # consistency of the pattern across starts
    j, d = (table.xs(g, level="group") for g in ("lowDCA_lowSB", "lowDCA_highSB"))
    print("joint > DCA at 60C hit rate:", int((j.hit_rate_60 > d.hit_rate_60).sum()), "/ 10")
    print("joint > DCA retention:", int((j.retention > d.retention).sum()), "/ 10")
    print("joint > DCA median delta:", int((j.median_delta_hits37 > d.median_delta_hits37).sum()), "/ 10")
    print("joint-DCA median delta shift per start:", (j.median_delta_hits37 - d.median_delta_hits37).round(2).tolist())
    print("Wilcoxon signed-rank on per-start hit_rate_60 (joint vs DCA):",
          stats.wilcoxon(j.hit_rate_60, d.hit_rate_60))
    print("Spearman across starts, joint hit60 vs DCA hit60:", stats.spearmanr(j.hit_rate_60, d.hit_rate_60))
    # per-start Mann-Whitney on Delta (hits at 37 C), joint vs DCA-favoured
    for s in STARTS:
        x = designs[(designs.start_id == s) & designs.hit_37 & (designs.group == "lowDCA_lowSB")]["delta"]
        y = designs[(designs.start_id == s) & designs.hit_37 & (designs.group == "lowDCA_highSB")]["delta"]
        print(s, len(x), len(y), "MWU p = %.2g" % stats.mannwhitneyu(x, y, alternative="greater").pvalue)
    si_fig08(designs, table)


if __name__ == "__main__":
    main()
