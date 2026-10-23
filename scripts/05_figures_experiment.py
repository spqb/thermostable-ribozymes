"""Step 5 - figures of the high-throughput experiment.

Reads   data/experiment.csv, results/tables/{experiment,lambert}_scores.csv (step 1)
Writes  results/figures/fig3_phase_space_start5           Fig. 3: designs of start sequence 5 in the (H_gen, S_sb) plane
        results/figures/fig4_activity_violins             Fig. 4: activity of the three design groups and noise
        results/figures/fig5_delta_activity               Fig. 5: loss of activity upon heating
        results/figures/si_fig4_overlap_activity          SI Fig. 4: overlap group, hit threshold
        results/figures/si_fig6_noise_raw_activity        SI Fig. 6: raw activity of the noise group

Activities are normalised by subtracting, for each temperature, the mean activity of the detected overlap-group
sequences. A design is a "hit" if its normalised activity is >= -1 (10% of the overlap activity).
"Not detected" activities are plotted at the sentinel value UNDETECTED.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch

from gpsopt import data
from gpsopt.analysis import INACTIVE_DELTA, delta_activity_table
from gpsopt.config import GROUP_COLORS, GROUP_LABELS, HIT_THRESHOLD, TABLES
from gpsopt.plotting import save_figure, use_paper_style

UNDETECTED = -5.0   # plotting value of sequences without detected activity


# ---------------------------------------------------------------------------------------------- Fig. 3
def fig3_phase_space(df, start_id=5):
    use_paper_style(20)
    lam = pd.read_csv(TABLES / "lambert_scores.csv")
    start = lam[lam["start_id"] == start_id].iloc[0]

    groups = ["lowDCA_highSB", "lowDCA_lowSB", "highDCA_lowSB"]
    sub = df[(df["start_id"] == start_id) & df["group"].isin(groups)].copy()
    hit_rate = {g: [100 * sub.loc[sub.group == g, f"hit_{T}"].mean() for T in (37, 60)] for g in groups}
    sub["Region"] = sub["group"].map(GROUP_LABELS)
    palette = dict(zip([GROUP_LABELS[g] for g in groups], sns.color_palette("deep", n_colors=3)))

    fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
    sns.kdeplot(data=sub, x="dca_energy", y="sb_37", hue="Region", palette=palette, fill=True, alpha=0.35, ax=ax,
                legend=False)
    sns.scatterplot(data=sub, x="dca_energy", y="sb_37", hue="Region", palette=palette, s=30, edgecolor="k",
                    alpha=1.0, ax=ax)
    ax.plot(start["dca_energy"], start["sb_37"], "r*", markersize=15, label="Start")
    ax.set_xlabel(r"$\mathcal{H}_{\mathrm{gen}}$")
    ax.set_ylabel(r"$\mathcal{S}_{\mathrm{sb}}$")

    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    offsets = {"lowDCA_highSB": (0.25, -0.05), "lowDCA_lowSB": (0.23, 0.02), "highDCA_lowSB": (0.13, -0.05)}
    inset_font = 18
    for g in groups:
        color = palette[GROUP_LABELS[g]]
        sg = sub[sub.group == g]
        xfrac = (sg["dca_energy"].median() - xmin) / (xmax - xmin)
        yfrac = (sg["sb_37"].median() - ymin) / (ymax - ymin)
        dx, dy = offsets[g]
        w = h = 0.22
        left = np.clip(xfrac + dx - w / 2, 0.04, 1.5)
        bottom = np.clip(yfrac + dy - h / 2, 0.02, 0.98 - h)
        iax = ax.inset_axes([left, bottom, w, h], transform=ax.transAxes)
        vals = hit_rate[g]
        bars = iax.bar([0, 1], vals, color=[color, color], alpha=0.9)
        iax.set_ylim(0, 100)
        iax.set_xticks([0, 1])
        iax.set_xticklabels([r"$37^\circ\mathrm{C}$", r"$60^\circ\mathrm{C}$"], fontsize=inset_font)
        iax.set_yticks([])
        iax.set_title(r"Hit rate \%", fontsize=inset_font, pad=1)
        iax.tick_params(length=0)
        for bar, v in zip(bars, vals):
            high = v >= 92
            iax.text(bar.get_x() + bar.get_width() / 2, max(v - 4, 2) if high else min(v + 2, 99), f"{v:.1f}",
                     ha="center", va="top" if high else "bottom", fontsize=inset_font, fontweight="bold",
                     color="white" if high else "black")
        for spine in iax.spines.values():
            spine.set_linewidth(0.6)

    handles, labels = ax.get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    ax.legend(unique.values(), unique.keys(), loc="upper right", frameon=True, fontsize=20)
    ax.set_xlim(xmin, xmax + 17)
    fig.tight_layout()
    save_figure(fig, f"fig3_phase_space_start{start_id}")
    plt.close(fig)
    print("Fig. 3 hit rates (37C, 60C):", {g: [round(v, 1) for v in r] for g, r in hit_rate.items()})


# ---------------------------------------------------------------------------------------------- Fig. 4
def fig4_violins(df):
    use_paper_style(24)
    bar_font = 20
    regions = ["lowDCA_lowSB", "lowDCA_highSB", "highDCA_lowSB", "noise"]
    labels = {"37": "37℃", "60": "60℃", **GROUP_LABELS}
    palette = {g: GROUP_COLORS[g] for g in regions}

    fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(16, 7), sharex=True, constrained_layout=True,
                                         gridspec_kw={"height_ratios": [4, 1], "hspace": 0.05}, dpi=300)
    positions, tick_labels, centers, bodies, fracs, texts, colors = [], [], [], [], [], [], []
    pos0 = 0
    for g in regions:
        sub = df[df["group"] == g]
        centers.append(pos0 + 0.5)
        for pos, T in [(pos0, 37), (pos0 + 1, 60)]:
            col = sub[f"norm_{T}"]
            detected = col.dropna().to_numpy()
            vp = ax_top.violinplot([detected], positions=[pos], widths=0.8, showmeans=False, showextrema=False,
                                   showmedians=False)
            for body in vp["bodies"]:
                body.set_facecolor(palette[g])
                body.set_edgecolor("black")
                body.set_alpha(0.85)
                bodies.append((body, pos))
            q = np.percentile(col.fillna(UNDETECTED).to_numpy(), [25, 50, 75])  # undetected count as lowest
            lines = ax_top.hlines(q, pos - 0.28, pos + 0.28, colors="white", linewidth=2.0, zorder=6)
            lines.set_path_effects([pe.Stroke(linewidth=3.2, foreground="black"), pe.Normal()])
            n_missing = int(col.isna().sum())
            fracs.append(n_missing / len(col))
            texts.append(f"{n_missing}/{len(col)}")
            colors.append(palette[g])
            positions.append(pos)
            tick_labels.append(labels[str(T)])
        pos0 += 3

    ymax = ax_top.get_ylim()[1]
    for body, pos in bodies:  # shade and hatch the "hit" part of each violin
        for kw in (dict(facecolor="black", alpha=0.12, edgecolor="none", zorder=3),
                   dict(facecolor="none", edgecolor="black", hatch="////", linewidth=0, zorder=4)):
            patch = plt.Rectangle((pos - 0.4, HIT_THRESHOLD), 0.8, ymax - HIT_THRESHOLD, **kw)
            patch.set_clip_path(body.get_paths()[0], transform=body.get_transform())
            ax_top.add_patch(patch)
    ax_top.set_ylim(bottom=-4.5)

    bars = ax_bot.bar(positions, fracs, width=0.8, color=colors, edgecolor="black", alpha=0.9)
    ax_bot.set_ylabel("Not detected\nfraction")
    ax_bot.set_xticks(positions)
    ax_bot.set_xticklabels(tick_labels)
    top = np.nanmax(fracs)
    ax_bot.set_ylim(0, top + max(0.03, 0.08 * top))
    ax_bot.margins(y=0)
    for bar, text in zip(bars, texts):
        ax_bot.annotate(text, xy=(bar.get_x() + bar.get_width() / 2, bar.get_height()), xytext=(0, 3),
                        textcoords="offset points", ha="center", va="bottom", fontsize=bar_font, clip_on=True)
    lo, hi = ax_bot.get_ylim()
    ax_bot.set_ylim(lo, hi + 0.2)
    for x in (2, 5, 8):
        for ax in (ax_top, ax_bot):
            ax.axvline(x=x, color="gray", linestyle=":", linewidth=1.5, alpha=0.8)
    ax_top.set_ylabel("Activity")
    ax_top.tick_params(axis="x", which="both", labelbottom=False)
    for c, g in zip(centers, regions):
        ax_bot.text(c, -0.55, labels[g], transform=ax_bot.get_xaxis_transform(), ha="center", va="top")
    ax_top.legend(handles=[
        Patch(facecolor="black", edgecolor="black", hatch="////", alpha=0.12, label="Hits"),
        plt.Line2D([0], [0], color="white", linewidth=2.0,
                   path_effects=[pe.Stroke(linewidth=3.2, foreground="black"), pe.Normal()],
                   label="25th, 50th, 75th percentiles")], loc="upper right")
    ax_bot.spines["top"].set_visible(False)
    save_figure(fig, "fig4_activity_violins")
    plt.close(fig)


# ---------------------------------------------------------------------------------------------- Fig. 5
def fig5_delta(df):
    font = 18
    use_paper_style(font)
    sns.set_theme(style="whitegrid", rc={"text.usetex": True, "font.size": font})
    hits = delta_activity_table(df)
    groups = ["lowDCA_lowSB", "lowDCA_highSB"]  # joint (orange), DCA only (blue)
    pal = {"lowDCA_lowSB": sns.color_palette("tab10")[1], "lowDCA_highSB": sns.color_palette("tab10")[0]}

    fig, ax = plt.subplots(figsize=(8, 4), dpi=300)
    vals_all = hits["delta_activity"].to_numpy()
    edges = np.linspace(vals_all.min(), vals_all.max(), 25)
    medians = {}
    for g in groups:
        vals = hits.loc[hits.group == g, "delta_activity"]
        sns.histplot(vals, bins=edges, stat="density", element="step", fill=True, alpha=0.28, linewidth=2.0,
                     color=pal[g], label=GROUP_LABELS[g], ax=ax)
        med = np.median(vals)
        medians[g] = med
        ax.axvline(med, color=pal[g], linestyle="-.", linewidth=2.0)
        joint = g == "lowDCA_lowSB"
        ax.text(med + 0.09 if joint else med - 0.09, ax.get_ylim()[1] * 1.15, f"median={med:.2f}", color=pal[g],
                va="top", ha="left" if joint else "right", fontsize=font,
                bbox=dict(facecolor="white", edgecolor=pal[g], alpha=0.85, boxstyle="round,pad=0.2"))

    # highlight the bin of designs with no detected activity at 60C
    vals_blue = hits.loc[hits.group == "lowDCA_highSB", "delta_activity"].to_numpy()
    hist, _ = np.histogram(vals_blue, bins=edges, density=True)
    k = np.searchsorted(edges, INACTIVE_DELTA, side="right") - 1
    if 0 <= k < len(hist) and hist[k] > 0:
        xc = (edges[k] + edges[k + 1]) / 2
        ax.bar(x=xc, height=hist[k], width=(edges[k + 1] - edges[k]) * 0.95, color=pal["lowDCA_highSB"],
               edgecolor=pal["lowDCA_highSB"], alpha=0.5, hatch="\\\\\\", linewidth=1.0, zorder=6)
        ax.annotate("No activity detected\nat 60°C", xy=(xc, hist[k]), xytext=(xc, hist[k] + 0.2),
                    arrowprops=dict(arrowstyle="->", color=pal["lowDCA_highSB"], linewidth=1.2), ha="center",
                    va="bottom", fontsize=font, color=pal["lowDCA_highSB"],
                    bbox=dict(facecolor="white", edgecolor=pal["lowDCA_highSB"], alpha=0.85, boxstyle="round,pad=0.2"))
    ax.axvline(0.0, color="black", linestyle="-", linewidth=1.0, alpha=0.3, zorder=0)
    ax.set_xlabel(r"$\Delta$ activity $= \mathrm{activity}_{60} - \mathrm{activity}_{37}$", fontsize=font)
    ax.set_ylabel("Density", fontsize=font)
    ax.tick_params(axis="both", which="major", labelsize=font)
    ax.legend(frameon=True, title="Region", fontsize=font, title_fontsize=font)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(False)
    ax.set_xlim(-4.5)
    ax.set_ylim(0, 1.2)
    fig.tight_layout()
    save_figure(fig, "fig5_delta_activity")
    plt.close(fig)
    sns.reset_orig()
    print("Fig. 5 medians:", {g: round(m, 2) for g, m in medians.items()})
    return hits


# ---------------------------------------------------------------------------------------------- SI figures
def si_fig4_overlap(df):
    use_paper_style(16)
    ov = df[df["group"] == "overlap"]
    fig, axs = plt.subplots(1, 2, figsize=(10, 4), dpi=300)
    for ax, T, color in zip(axs, (37, 60), ("cornflowerblue", "firebrick")):
        mean = ov[f"activity_{T}"].mean()
        act = ov[f"activity_{T}"].fillna(UNDETECTED).to_numpy() - mean
        ax.hist(act, bins=30, alpha=0.7, color=color)
        ax.axvline(HIT_THRESHOLD, color="black", linestyle="-.", linewidth=2, label="Hit threshold" if T == 60 else None)
        ax.set_xlabel(f"Activity {T}°C")
        ax.set_ylabel("Frequency")
    axs[1].legend()
    fig.tight_layout()
    for ax, T in zip(axs, (37, 60)):
        mean = ov[f"activity_{T}"].mean()
        counts, edges = np.histogram(ov[f"activity_{T}"].fillna(UNDETECTED).to_numpy() - mean, bins=30)
        x_no = UNDETECTED - mean + 0.15
        y = counts[np.clip(np.digitize([x_no], edges)[0] - 1, 0, len(counts) - 1)]
        ax.annotate("No activity detected", xy=(x_no, y), xytext=(x_no + 0.1, y + 0.3 * counts.max()),
                    arrowprops=dict(arrowstyle="->", lw=1.5, color="black"), fontsize=14, ha="left")
    save_figure(fig, "si_fig4_overlap_activity")
    plt.close(fig)


def si_fig6_noise(df):
    use_paper_style(14)
    noise = df[df["group"] == "noise"]
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.5))
    for ax, T, color in zip(axs, (37, 60), ("skyblue", "lightcoral")):
        ax.hist(noise[f"activity_{T}"].dropna(), bins=15, color=color, edgecolor="black")
        ax.set_title(f"Raw Activity at {T}°C")
        ax.set_xlabel("Activity")
        ax.set_ylabel("Frequency")
    fig.tight_layout()
    save_figure(fig, "si_fig6_noise_raw_activity")
    plt.close(fig)


def main():
    df = data.load_experiment_scored()
    for T in (37, 60):
        thr = df.loc[df.group == "noise", f"norm_{T}"].quantile(0.999)
        print(f"{T}C: 99.9th percentile of noise = {thr:.2f} (hit threshold {HIT_THRESHOLD})")
    fig3_phase_space(df)
    fig4_violins(df)
    fig5_delta(df)
    si_fig4_overlap(df)
    si_fig6_noise(df)


if __name__ == "__main__":
    main()
