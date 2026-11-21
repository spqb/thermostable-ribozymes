"""Step 13 - SI figure justifying the two activity cutoffs of Lambert et al. (SI section S6).

Reads   data/lambert_dataset.csv (through thermodca.data.load_lambert)
Writes  results/figures/si_fig_lambert_cutoffs  activity against the distance from Azoarcus, at both
                                                temperatures, with the cutoffs -2.76 and -2/ln(10)

The two cutoffs are taken from Ref. Lambert et al. and are used throughout this work to call a
variant of that dataset active. They are illustrative, but they are not arbitrary: sequences far
enough from the Azoarcus wild type cannot splice, so the activity they still show measures the
background of the assay, and both cutoffs sit at the top edge of that background. The figure shows
this twice: as the activity of every variant against its distance from the wild type (top row), and
as the distribution of the activity of all variants against that of the distant ones (bottom row).

Undetected activities (-inf in the data) are censored, not zero, and are left out, as are missing and
untested ones (NaN): the figure is about the distribution of the detected activities.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from thermodca import data
from thermodca.plotting import save_figure, use_paper_style

# Cutoffs of Ref. Lambert et al., in log10 units. At 60 C the published value is -2 in natural-log
# units, i.e. -2/ln(10).
CUTOFF = {37: -2.76, 60: -2 / np.log(10)}
TCOLOR = {37: "cornflowerblue", 60: "firebrick"}
FAR = 60          # substitutions beyond which no variant of the dataset is expected to splice
BIN_WIDTH = 10    # distance bins of the running percentile


def _running_percentile(dist, act, percentile=99, bin_width=BIN_WIDTH):
    """Per-bin percentile of the detected activity, as (bin edges, values) ready for a step plot."""
    edges = np.arange(0, dist.max() + bin_width, bin_width)
    values = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        inside = act[(dist >= lo) & (dist < hi)]
        values.append(np.percentile(inside, percentile) if len(inside) >= 20 else np.nan)
    return edges, np.array(values)


def _density_scatter(ax, dist, act, color, bin_width=5):
    """Scatter with a per-bin transparency, so that the sparse distant variants stay visible.

    The dataset is extremely uneven along the distance axis: tens of thousands of variants within
    the first 40 substitutions and a handful beyond 60. A single alpha therefore either saturates
    the near bins or hides the distant ones, which are exactly the bins the figure is about. Each
    distance bin instead gets an opacity and a marker size set by how many variants it holds.
    """
    rng = np.random.default_rng(0)
    edges = np.arange(0, dist.max() + bin_width, bin_width)
    index = np.clip(np.digitize(dist, edges) - 1, 0, len(edges) - 1)
    for b in np.unique(index):
        inside = index == b
        n = int(inside.sum())
        ax.scatter(dist[inside] + rng.uniform(-0.35, 0.35, n), act[inside],
                   s=float(np.clip(1200 / n, 4, 16)), color=color,
                   alpha=float(np.clip(100 / n, 0.10, 0.85)), edgecolor="none")


def si_fig_lambert_cutoffs(lam):
    use_paper_style(15)
    fig, axs = plt.subplots(2, 2, figsize=(13, 8.6), dpi=300,
                            gridspec_kw={"hspace": 0.33, "wspace": 0.22})
    stats_out = {}

    for col, T in enumerate((37, 60)):
        cut = CUTOFF[T]
        sub = lam[np.isfinite(lam[f"activity_{T}"])]   # detected activities only
        dist = sub["distance_azo"].to_numpy(float)
        act = sub[f"activity_{T}"].to_numpy(float)
        far = act[dist > FAR]

        # --- top row: activity against distance
        ax = axs[0, col]
        ax.axvspan(FAR, dist.max() + 3, color="0.88", zorder=0)
        _density_scatter(ax, dist, act, TCOLOR[T])
        edges, values = _running_percentile(dist, act)
        ax.step(edges, np.append(values, values[-1]), where="post", color="black", lw=2,
                label=r"$99$th percentile per bin")
        ax.axhline(cut, color="black", ls="-.", lw=2, label=f"Cutoff ${cut:.2f}$")
        ax.set_xlabel(r"Hamming distance from \textit{Azoarcus}")
        ax.set_ylabel(f"Activity at ${T}^\\circ$C")
        ax.set_xlim(0, dist.max() + 3)
        ax.set_title(f"${T}^\\circ$C, $n = {len(sub)}$ with detected activity")
        ax.legend(fontsize=12, loc="upper right", framealpha=0.95)

        # --- bottom row: the distribution of the same activities, with the cutoff
        ax = axs[1, col]
        ax.hist(act, bins=np.linspace(act.min(), act.max(), 45), density=True,
                color=TCOLOR[T], alpha=0.65, edgecolor="none")
        ax.axvline(cut, color="black", ls="-.", lw=2, label=f"Cutoff ${cut:.2f}$")
        pct = stats.percentileofscore(far, cut)
        ax.set_xlabel(f"Activity at ${T}^\\circ$C")
        ax.set_ylabel("Density")
        ax.set_title(f"all $n = {len(act)}$ detected activities")
        ax.legend(fontsize=12, loc="upper left", framealpha=0.95)

        stats_out[T] = dict(n=len(sub), n_far=len(far), pct=pct,
                            above=100 * (far > cut).mean(), far_max=far.max())

    for letter, ax in zip("abcd", axs.ravel()):   # row-wise: a,b on top, c,d below
        ax.text(-0.12, 1.06, rf"\textbf{{{letter}}}", transform=ax.transAxes, fontsize=18,
                va="bottom", ha="left")

    fig.tight_layout()
    save_figure(fig, "si_fig_lambert_cutoffs")
    plt.close(fig)
    return stats_out


def main():
    lam = data.load_lambert()
    out = si_fig_lambert_cutoffs(lam)
    for T, s in out.items():
        print(f"{T}C: cutoff {CUTOFF[T]:.2f} is the {s['pct']:.1f}th percentile of the {s['n_far']} "
              f"detected variants beyond {FAR} substitutions ({s['above']:.1f}% above it, "
              f"largest {s['far_max']:.2f}); {s['n']} variants with a detected activity")
    for T in (37, 60):
        sub = lam[np.isfinite(lam[f"activity_{T}"])]
        print(f"  {T}C, fraction above the cutoff by distance bin:", {
            f"{lo}-{lo + 20}": round(100 * (s[f'activity_{T}'] > CUTOFF[T]).mean(), 1)
            for lo in range(0, 160, 20)
            if len(s := sub[(sub.distance_azo >= lo) & (sub.distance_azo < lo + 20)]) > 0})


if __name__ == "__main__":
    main()
