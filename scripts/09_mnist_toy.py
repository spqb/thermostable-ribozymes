"""Step 9 - toy example of the effective Hamiltonian on MNIST (Fig. 1 of the paper).

A restricted Boltzmann machine trained on the zeros and eights of binarised MNIST defines H_gen; a classifier
gives p(8|x) and the steering potential is S_class = -log p(8|x). Starting from one zero ("Start"), samples at a
a fixed Hamming distance populate the four regions of the (H_gen, S_class) plane.

Reads   data/mnist_toy_samples.npz   precomputed samples (images, H_gen, p(8)) of the four regions and the start image
Writes  results/figures/fig1_mnist_phase_space   (y position -p(8|x), tick labels S_class = -log p(8|x))
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.offsetbox import AnnotationBbox, OffsetImage

from thermodca.config import DATA
from thermodca.plotting import save_figure, use_paper_style

# key in the npz file, plot colour, selection of the samples shown, label, offset (points) of the example image
REGIONS = [
    ("lowH_highS", "C0", lambda p: p < 0.3, r"$\downarrow \mathcal{H}_\mathrm{gen} \, \uparrow \mathcal{S}_\mathrm{class}$", (85, -30)),
    ("lowH_lowS", "C1", lambda p: p > 0.8, r"$\downarrow \mathcal{H}_\mathrm{gen} \, \downarrow \mathcal{S}_\mathrm{class}$", (85, 70)),
    ("highH_highS", "C4", lambda p: p < 0.3, r"$\uparrow \mathcal{H}_\mathrm{gen} \, \uparrow \mathcal{S}_\mathrm{class}$", (-95, 10)),
    ("highH_lowS", "C2", lambda p: p > 0.7, r"$\uparrow \mathcal{H}_\mathrm{gen} \, \downarrow \mathcal{S}_\mathrm{class}$", (-95, 70)),
]
EXAMPLE = 1  # index of the sample whose image is displayed


def digit_image(flat, pad=2):
    """Black digit on white background with a white frame."""
    return np.pad(1.0 - flat.reshape(28, 28).astype(float), pad, mode="constant", constant_values=1.0)


def main():
    use_paper_style(16, **{"xtick.labelsize": 22, "ytick.labelsize": 22, "axes.labelsize": 20})
    d = np.load(DATA / "mnist_toy_samples.npz")
    fig, ax = plt.subplots(figsize=(8, 8), dpi=300)
    xs, ys = [], []
    arrow = dict(arrowstyle="->", mutation_scale=10, linewidth=1.0, color="black")
    box = dict(facecolor="white", edgecolor="black", boxstyle="square,pad=0.1", linewidth=0.8)
    for key, color, keep, label, offset in REGIONS:
        mask = keep(d[f"{key}__p8"])
        x = d[f"{key}__Hgen"][mask]
        y = -d[f"{key}__p8"][mask]  # plotted position: -p(8|x); the ticks are labelled with S_class = -log p(8|x)
        xs.append(x), ys.append(y)
        sns.kdeplot(x=x, y=y, fill=True, levels=5, thresh=0.05, bw_adjust=0.9, cut=2, color=color, alpha=0.30, ax=ax,
                    zorder=1)
        ax.scatter(x, y, s=40, edgecolors="black", linewidths=0.3, label=label, zorder=3, c=color)
        ax.scatter(x[EXAMPLE], y[EXAMPLE], color=color, edgecolors="black", linewidths=0.6, s=100, zorder=4)
        image = OffsetImage(digit_image(d[f"{key}__samples"][mask][EXAMPLE]), zoom=3, cmap="gray",
                            interpolation="nearest")
        ax.add_artist(AnnotationBbox(image, (x[EXAMPLE], y[EXAMPLE]), frameon=True, bboxprops=box, xybox=offset,
                                     xycoords="data", boxcoords="offset points", arrowprops=arrow))
    xs, ys = np.concatenate(xs), np.concatenate(ys)
    pad = max(8.0, 0.05 * (xs.max() - xs.min()))
    ax.set_xlim(xs.min() - pad, xs.max() + pad)
    ax.set_ylim(-1.1, 0.3)
    # linear axis in -p with equally spaced ticks; each tick is labelled with S_class = -log p at that position
    p_ticks = np.linspace(1.0, 0.0, 6)
    ax.set_yticks(-p_ticks)
    ax.set_yticklabels([r"$\infty$" if p == 0 else "$0$" if p == 1 else f"${-np.log(p):.2f}$" for p in p_ticks])
    ax.set_xlabel(r"$\mathcal{H}_\mathrm{gen}$")
    ax.set_ylabel(r"$\mathcal{S}_\mathrm{class} = -\log \mathrm{Prob}(\mathrm{class}=8)$")
    ax.legend(loc="lower center", bbox_to_anchor=(0.27, 1.03), ncol=2, frameon=True, borderaxespad=0.0,
              columnspacing=1.0, handletextpad=0.3, fontsize=14)
    ax_start = ax.inset_axes([0.80, 1.015, 0.2, 0.2], transform=ax.transAxes)
    ax_start.imshow(digit_image(d["start"]), cmap="gray", vmin=0, vmax=1, interpolation="nearest")
    ax_start.set_xticks([]), ax_start.set_yticks([])
    for spine in ax_start.spines.values():
        spine.set_linewidth(0.8)
    ax_start.set_title("Start", fontsize=18, pad=1.5)
    # caption + arrow pointing to the Start image: it is not repeated in the caption
    ax.text(0.705, 1.16, r"$\leq 128$ pixel flips" "\nfrom", transform=ax.transAxes, ha="center", va="center", fontsize=14.5,
            linespacing=1.0)
    ax.annotate("", xy=(0.795, 1.075), xycoords="axes fraction", xytext=(0.655, 1.075), textcoords="axes fraction",
                annotation_clip=False, arrowprops=dict(arrowstyle="-|>", color="0.2", linewidth=1.8, mutation_scale=22,
                                                       shrinkA=0, shrinkB=0))
    fig.tight_layout()
    save_figure(fig, "fig1_mnist_phase_space")
    plt.close(fig)


if __name__ == "__main__":
    main()
