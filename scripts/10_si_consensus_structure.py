"""Step 10 - SI figure of the Azoarcus consensus secondary structure (SI section "Model and thermodynamic scores").

Reads   data/azoarcus.csv (sequence and consensus dot-bracket, 197 nt)
Writes  results/figures/si_fig_consensus_structure  consensus structure with the contact set and the
                                                    positions covered by the DCA model

Layout: ViennaRNA's naview coordinates of the nested part of the structure (ViennaRNA cannot lay out a
pseudoknot); the six pseudoknot pairs are drawn as dashed chords on top of that layout. Nucleotides are
coloured by whether the DCA model covers them: the constant flanks GUG (construct positions 1-3) and the
terminal G (position 197) are outside the 193 aligned columns.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt
import numpy as np
import RNA
from matplotlib.lines import Line2D

from thermodca import data, rna
from thermodca.config import FLANK_3, FLANK_5
from thermodca.plotting import save_figure, use_paper_style

MODELLED = "#4575b4"   # positions 4-196: the 193 columns of the DCA model
FLANK = "#d73027"      # constant flanks, not modelled
NESTED = "#6b6b6b"     # nested consensus pairs
PSEUDO = "#1a9850"     # pseudoknot pairs


def layout(structure):
    """2-D coordinates (197, 2) of the nested structure, as used by ViennaRNA's RNAplot."""
    xy = RNA.naview_xy_coordinates(rna.remove_pseudoknot(structure))
    return np.array([[p.X, p.Y] for p in xy])[: len(structure)]


def main():
    azo = data.load_azoarcus()
    seq, structure = azo["sequence"], azo["dot_bracket"]
    n = len(structure)
    xy = layout(structure)

    nested = rna.structure_pairs(structure, nested_only=True)
    pseudo = sorted(set(rna.structure_pairs(structure)) - set(nested))
    modelled = np.zeros(n, dtype=bool)
    modelled[len(FLANK_5) : n - len(FLANK_3)] = True

    use_paper_style(16)
    fig, ax = plt.subplots(figsize=(9, 9), dpi=300)

    ax.plot(xy[:, 0], xy[:, 1], color="black", linewidth=0.7, zorder=1)
    for i, j in nested:
        ax.plot(*xy[[i, j]].T, color=NESTED, linewidth=0.9, zorder=2)
    for i, j in pseudo:
        ax.plot(*xy[[i, j]].T, color=PSEUDO, linewidth=1.4, linestyle=(0, (3, 2)), zorder=3)

    for mask, color in [(modelled, MODELLED), (~modelled, FLANK)]:
        ax.scatter(xy[mask, 0], xy[mask, 1], s=26, c=color, edgecolors="white", linewidths=0.4, zorder=4)

    # position numbers every 20 nt, pushed outwards along the local normal of the backbone
    centre = xy.mean(axis=0)
    for pos in range(19, n, 20):
        d = xy[min(pos + 1, n - 1)] - xy[max(pos - 1, 0)]
        normal = np.array([-d[1], d[0]])
        normal = normal / (np.linalg.norm(normal) or 1.0)
        if normal @ (xy[pos] - centre) < 0:
            normal = -normal
        ax.annotate(str(pos + 1), xy[pos] + 20 * normal, va="center", ha="center", fontsize=11, color="#404040")

    for pos, text, dx in [(0, r"5$'$", -18), (n - 1, r"3$'$", 18)]:
        ax.annotate(text, xy[pos], xytext=(xy[pos, 0] + dx, xy[pos, 1]), va="center", ha="center", fontsize=18)

    handles = [
        Line2D([], [], marker="o", color="none", markerfacecolor=MODELLED, markersize=9,
               label=r"modelled positions ($4$--$196$)"),
        Line2D([], [], marker="o", color="none", markerfacecolor=FLANK, markersize=9,
               label=r"constant flanks ($1$--$3$, $197$)"),
        Line2D([], [], color=NESTED, linewidth=1.2, label=rf"nested consensus pairs (${len(nested)}$)"),
        Line2D([], [], color=PSEUDO, linewidth=1.6, linestyle=(0, (3, 2)),
               label=rf"pseudoknot pairs (${len(pseudo)}$)"),
    ]
    ax.legend(handles=handles, loc="upper left", frameon=False, fontsize=14, handlelength=1.8)

    ax.set_aspect("equal")
    ax.axis("off")
    fig.tight_layout()
    save_figure(fig, "si_fig_consensus_structure", bbox_inches="tight")
    plt.close(fig)

    print(f"{n} nt, {len(nested)} nested + {len(pseudo)} pseudoknot consensus pairs")
    print("pseudoknot pairs (1-based construct positions):",
          ", ".join(f"{i + 1}-{j + 1}" for i, j in pseudo))
    print("GC content of the consensus pairs:", rna.base_pair_counts(seq, structure))


if __name__ == "__main__":
    main()
