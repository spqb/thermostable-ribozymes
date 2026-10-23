"""Matplotlib style shared by all figures (LaTeX text, same look as the figures of the paper)."""
import matplotlib.pyplot as plt

from .config import FIGURES


def use_paper_style(font_size=16, **extra):
    """LaTeX rendering with Computer Modern fonts; ``extra`` are additional rcParams."""
    plt.rcdefaults()
    plt.rcParams.update({"text.usetex": True, "font.size": font_size, **extra})


def save_figure(fig, name, formats=("pdf", "png"), **kwargs):
    """Save ``fig`` to results/figures/<name>.<fmt>."""
    FIGURES.mkdir(parents=True, exist_ok=True)
    for fmt in formats:
        fig.savefig(FIGURES / f"{name}.{fmt}", **kwargs)
