"""Analysis helpers shared by several scripts."""
import numpy as np

from .config import FLANK_3, FLANK_5  # noqa: F401

INACTIVE_DELTA = -3.5  # delta-activity assigned to hits at 37C with no detected activity at 60C


def delta_activity_table(df):
    """Designs of the two low-H_gen groups with a hit at 37C and their loss of activity upon heating.

    delta = (normalised activity at 60C) - (normalised activity at 37C); designs with no detected activity
    at 60C get delta = INACTIVE_DELTA.
    """
    hits = df[df["group"].isin(["lowDCA_lowSB", "lowDCA_highSB"]) & df["hit_37"]].copy()
    hits["delta_activity"] = hits["norm_60"] - hits["norm_37"]
    hits.loc[hits["norm_60"].isna(), "delta_activity"] = INACTIVE_DELTA
    return hits


def mutation_frequency(seqs, reference):
    """Per-position fraction of ``seqs`` differing from ``reference``."""
    arr = np.frombuffer("".join(seqs).encode(), dtype="S1").reshape(len(seqs), -1)
    return (arr != np.frombuffer(reference.encode(), dtype="S1")).mean(axis=0)
