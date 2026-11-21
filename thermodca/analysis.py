"""Analysis helpers shared by several scripts.

Activity convention of the experiment (see ``thermodca.data.load_experiment``): undetected activities are
-inf and enter every count, quantile and hit rate; missing activities (no reference read) are NaN and are
excluded from every statistic.
"""
import numpy as np

from .config import FLANK_3, FLANK_5  # noqa: F401

# Display position of the designs with a hit at 37C and no detected activity at 60C (Delta = -inf) in Fig. 5.
# It is a plotting convention only: the Delta of these designs is -inf in every statistic.
INACTIVE_DELTA = -3.5


def quantile(x, q):
    """Linear-interpolation quantile of ``x`` with undetected values (-inf) included and NaN (missing) dropped.

    Returns -inf when the quantile falls among the undetected values (numpy's interpolation would give NaN
    when both neighbours are -inf). Equals ``np.quantile`` on any input without -inf.
    """
    v = np.sort(np.asarray(x, dtype=float))
    v = v[~np.isnan(v)]
    if not len(v):
        return np.nan
    h = (len(v) - 1) * q
    lo = int(np.floor(h))
    if np.isneginf(v[lo]):
        return -np.inf
    return float(np.quantile(v, q))


def hit_rate(g, T):
    """Percentage of hits at temperature T among the measured sequences of ``g`` (missing ones excluded)."""
    measured = g[f"measured_{T}"].to_numpy(bool)
    return 100 * g[f"hit_{T}"].to_numpy(bool)[measured].mean() if measured.any() else np.nan


def delta_activity_table(df):
    """Designs of the two low-H_gen groups with a hit at 37C and their loss of activity upon heating.

    delta_activity = (normalised activity at 60C) - (normalised activity at 37C); it is -inf for designs with
    no detected activity at 60C. Designs with a missing activity at 60C are dropped. ``delta_display`` is
    the same with -inf replaced by INACTIVE_DELTA, for histograms only.
    """
    hits = df[df["group"].isin(["lowDCA_lowSB", "lowDCA_highSB"]) & df["hit_37"] & df["measured_60"]].copy()
    hits["delta_activity"] = hits["norm_60"] - hits["norm_37"]
    hits["delta_display"] = hits["delta_activity"].where(np.isfinite(hits["delta_activity"]), INACTIVE_DELTA)
    return hits


def mutation_frequency(seqs, reference):
    """Per-position fraction of ``seqs`` differing from ``reference``."""
    arr = np.frombuffer("".join(seqs).encode(), dtype="S1").reshape(len(seqs), -1)
    return (arr != np.frombuffer(reference.encode(), dtype="S1")).mean(axis=0)
