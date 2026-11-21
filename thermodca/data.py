"""Loaders for the tables in data/."""
import numpy as np
import pandas as pd

from .config import DATA, HIT_THRESHOLD, TABLES


def load_azoarcus():
    """Azoarcus reference: dict(name, sequence (197 nt), dot_bracket)."""
    row = pd.read_csv(DATA / "azoarcus.csv").iloc[0]
    return row.to_dict()


def check_activity_convention(df, temps=(37, 60), key="id", source=""):
    """Check that ``activity_T`` follows the read counts ``ref_count_T`` / ``sel_count_T``.

      ref_count > 0, sel_count > 0   finite activity log10(f_sel / f_ref)
      ref_count > 0, sel_count = 0   undetected: -inf (inactive, kept in every count, quantile and hit rate)
      ref_count = 0                  missing: NaN (no measurement, excluded from every statistic)
      ref_count empty                not tested: NaN (excluded, like missing)
    """
    for T in temps:
        act, ref, sel = df[f"activity_{T}"], df[f"ref_count_{T}"], df[f"sel_count_{T}"]
        expected = np.where(~(ref > 0), "missing", np.where(sel > 0, "finite", "undetected"))
        found = np.where(act.isna(), "missing", np.where(np.isneginf(act), "undetected", "finite"))
        bad = expected != found
        if bad.any():
            raise ValueError(f"{source}: activity_{T} does not follow the read counts for "
                             f"{df.loc[bad, key].tolist()[:5]} ({bad.sum()} rows)")


def load_lambert():
    """All sequences of Lambert et al. with their activities and read counts (log10 scale).

    Activities follow the convention of ``check_activity_convention``: -inf = undetected, NaN = missing
    (no reference read) or not tested (60C only, empty read counts). ``measured_37`` / ``measured_60`` are
    False where the activity is missing or not tested.
    ``matched`` marks the 3,429 sequences used for Fig. 2 and SI Figs. 2-3: tested at both temperatures,
    activity at 37C above the cutoff -2.76 and at least 11 reference reads at 60C.
    ``distance_azo`` is added here.
    """
    df = pd.read_csv(DATA / "lambert_dataset.csv")
    check_activity_convention(df, key="name", source="lambert_dataset.csv")
    for T in (37, 60):
        df[f"measured_{T}"] = df[f"activity_{T}"].notna()
    df["distance_azo"] = hamming(df["sequence"].tolist(), load_azoarcus()["sequence"])
    return df


def load_starting_sequences():
    """The ten starting sequences (rows of the Lambert dataset with a start_id), indexed by start_id."""
    df = load_lambert()
    df = df[df["start_id"].notna()].copy()
    df["start_id"] = df["start_id"].astype(int)
    return df.sort_values("start_id").set_index("start_id")


def load_experiment():
    """Designed sequences, noise/overlap controls and measured activities.

    Activity convention (main text, Methods; SI Sec. S6.3), checked against the read counts by
    ``check_activity_convention``: -inf = undetected, NaN = missing (no reference read).
    """
    df = pd.read_csv(DATA / "experiment.csv")
    check_activity_convention(df, key="id", source="experiment.csv")
    return df


def normalize_activity(df):
    """Shift activities so that the mean over the detected overlap sequences is zero, per temperature.

    Adds columns ``norm_37`` and ``norm_60`` (-inf where undetected, NaN where missing) and
    ``measured_37`` / ``measured_60`` (False where missing).
    """
    df = df.copy()
    overlap = df[df["group"] == "overlap"]
    for T in (37, 60):
        ref = overlap[f"activity_{T}"]
        df[f"norm_{T}"] = df[f"activity_{T}"] - ref[np.isfinite(ref)].mean()
        df[f"measured_{T}"] = df[f"activity_{T}"].notna()
    return df


def hamming(seqs, ref):
    """Hamming distance of each string in ``seqs`` from the string ``ref``."""
    a = np.frombuffer("".join(seqs).encode(), dtype="S1").reshape(len(seqs), -1)
    b = np.frombuffer(ref.encode(), dtype="S1")
    return (a != b).sum(axis=1)


def load_experiment_scored():
    """Experiment table merged with the scores of step 1 and with the overlap-normalised activities.

    Requires results/tables/experiment_scores.csv (scripts/01_score_sequences.py).
    Adds ``hit_37`` / ``hit_60``: normalised activity >= HIT_THRESHOLD (False for undetected and missing
    sequences; compute hit rates over the measured ones only, see ``thermodca.analysis.hit_rate``).
    """
    exp = load_experiment()
    scores = pd.read_csv(TABLES / "experiment_scores.csv").drop(columns=["group", "start_id"])
    df = normalize_activity(exp.merge(scores, on="id", validate="one_to_one"))
    for T in (37, 60):
        df[f"hit_{T}"] = df[f"norm_{T}"] >= HIT_THRESHOLD
    return df
