"""Loaders for the tables in data/."""
import numpy as np
import pandas as pd

from .config import DATA, HIT_THRESHOLD, TABLES


def load_azoarcus():
    """Azoarcus reference: dict(name, sequence (197 nt), dot_bracket)."""
    row = pd.read_csv(DATA / "azoarcus.csv").iloc[0]
    return row.to_dict()


def load_lambert():
    """All sequences of Lambert et al. with their activities (log10 scale, NaN = no detected activity or not tested).

    ``matched`` marks the 3,429 sequences of the matched dataset used for the analyses of Fig. 2 and SI Figs. 2-3.
    ``distance_azo`` is added here.
    """
    df = pd.read_csv(DATA / "lambert_dataset.csv")
    df["distance_azo"] = hamming(df["sequence"].tolist(), load_azoarcus()["sequence"])
    return df


def load_starting_sequences():
    """The ten starting sequences (rows of the Lambert dataset with a start_id), indexed by start_id."""
    df = load_lambert()
    df = df[df["start_id"].notna()].copy()
    df["start_id"] = df["start_id"].astype(int)
    return df.sort_values("start_id").set_index("start_id")


def load_experiment():
    """Designed sequences, noise/overlap controls and measured activities (NaN = not detected)."""
    return pd.read_csv(DATA / "experiment.csv")


def normalize_activity(df):
    """Shift activities so that the mean over the detected overlap sequences is zero, per temperature.

    Adds columns ``norm_37`` and ``norm_60`` (NaN where no activity was detected).
    """
    df = df.copy()
    overlap = df[df["group"] == "overlap"]
    for T in (37, 60):
        df[f"norm_{T}"] = df[f"activity_{T}"] - overlap[f"activity_{T}"].mean()
    return df


def hamming(seqs, ref):
    """Hamming distance of each string in ``seqs`` from the string ``ref``."""
    a = np.frombuffer("".join(seqs).encode(), dtype="S1").reshape(len(seqs), -1)
    b = np.frombuffer(ref.encode(), dtype="S1")
    return (a != b).sum(axis=1)


def load_experiment_scored():
    """Experiment table merged with the scores of step 1 and with the overlap-normalised activities.

    Requires results/tables/experiment_scores.csv (scripts/01_score_sequences.py).
    Adds ``hit_37`` / ``hit_60``: normalised activity >= HIT_THRESHOLD.
    """
    exp = load_experiment()
    scores = pd.read_csv(TABLES / "experiment_scores.csv").drop(columns=["group", "start_id"])
    df = normalize_activity(exp.merge(scores, on="id", validate="one_to_one"))
    for T in (37, 60):
        df[f"hit_{T}"] = df[f"norm_{T}"] >= HIT_THRESHOLD
    return df
