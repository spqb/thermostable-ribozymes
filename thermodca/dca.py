"""DCA model of the group I intron family (adabmDCA >= 2.0 conventions)."""
import numpy as np
import torch
from adabmDCA import compute_energy, encode_sequence, load_params, one_hot

from .config import DATA, FLANK_3, FLANK_5, TOKENS


def default_device():
    return torch.device("cuda" if torch.cuda.is_available() else "cpu")


def load_model(path=DATA / "dca_model_lambert.dat", device=None):
    """Load the Potts model (L = 193, alphabet '-ACGU')."""
    return load_params(str(path), tokens=TOKENS, device=device or default_device())


def strip_flanks(seqs):
    """Remove the constant 5'-GUG and 3'-G flanks (197 nt -> 193 nt)."""
    return [s[len(FLANK_5) : len(s) - len(FLANK_3)] for s in seqs]


def add_flanks(seqs):
    return [FLANK_5 + s + FLANK_3 for s in seqs]


def to_onehot(seqs, device=None, flanked=True):
    """One-hot encode sequences (N, 193, 5). ``flanked``: input has the 197-nt form."""
    seqs = list(seqs)
    if flanked:
        seqs = strip_flanks(seqs)
    enc = torch.tensor(encode_sequence(seqs, TOKENS), dtype=torch.long)
    return one_hot(enc, num_classes=len(TOKENS)).to(device or default_device())


@torch.no_grad()
def dca_energy(seqs, params, flanked=True, batch_size=4096):
    """H_gen of each sequence (lower = more compatible with the family)."""
    seqs = list(seqs)
    out = []
    for i in range(0, len(seqs), batch_size):
        x = to_onehot(seqs[i : i + batch_size], device=params["bias"].device, flanked=flanked)
        out.append(compute_energy(x, params).cpu().numpy())
    return np.concatenate(out)
