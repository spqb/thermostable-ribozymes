"""RNA secondary-structure utilities based on ViennaRNA (Turner nearest-neighbour model)."""
import os
from concurrent.futures import ProcessPoolExecutor
from functools import partial

import numpy as np
import RNA

_BRACKETS = {"(": ")", "[": "]", "{": "}", "<": ">"}


def structure_pairs(structure, nested_only=False):
    """Base pairs (i, j), i < j, 0-based, of a dot-bracket string (pseudoknot brackets included).

    nested_only: consider only round brackets (the structure without the pseudoknot).
    """
    opening = {"(": ")"} if nested_only else _BRACKETS
    stacks = {close: [] for close in opening.values()}
    pairs = []
    for pos, char in enumerate(structure):
        if char in opening:
            stacks[opening[char]].append(pos)
        elif char in stacks:
            pairs.append((stacks[char].pop(), pos))
    return sorted(pairs)


def contact_matrix(structure, nested_only=False):
    """Symmetric binary contact matrix of a dot-bracket string."""
    cm = np.zeros((len(structure),) * 2, dtype=int)
    for i, j in structure_pairs(structure, nested_only):
        cm[i, j] = cm[j, i] = 1
    return cm


def remove_pseudoknot(structure):
    return structure.replace("[", ".").replace("]", ".")


def _fold_compound(seq, temperature):
    md = RNA.md()
    md.temperature = temperature
    return RNA.fold_compound(seq.upper().replace("T", "U"), md)


def _sb_one(seq, pairs, temperature):
    fc = _fold_compound(seq, temperature)
    fc.pf()
    bpp = fc.bpp()  # 1-based, upper triangle
    return -sum(bpp[i + 1][j + 1] for i, j in pairs)


def _pmap(fn, seqs, n_workers):
    seqs = list(seqs)
    n_workers = n_workers or os.cpu_count()
    if n_workers <= 1 or len(seqs) < 2 * n_workers:
        return [fn(s) for s in seqs]
    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        return list(ex.map(fn, seqs, chunksize=max(1, len(seqs) // (8 * n_workers))))


def sb_potential(seqs, structure, temperature=37.0, n_workers=None):
    """Structure-based potential S_sb = -sum_{(i,j) in consensus} p_ij(x, T).

    p_ij is the equilibrium pairing probability in the full Turner ensemble (no structural constraint).
    More negative values mean that more consensus contacts are formed at equilibrium.
    """
    fn = partial(_sb_one, pairs=structure_pairs(structure), temperature=temperature)
    return np.array(_pmap(fn, seqs, n_workers))


def _fe_one(seq, structure, temperature):
    return _fold_compound(seq, temperature).eval_structure(structure)


def consensus_free_energy(seqs, structure, temperature=37.0, pseudoknot="remove", n_workers=None):
    """Turner free energy (kcal/mol) of each sequence forced in the consensus structure.

    pseudoknot:
        "remove": the pseudoknot pairs are left unpaired (ViennaRNA cannot evaluate pseudoknots).
        "merge":  square brackets are replaced by round ones.
                  This yields a different (nested) pairing pattern; it is kept only to reproduce the
                  Delta G panels of the supplementary figures.
    """
    if pseudoknot == "remove":
        structure = remove_pseudoknot(structure)
    elif pseudoknot == "merge":
        structure = structure.replace("[", "(").replace("]", ")")
    else:
        raise ValueError(pseudoknot)
    fn = partial(_fe_one, structure=structure, temperature=temperature)
    return np.array(_pmap(fn, seqs, n_workers))


_H_BONDS = {"GC": 3, "CG": 3, "AU": 2, "UA": 2, "GU": 2, "UG": 2}


def hydrogen_bonds(seq, structure, nested_only=True):
    """Number of hydrogen bonds of the consensus pairs formed in ``seq`` (GC: 3; AU, GU: 2)."""
    return sum(_H_BONDS.get(seq[i] + seq[j], 0) for i, j in structure_pairs(structure, nested_only))


def base_pair_counts(seq, structure, nested_only=True):
    """Number of G-C and Watson-Crick (G-C + A-U) pairs formed in ``seq`` among the consensus pairs."""
    gc = au = 0
    for i, j in structure_pairs(structure, nested_only):
        pair = {seq[i], seq[j]}
        gc += pair == {"G", "C"}
        au += pair == {"A", "U"}
    return {"gc": gc, "wc": gc + au}
