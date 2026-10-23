"""Sampling of designs from the effective Hamiltonian (Eq. 6 of the paper).

    H_eff(x | x_start) = beta * H_gen(x) + lambda_S * S_sb(x) + lambda_D * d_H(x_start, x)

The "fast" terms (H_gen and the Hamming-distance penalty) are sampled with Gibbs sweeps, i.e. they are
absorbed in a modified Potts model; the expensive "slow" term S_sb (ViennaRNA) is applied as a Metropolis
accept/reject step on whole blocks of fast sweeps, which preserves detailed balance for the joint Hamiltonian.
"""
import numpy as np
import pandas as pd
import torch
from adabmDCA import decode_sequence, gibbs_sampling

from .config import (AZO_RANGE, H_HIGH_THRESHOLD, MUT_RANGE, SB_HIGH_THRESHOLD, SB_LOW_THRESHOLD, TOKENS)
from . import data, dca, rna

GAP_FIELD = -1e4  # the gap state is never proposed


def biased_model(params, start_onehot, beta, lambda_D):
    """Potts model of the fast part of H_eff: beta * H_gen + lambda_D * d_H(x_start, .)."""
    bias = beta * params["bias"] + lambda_D * start_onehot
    bias[:, TOKENS.index("-")] = GAP_FIELD
    return {"bias": bias, "coupling_matrix": beta * params["coupling_matrix"]}


@torch.no_grad()
def sample_designs(start_seq, params, structure, beta, lambda_S, lambda_D, n_seqs, n_slow=20, n_fast=10,
                   n_workers=None):
    """Run ``n_seqs`` independent chains from ``start_seq`` (197 nt) and return the final sequences.

    One slow step = ``n_fast`` Gibbs sweeps (L attempted mutations each) followed by a Metropolis test on S_sb.
    """
    x0 = dca.to_onehot([start_seq], device=params["bias"].device)[0]
    fast = biased_model(params, x0, beta, lambda_D)
    chains = x0.unsqueeze(0).repeat(n_seqs, 1, 1)

    def slow_score(c):
        return rna.sb_potential(dca.add_flanks(decode_sequence(c.argmax(-1).cpu().numpy(), TOKENS)), structure,
                                37.0, n_workers)

    s0 = slow_score(chains)
    for _ in range(n_slow):
        proposal = gibbs_sampling(chains, fast, nsweeps=n_fast)
        s1 = slow_score(proposal)
        accept = torch.from_numpy(np.random.rand(n_seqs) < np.exp(-lambda_S * (s1 - s0))).to(chains.device)
        chains[accept] = proposal[accept]
        s0 = np.where(accept.cpu().numpy(), s1, s0)
    return dca.add_flanks(decode_sequence(chains.argmax(-1).cpu().numpy(), TOKENS))


def score_candidates(seqs, params, structure, start_seq, azo_seq, n_workers=None):
    df = pd.DataFrame({"sequence": seqs})
    df["dca_energy"] = dca.dca_energy(seqs, params)
    df["sb_37"] = rna.sb_potential(seqs, structure, 37.0, n_workers)
    df["dist_start"] = data.hamming(seqs, start_seq)
    df["dist_azo"] = data.hamming(seqs, azo_seq)
    return df


def select_designs(candidates, group, start_dca_energy, n, seed=42, mut_range=MUT_RANGE, azo_range=AZO_RANGE):
    """Candidates falling in the requested region of the (H_gen, S_sb) plane (Table I) within the budgets."""
    c = candidates.drop_duplicates("sequence")
    c = c[c["dist_start"].between(*mut_range) & c["dist_azo"].between(*azo_range)]
    if group == "lowDCA_lowSB":
        c = c[(c["dca_energy"] < start_dca_energy) & (c["sb_37"] < SB_LOW_THRESHOLD)]
    elif group == "lowDCA_highSB":
        c = c[(c["dca_energy"] < start_dca_energy) & (c["sb_37"] > SB_HIGH_THRESHOLD)]
    elif group == "highDCA_lowSB":
        c = c[(c["dca_energy"] > H_HIGH_THRESHOLD) & (c["sb_37"] < SB_LOW_THRESHOLD)]
    else:
        raise ValueError(group)
    return c.sample(n=min(n, len(c)), random_state=seed)
