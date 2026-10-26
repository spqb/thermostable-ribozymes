"""Paths, constants and naming conventions shared by all scripts."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RESULTS = ROOT / "results"
TABLES = RESULTS / "tables"
FIGURES = RESULTS / "figures"
DESIGN = RESULTS / "design"
MD = RESULTS / "md_candidates"

# Alphabet of the DCA model: default RNA alphabet of adabmDCA 2.0 (gap first).
TOKENS = "-ACGU"
# Constant flanking nucleotides added to the 193-nt aligned intron: 5'-GUG-[model]-G-3' (197 nt).
FLANK_5, FLANK_3 = "GUG", "G"

N_STARTS = 10

# Experimental groups. Naming follows the paper: "low"/"high" refer to the generative energy H_gen
# and to the structure-based potential S_sb = -sum_{(i,j) in consensus} p_ij (more negative = more
# consensus contacts formed). E.g. lowDCA_lowSB = low H_gen and low S_sb (jointly optimised).
GROUPS = ["lowDCA_lowSB", "lowDCA_highSB", "highDCA_lowSB", "noise", "overlap"]
DESIGN_GROUPS = GROUPS[:3]

GROUP_LABELS = {
    "lowDCA_lowSB": r"$\downarrow \mathcal{H}_\mathrm{gen} \downarrow \mathcal{S}_\mathrm{sb}$",
    "lowDCA_highSB": r"$\downarrow \mathcal{H}_\mathrm{gen} \uparrow \mathcal{S}_\mathrm{sb}$",
    "highDCA_lowSB": r"$\uparrow \mathcal{H}_\mathrm{gen} \downarrow \mathcal{S}_\mathrm{sb}$",
    "noise": "Noise",
    "overlap": "Overlap",
}
# tab10 indices (orange, blue, green, red) used in the paper figures
GROUP_COLORS = {
    "lowDCA_lowSB": "tab:orange",
    "lowDCA_highSB": "tab:blue",
    "highDCA_lowSB": "tab:green",
    "noise": "tab:red",
    "overlap": "firebrick",
}

# Sampling parameters of the three design regions, in the notation of Eq. (6):
#   H_eff = beta * H_gen + lambda_S * S_sb + lambda_D * d_H(x_start, x)
# (beta=0 for the SB-only group, which ignores H_gen). lambda_D = 2.0 is assumed for the three regions;
# see docs/design_parameters.md.
DESIGN_PARAMS = {
    "lowDCA_lowSB": dict(beta=1.5, lambda_S=2.0, lambda_D=2.0),
    "lowDCA_highSB": dict(beta=1 / 0.8, lambda_S=-3.0, lambda_D=2.0),
    "highDCA_lowSB": dict(beta=0.0, lambda_S=3.5, lambda_D=2.0),
}
# Region definitions (Table I of the paper): thresholds on H_gen (start-relative or absolute) and S_sb.
SB_LOW_THRESHOLD = -33.0   # "low S_sb":  S_sb < -33
SB_HIGH_THRESHOLD = -10.0  # "high S_sb": S_sb > -10
H_HIGH_THRESHOLD = 110.0   # "high H_gen": H_gen > 110
N_PER_REGION = {"lowDCA_lowSB": 150, "lowDCA_highSB": 50, "highDCA_lowSB": 50}
MUT_RANGE = (5, 10)        # substitutions w.r.t. the starting sequence
AZO_RANGE = (25, 30)       # substitutions w.r.t. the Azoarcus wild type

# Hit definition: normalised activity >= -1, i.e. at least 10% of the overlap-group activity.
HIT_THRESHOLD = -1.0
