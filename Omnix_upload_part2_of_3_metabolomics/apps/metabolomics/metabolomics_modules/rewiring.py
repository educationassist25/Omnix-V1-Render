"""
rewiring.py - Differential correlation ("Rewiring Map") engine for MetaboAI Pro.

Fold-change statistics ask "which metabolites go up or down?". This module asks a
different question: "which metabolite *relationships* change between two groups?"
Two metabolites can keep the same average level while their coordination breaks,
appears, or reverses -- a pattern that volcano plots, heatmaps and PCA never show.

Pipeline (all functions are pure: DataFrames/arrays in, DataFrames/dicts out)
---------------------------------------------------------------------------
1. prepare_data()             align samples to two groups, log2 if needed, filter
                              missingness, impute (log2 half-minimum), drop
                              zero-variance features, cap to the most variable
                              features, attach pathway labels (or cluster).
2. differential_correlation() per-group Pearson or Spearman correlation for every
                              feature pair, Fisher z-test on r_B - r_A,
                              Benjamini-Hochberg FDR (statsmodels "fdr_bh", same
                              as the rest of MetaboAI Pro).
3. classify_pairs()           lost / gained / flipped / intact.
4. feature_statistics()       per-feature Welch t-test log2 fold change + BH FDR,
                              used to show which rewired metabolites a fold-change
                              analysis would miss.
5. permutation_test()         global "is there more rewiring than chance?" test:
                              mean |z| over all pairs vs. label permutations
                              (restricted within batch when a batch column is set).
6. build_hypotheses()         groups rewired pairs into pathway modules and scores
                              each one: bootstrap stability, agreement with the
                              other correlation method, within-batch replication.
                              Produces ranked "hypothesis cards" with a claim, an
                              interpretation and a suggested follow-up experiment.

run_rewiring() chains 2-6 and returns a RewiringResult that also knows how to
serialize itself for the interactive view (to_payload) and for downloads.

Statistical notes
-----------------
* Fisher z standard error: sqrt(1/(n_A-3) + 1/(n_B-3)) for Pearson. For Spearman
  the variance is inflated by 1.06 (Fieller, Hartley & Pearson 1957), the usual
  large-sample approximation.
* Differential correlation is data-hungry: with n per group the SE of z_B - z_A is
  sqrt(2/(n-3)), so n=10 per group can only detect |delta z| > ~1.5. The engine
  refuses n < min_per_group (default 5) and the UI warns below 15.
* Interpretations are rule-based text keyed on pathway names. They are hypotheses
  to test, never conclusions; every card links back to the pairs that support it.
"""

from __future__ import annotations

import json
import math
import re
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

CLASSES = ("lost", "gained", "flipped", "intact")
REWIRED = ("lost", "gained", "flipped")
SPEARMAN_VAR_INFLATION = 1.06
UNANNOTATED = "Unannotated"
NOT_IN_LIBRARY = "Not in library"


class RewiringInputError(ValueError):
    """Raised with a user-facing message when the input can't be analysed."""


@dataclass
class RewiringConfig:
    method: str = "pearson"            # "pearson" | "spearman"
    fdr_method: str = "permutation"    # "permutation" (needs n_permutations > 0) | "bh"
    fdr: float = 0.05                  # cutoff for a rewired pair (on q or p, see cutoff_on); up to 1
    cutoff_on: str = "fdr"             # "fdr" (q-value) | "p" (p-value: permutation p when permutations run)
    min_delta: float = 0.45            # minimum |r_B - r_A| for a rewired pair
    flip_min_abs_r: float = 0.4        # both |r| >= this and opposite sign -> flipped
    intact_min_abs_r: float = 0.6      # both |r| >= this, same sign, not rewired -> intact
    n_permutations: int = 200
    n_bootstrap: int = 200
    robust_delta: float = 0.25         # bootstrap keeps direction with |delta r| >= this
    stability_cutoff: float = 0.6
    influence_frac: float = 0.2        # outlier check removes at most this share of a group
    batch_min_per_group: int = 5
    max_modules_scored: int = 30
    max_pairs_per_module: int = 25
    pair_scope: str = "all"            # "all" | "reaction" (<= max_steps) | "direct" (1 step) -- pairs that are tested
    max_steps: int = 2                 # reaction distance for "reaction" scope and for annotation
    pathway_min_size: int = 3          # measured members for a pathway-level (set) test
    seed: int = 0

    def __post_init__(self):
        if self.method not in ("pearson", "spearman"):
            raise RewiringInputError("method must be 'pearson' or 'spearman'.")
        if not (0 < self.fdr <= 1):
            raise RewiringInputError("The rewired-pair cutoff must be above 0 and at most 1.")
        if self.cutoff_on not in ("fdr", "p"):
            raise RewiringInputError("cutoff_on must be 'fdr' or 'p'.")
        if self.fdr_method not in ("permutation", "bh"):
            raise RewiringInputError("fdr_method must be 'permutation' or 'bh'.")
        if not (0 <= self.min_delta <= 2):
            raise RewiringInputError("min_delta must be between 0 and 2.")
        if self.pair_scope not in ("all", "reaction", "direct"):
            raise RewiringInputError("pair_scope must be 'all', 'reaction' or 'direct'.")


# ---------------------------------------------------------------------------
# 1. Data preparation
# ---------------------------------------------------------------------------
@dataclass
class PreparedData:
    X: np.ndarray                  # features x samples, log2 scale, no NaN
    features: list
    pathways: list
    samples: list
    group: np.ndarray              # 0 = group A (reference), 1 = group B
    group_labels: tuple
    batch: Optional[np.ndarray]
    batch_col: Optional[str]
    logged: bool
    pathway_source: str            # "annotation" | "cluster" | "mixed"
    scale_label: str = "log2"      # axis label for values ("log2" unless the user forced no transform)
    imputed: Optional[np.ndarray] = None   # features x samples, True where a value was imputed
    notes: list = field(default_factory=list)
    pathway_library: str = ""              # library/column the ring grouping came from
    sets: dict = field(default_factory=dict)      # set name -> feature indices (pathway-level test)
    gem_ids: list = field(default_factory=list)   # Human-GEM metabolite per feature ('' if unmatched)
    id_map: Optional[pd.DataFrame] = None         # ID standardization + Human-GEM + pathway per feature

    @property
    def idx_a(self) -> np.ndarray:
        return np.flatnonzero(self.group == 0)

    @property
    def idx_b(self) -> np.ndarray:
        return np.flatnonzero(self.group == 1)

    @property
    def n_a(self) -> int:
        return int((self.group == 0).sum())

    @property
    def n_b(self) -> int:
        return int((self.group == 1).sum())


def _display_name(name) -> str:
    text = str(name)
    return text.split("::", 1)[1] if "::" in text else text


def _match_pathways(features, pathways: Optional[pd.Series]) -> list:
    """Look each feature up in a Metabolite -> pathway Series: exact name first, then the
    part after a multi-dataset 'Dataset::' prefix (mirrors the Heatmap's row-annotation
    matching), then a case-insensitive match."""
    if pathways is None:
        return [None] * len(features)
    s = pathways.dropna()
    s = s[s.astype(str).str.strip() != ""]
    s.index = s.index.map(str)
    exact = s.to_dict()
    lower = {k.lower(): v for k, v in exact.items()}
    out = []
    for f in features:
        f = str(f)
        v = exact.get(f)
        if v is None:
            v = exact.get(_display_name(f))
        if v is None:
            v = lower.get(_display_name(f).lower())
        out.append(str(v) if v is not None else None)
    return out


def _cluster_features(X: np.ndarray, group: np.ndarray, seed: int) -> list:
    """Group features by correlation profile when no pathway annotation exists.
    Each group is centred separately first so a mean shift between groups doesn't
    masquerade as correlation."""
    from sklearn.cluster import KMeans

    m = X.shape[0]
    Xc = X.copy()
    for g in (0, 1):
        idx = group == g
        Xc[:, idx] -= Xc[:, idx].mean(axis=1, keepdims=True)
    sd = Xc.std(axis=1, ddof=1, keepdims=True)
    sd[sd == 0] = 1
    Z = Xc / sd
    C = np.clip(Z @ Z.T / (Z.shape[1] - 1), -1, 1)
    k = int(np.clip(round(math.sqrt(m / 2)), 2, 8))
    k = min(k, m)
    labels = KMeans(n_clusters=k, n_init=10, random_state=seed).fit_predict(C)
    sizes = pd.Series(labels).value_counts()
    rank = {lab: i + 1 for i, lab in enumerate(sizes.index)}   # Cluster 1 = largest
    return [f"Cluster {rank[l]}" for l in labels]


def prepare_data(
    data: pd.DataFrame,
    meta: pd.DataFrame,
    group_col: str,
    group_a,
    group_b,
    *,
    pathways: Optional[pd.Series] = None,
    batch_col: Optional[str] = None,
    log_transform: str = "auto",
    max_features: int = 300,
    max_missing_frac: float = 0.3,
    min_per_group: int = 3,
    sample_median_center: bool = False,
    seed: int = 0,
    pathway_library: Optional[str] = None,
    row_annotations: Optional[pd.DataFrame] = None,
    annotation_column: Optional[str] = None,
) -> PreparedData:
    """
    data     : features x samples (rows = metabolites). Raw peak areas or log2 values.
    meta     : one row per sample, indexed by sample name (as in MetaboAI Pro after
               validate_metadata), or with a 'Sample' column.
    group_a  : reference group (e.g. Control); group_b : comparison group.
    pathways : optional Series Metabolite -> pathway/class label (legacy; prefer pathway_library).
    pathway_library : how metabolites are grouped and which sets the pathway-level test uses:
               "KEGG" / "SMPDB" / "LIPID MAPS" (the MSEA libraries, same ID standardization as the
               Pathway Analysis page), "Human-GEM subsystems", "annotation" (annotation_column of
               row_annotations) or "cluster". row_annotations may carry HMDB/KEGG/ChEBI/PubChem
               columns; they are used for ID standardization exactly as on the MSEA page.
    min_per_group : 3 is the floor (a correlation needs 3 points); with 3 per group only
               permutation-based tests are available and their resolution is 1/20.
    log_transform : "auto" (log2 when values look like raw intensities: all >= 0
               and max > 100), "yes", or "no".
    sample_median_center : subtract each sample's median (log scale). Sample-wide
               dilution/loading differences make every pair look positively
               correlated; centring removes that before correlations are computed.
    """
    notes = []
    if data is None or data.empty:
        raise RewiringInputError("The data matrix is empty.")
    if meta is None or meta.empty:
        raise RewiringInputError("Sample metadata is required to assign groups.")
    meta = meta.copy()
    if "Sample" in meta.columns and meta.index.name != "Sample":
        meta = meta.set_index("Sample")
    meta.index = meta.index.map(str)
    if group_col not in meta.columns:
        raise RewiringInputError(f"Metadata has no column '{group_col}'.")
    if str(group_a) == str(group_b):
        raise RewiringInputError("Choose two different groups to compare.")

    df = data.copy()
    df.columns = df.columns.map(str)
    df = df.apply(pd.to_numeric, errors="coerce")

    is_qc = pd.Series(False, index=meta.index)
    if "IsQC" in meta.columns:
        is_qc = meta["IsQC"].astype(str).str.lower().isin(["true", "1", "yes", "qc"])
    labels = meta[group_col].astype(str)
    in_meta = [c for c in df.columns if c in meta.index]
    if not in_meta:
        raise RewiringInputError("None of the data columns match sample names in the metadata.")
    samples = [c for c in in_meta if not is_qc.get(c, False) and labels[c] in (str(group_a), str(group_b))]
    group = np.array([0 if labels[s] == str(group_a) else 1 for s in samples], dtype=int)
    n_a, n_b = int((group == 0).sum()), int((group == 1).sum())
    if n_a == 0 or n_b == 0:
        missing = str(group_a) if n_a == 0 else str(group_b)
        raise RewiringInputError(f"No samples found for group '{missing}' in column '{group_col}'.")
    if min(n_a, n_b) < max(3, min_per_group):
        raise RewiringInputError(
            f"Differential correlation needs at least {max(3, min_per_group)} samples per group "
            f"(found {n_a} {group_a} and {n_b} {group_b})."
        )
    df = df[samples]

    # log2 decision
    vals = df.values[np.isfinite(df.values)]
    if vals.size == 0:
        raise RewiringInputError("The data matrix has no numeric values.")
    if log_transform == "auto":
        do_log = bool(vals.min() >= 0 and vals.max() > 100)
    else:
        do_log = log_transform == "yes"
    if do_log:
        df = df.where(df > 0)
        df = np.log2(df)
        notes.append("Values looked like raw peak areas, so they were log2-transformed (zeros treated as missing).")

    # missingness filter (per group) + log2 half-minimum imputation
    miss_a = df.iloc[:, group == 0].isna().mean(axis=1)
    miss_b = df.iloc[:, group == 1].isna().mean(axis=1)
    keep = (miss_a <= max_missing_frac) & (miss_b <= max_missing_frac)
    n_drop = int((~keep).sum())
    if n_drop:
        notes.append(f"{n_drop} feature(s) dropped for >{max_missing_frac:.0%} missing values in a group.")
    df = df[keep]
    imputed_mask = df.isna()
    n_missing = int(imputed_mask.values.sum())
    if n_missing:
        fill = df.min(axis=1) - 1.0          # log2(min/2)
        df = df.apply(lambda row: row.fillna(fill[row.name]), axis=1)
        notes.append(f"{n_missing} remaining missing value(s) imputed at half the feature minimum.")

    if sample_median_center:
        df = df - df.median(axis=0)
        notes.append("Each sample was median-centred (log scale) to remove dilution/loading differences.")

    # zero variance inside either group -> correlation undefined
    X = df.values.astype(float)
    sd_a = X[:, group == 0].std(axis=1, ddof=1)
    sd_b = X[:, group == 1].std(axis=1, ddof=1)
    ok = (sd_a > 1e-10) & (sd_b > 1e-10)
    if (~ok).sum():
        notes.append(f"{int((~ok).sum())} feature(s) dropped for having no variation within a group.")
    df = df[ok]
    if len(df) > max_features:
        var = df.var(axis=1, ddof=1)
        top = var.sort_values(ascending=False).index[:max_features]
        df = df.loc[[i for i in df.index if i in set(top)]]
        notes.append(f"Analysing the {max_features} most variable of {int(ok.sum())} features.")
    if len(df) < 3:
        raise RewiringInputError("Fewer than 3 usable features remain after filtering.")

    imputed = imputed_mask.loc[df.index].values.astype(bool)
    features = [str(f) for f in df.index]
    X = df.values.astype(float)

    from metabolomics_modules import rewiring_reactions as rr
    sets, id_map, library_label = {}, None, ""
    if pathway_library is not None:
        lib = pathway_library
        fm = rr.map_features(features, row_annotations,
                             library=rr.LIB_ANNOTATION if lib == "annotation" else lib,
                             annotation_column=annotation_column)
        id_map, gem_ids = fm.table, list(fm.table["GEM_ID"])
        notes.extend(fm.notes)
        pw = list(fm.pathway)
        sets = dict(fm.sets)
        library_label = annotation_column if lib == "annotation" else lib
        n_annot = sum(p is not None for p in pw)
        if lib == rr.LIB_CLUSTER or n_annot == 0:
            pw = _cluster_features(X, group, seed)
            sets = {}
            for k, p_ in enumerate(pw):
                sets.setdefault(p_, []).append(k)
            source = "cluster"
            library_label = "correlation clusters"
            if lib != rr.LIB_CLUSTER:
                notes.append(f"No feature matched a '{lib}' set, so features were grouped into correlation clusters.")
        else:
            pw = [p_ if p_ is not None else NOT_IN_LIBRARY for p_ in pw]
            source = "library" if lib not in ("annotation",) else "annotation"
            if n_annot < len(features):
                notes.append(f"{len(features) - n_annot} feature(s) are in no '{library_label}' set and are grouped as "
                             f"'{NOT_IN_LIBRARY}'.")
    else:
        try:
            fm = rr.map_features(features, None, library=rr.LIB_CLUSTER)
            id_map, gem_ids = fm.table, list(fm.table["GEM_ID"])
        except Exception:   # reaction library unavailable: analysis still runs without reaction links
            id_map, gem_ids = None, [""] * len(features)
        pw = _match_pathways(features, pathways)
        n_annot = sum(p is not None for p in pw)
        if n_annot == 0:
            pw = _cluster_features(X, group, seed)
            source = "cluster"
            library_label = "correlation clusters"
            notes.append("No pathway annotation matched, so features were grouped into correlation clusters.")
        else:
            pw = [p if p is not None else UNANNOTATED for p in pw]
            source = "annotation" if n_annot == len(features) else "mixed"
            library_label = "annotation"
            if source == "mixed":
                notes.append(f"{len(features) - n_annot} feature(s) had no pathway annotation and are grouped as '{UNANNOTATED}'.")
        for k, p_ in enumerate(pw):
            if p_ != UNANNOTATED:
                sets.setdefault(p_, []).append(k)

    batch = None
    if batch_col:
        if batch_col not in meta.columns:
            raise RewiringInputError(f"Metadata has no batch column '{batch_col}'.")
        batch = meta.loc[samples, batch_col].astype(str).values

    if min(n_a, n_b) <= 3:
        notes.append(
            f"{min(n_a, n_b)} samples in the smaller group: each correlation rests on 3 points, the analytic "
            "Fisher-z test is undefined, and every p-value comes from label permutations "
            f"({math.comb(n_a + n_b, n_a)} possible label splits, so no p-value can be below "
            f"{(2 if n_a == n_b else 1) / math.comb(n_a + n_b, n_a):.3g}). Nothing can reach FDR q-value < 0.05 at this size; "
            "the ranked correlation changes and pathway statistics are for exploration and for planning a larger study."
        )
    elif min(n_a, n_b) < 15:
        notes.append(
            f"{min(n_a, n_b)} samples in the smaller group: single pairs need large changes (roughly |Δr| > 0.8) "
            "to reach significance; the pathway-level test pools pairs and has more power."
        )

    if id_map is not None:
        id_map = id_map.copy()
        id_map["Ring_Group"] = pw
    return PreparedData(
        X=X, features=features, pathways=pw, samples=samples, group=group, imputed=imputed,
        group_labels=(str(group_a), str(group_b)), batch=batch, batch_col=batch_col,
        logged=do_log, scale_label="log2",   # "no" means the input is already on a log2 scale
        pathway_source=source, notes=notes, pathway_library=library_label,
        sets={k: sorted(v) for k, v in sets.items()}, gem_ids=[g or "" for g in gem_ids], id_map=id_map,
    )


# ---------------------------------------------------------------------------
# 2-4. Correlation, classification, per-feature statistics
# ---------------------------------------------------------------------------
def _standardize(Xs: np.ndarray, method: str) -> np.ndarray:
    if method == "spearman":
        Xs = stats.rankdata(Xs, axis=1)
    Xs = Xs - Xs.mean(axis=1, keepdims=True)
    sd = Xs.std(axis=1, ddof=1, keepdims=True)
    sd[sd == 0] = 1
    return Xs / sd


def correlation_matrix(X: np.ndarray, idx: np.ndarray, method: str = "pearson") -> np.ndarray:
    Z = _standardize(X[:, idx], method)
    C = Z @ Z.T / (Z.shape[1] - 1)
    np.fill_diagonal(C, 1.0)
    return np.clip(C, -1.0, 1.0)


def fisher_se(n_a: int, n_b: int, method: str) -> float:
    v = SPEARMAN_VAR_INFLATION if method == "spearman" else 1.0
    return math.sqrt(v / (n_a - 3) + v / (n_b - 3))


def pair_scale(n_a: int, n_b: int, method: str) -> Optional[float]:
    """Fisher-z standard error of z_B - z_A, or None when a group has <= 3 samples (the
    analytic SE is undefined; the statistic is then the raw difference of Fisher z values and
    all inference is by permutation, where a constant scale changes nothing)."""
    return fisher_se(n_a, n_b, method) if min(n_a, n_b) > 3 else None


def _fisher_z(r):
    return np.arctanh(np.clip(r, -0.999999, 0.999999))


def differential_correlation(prep: PreparedData, method: str = "pearson") -> pd.DataFrame:
    CA = correlation_matrix(prep.X, prep.idx_a, method)
    CB = correlation_matrix(prep.X, prep.idx_b, method)
    iu, ju = np.triu_indices(len(prep.features), k=1)
    ra, rb = CA[iu, ju], CB[iu, ju]
    se = pair_scale(prep.n_a, prep.n_b, method)
    dz = _fisher_z(rb) - _fisher_z(ra)
    if se is not None:
        z = dz / se
        p = 2 * stats.norm.sf(np.abs(z))
    else:
        z = dz
        p = np.full(len(dz), np.nan)
    f = np.array(prep.features, dtype=object)
    pw = np.array(prep.pathways, dtype=object)
    return pd.DataFrame({
        "i": iu, "j": ju,
        "feature_a": f[iu], "feature_b": f[ju],
        "pathway_a": pw[iu], "pathway_b": pw[ju],
        "r_a": ra, "r_b": rb, "delta_r": rb - ra,
        "z": z, "p_value": p, "q_value": _bh(p),
    })


def _bh(p, mask=None) -> np.ndarray:
    """Benjamini-Hochberg over the finite p-values inside `mask`; NaN elsewhere."""
    p = np.asarray(p, dtype=float)
    ok = np.isfinite(p) & (np.ones(len(p), bool) if mask is None else np.asarray(mask, bool))
    q = np.full(len(p), np.nan)
    if ok.any():
        q[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return q


def annotate_reactions(prep: PreparedData, pairs: pd.DataFrame, max_steps: int = 2) -> pd.DataFrame:
    """Reaction link for every pair from the Human-GEM library: steps (0 = none, 1 = direct,
    2 = through one unmeasured intermediate), direction relative to feature_a -> feature_b
    ('a>b', 'b>a', 'both'), intermediate, reactions, subsystems, EC numbers."""
    out = pairs.copy()
    n = len(out)
    cols = {"rxn_steps": np.zeros(n, int), "rxn_dir": [""] * n, "rxn_via": [""] * n, "rxn_names": [""] * n,
            "rxn_ids": [""] * n, "rxn_subsystems": [""] * n, "rxn_ec": [""] * n}
    try:
        from metabolomics_modules import rewiring_reactions as rr
        links = rr.reaction_library().links([g or None for g in prep.gem_ids], max_steps=max_steps) if any(prep.gem_ids) else {}
    except Exception:
        links = {}
    if links:
        pos = {(int(i), int(j)): k for k, (i, j) in enumerate(zip(out["i"], out["j"]))}
        for (i, j), L in links.items():
            k = pos.get((i, j))
            if k is None:
                continue
            cols["rxn_steps"][k] = L["steps"]
            cols["rxn_dir"][k] = {"i>j": "a>b", "j>i": "b>a", "both": "both"}[L["direction"]]
            cols["rxn_via"][k] = L["via_name"] or ""
            cols["rxn_names"][k] = L["reaction_names"]
            cols["rxn_ids"][k] = L["reactions"]
            cols["rxn_subsystems"][k] = L["subsystems"]
            cols["rxn_ec"][k] = L["ec"]
    for c, v in cols.items():
        out[c] = v
    return out


def tested_mask(pairs: pd.DataFrame, scope: str, max_steps: int = 2) -> np.ndarray:
    steps = pairs["rxn_steps"].values
    if scope == "direct":
        return steps == 1
    if scope == "reaction":
        return (steps >= 1) & (steps <= max_steps)
    return np.ones(len(pairs), bool)


def classify_pairs(pairs: pd.DataFrame, cfg: RewiringConfig) -> pd.DataFrame:
    out = pairs.copy()
    metric = "p_used" if (cfg.cutoff_on == "p" and "p_used" in out) else ("p_value" if cfg.cutoff_on == "p" else "q_value")
    ra, rb, dr, q = out["r_a"].values, out["r_b"].values, out["delta_r"].values, out[metric].values
    tested = out["tested"].values.astype(bool) if "tested" in out else np.ones(len(out), bool)
    # "<=" so that a cutoff of 1 keeps every tested pair that passes |delta r|
    sig = tested & np.isfinite(q) & (np.nan_to_num(q, nan=2.0) <= cfg.fdr) & (np.abs(dr) >= cfg.min_delta)
    flipped = sig & (np.sign(ra) != np.sign(rb)) & (np.abs(ra) >= cfg.flip_min_abs_r) & (np.abs(rb) >= cfg.flip_min_abs_r)
    lost = sig & ~flipped & (np.abs(ra) > np.abs(rb))
    gained = sig & ~flipped & ~lost
    intact = ~sig & (np.abs(ra) >= cfg.intact_min_abs_r) & (np.abs(rb) >= cfg.intact_min_abs_r) & (np.sign(ra) == np.sign(rb))
    cls = np.full(len(out), None, dtype=object)
    cls[intact] = "intact"
    cls[lost] = "lost"
    cls[gained] = "gained"
    cls[flipped] = "flipped"
    out["class"] = cls
    return out


def feature_statistics(prep: PreparedData) -> pd.DataFrame:
    A, B = prep.X[:, prep.idx_a], prep.X[:, prep.idx_b]
    t, p = stats.ttest_ind(B, A, axis=1, equal_var=False)
    p = np.nan_to_num(p, nan=1.0)
    q = multipletests(p, method="fdr_bh")[1]
    return pd.DataFrame({
        "feature": prep.features, "pathway": prep.pathways,
        "mean_a": A.mean(axis=1), "mean_b": B.mean(axis=1),
        "log2fc": B.mean(axis=1) - A.mean(axis=1),
        "t": t, "p_value": p, "q_value": q,
    })


# ---------------------------------------------------------------------------
# 5. Global permutation test
# ---------------------------------------------------------------------------
def _pair_z(X, idx_a, idx_b, method, se, iu, ju) -> np.ndarray:
    CA = correlation_matrix(X, idx_a, method)
    CB = correlation_matrix(X, idx_b, method)
    dz = _fisher_z(CB[iu, ju]) - _fisher_z(CA[iu, ju])
    return dz / se if se else dz


def _group_centered(prep: PreparedData) -> np.ndarray:
    """Remove each group's mean from every feature. Permuting labels on uncentred data
    would mix two group means inside each permuted group, and a mean difference alone
    induces correlation between metabolites that both change -- inflating the null."""
    Xc = prep.X.copy()
    for g in (0, 1):
        idx = prep.group == g
        Xc[:, idx] -= Xc[:, idx].mean(axis=1, keepdims=True)
    return Xc


def exact_p_floor(n_a: int, n_b: int) -> float:
    """Smallest p-value an exact permutation test can give. The statistics used here (z^2) are
    unchanged when the two labels are swapped, so with equal group sizes every split and its
    mirror image tie: the floor is 2 / C(n, n_a) (3 + 3 samples: 2/20 = 0.10)."""
    total = math.comb(n_a + n_b, n_a)
    return (2 if n_a == n_b else 1) / total


def label_arrangements(prep: PreparedData, n_perm: int, seed: int = 0):
    """Group-label arrangements for the permutation null: every distinct split when there are
    no more than n_perm of them and no batch strata (exact test), otherwise n_perm random splits
    (within batch when a batch column is set). Returns (list of group-0 index arrays, exact)."""
    import itertools
    n = len(prep.group)
    n_a = prep.n_a
    total = math.comb(n, n_a)
    if prep.batch is None and total <= n_perm:
        return [np.array(c) for c in itertools.combinations(range(n), n_a)], True
    rng = np.random.default_rng(seed)
    strata = prep.batch if prep.batch is not None else np.zeros(n, dtype=object)
    out = []
    for _ in range(n_perm):
        g = prep.group.copy()
        for st in pd.unique(strata):
            idx = np.flatnonzero(strata == st)
            g[idx] = rng.permutation(g[idx])
        out.append(np.flatnonzero(g == 0))
    return out, False


def permutation_test(prep: PreparedData, method: str = "pearson", n_perm: int = 200, seed: int = 0,
                     keep_null: bool = False, tested: Optional[np.ndarray] = None,
                     sets: Optional[dict] = None) -> dict:
    """Label-permutation null for the whole analysis.

    Group labels are permuted (within batch when a batch column is set) on group-centred
    data, and every pair's statistic is recomputed each time. This gives
      * a global test: mean z^2 over the tested pairs; p = (1 + #null >= observed) / (1 + B),
        or, when every distinct label split is enumerated (small groups), the exact p =
        #splits with statistic >= observed / #splits (the observed split is one of them);
      * a pooled null distribution of |z| over the tested pairs, for the permutation FDR;
      * a pathway-level (set) test: for every set, mean z^2 over the pairs inside it,
        against the same permutations -- the differential-correlation analogue of MSEA.
    z is the Fisher-z difference scaled by its SE, or unscaled when a group has <= 3
    samples; a constant scale does not change any permutation p-value."""
    se = pair_scale(prep.n_a, prep.n_b, method)
    m = prep.X.shape[0]
    iu, ju = np.triu_indices(m, k=1)
    tested = np.ones(len(iu), bool) if tested is None else np.asarray(tested, bool)
    Xc = _group_centered(prep)
    # observed statistic from the same group-centred matrix the null uses: identical in exact
    # arithmetic (correlation ignores a group's mean), and bit-identical in floating point, so
    # the observed split is always counted as one of the enumerated splits (with |r| -> 1 in
    # tiny groups, Fisher z would otherwise amplify rounding differences)
    z_obs = _pair_z(Xc, prep.idx_a, prep.idx_b, method, se, iu, ju)
    obs = float(np.mean(z_obs[tested] ** 2)) if tested.any() else float("nan")
    # pair positions inside each set
    set_pos = {}
    if sets:
        pos = np.full((m, m), -1, dtype=np.int64)
        pos[iu, ju] = np.arange(len(iu))
        for name, idx in sets.items():
            idx = np.asarray(sorted(idx))
            if len(idx) >= 2:
                a, b = np.triu_indices(len(idx), k=1)
                set_pos[name] = pos[idx[a], idx[b]]
    set_obs = {k: float(np.mean(z_obs[v] ** 2)) for k, v in set_pos.items()}
    arrangements, exact = label_arrangements(prep, n_perm, seed)
    B = len(arrangements)
    all_idx = np.arange(len(prep.group))
    null = np.empty(B)
    set_ge = {k: 0 for k in set_pos}
    pooled = np.empty((B, int(tested.sum())), dtype=np.float32) if keep_null else None
    for k, ia in enumerate(arrangements):
        ib = np.setdiff1d(all_idx, ia, assume_unique=True)
        zk = _pair_z(Xc, ia, ib, method, se, iu, ju)
        null[k] = float(np.mean(zk[tested] ** 2)) if tested.any() else 0.0
        for name, v in set_pos.items():
            if float(np.mean(zk[v] ** 2)) >= set_obs[name] * (1 - 1e-9):
                set_ge[name] += 1
        if keep_null:
            pooled[k] = np.abs(zk[tested])
    ge = int(np.sum(null >= obs * (1 - 1e-9)))
    p = ge / B if exact else (1 + ge) / (1 + B)
    out = {
        "p_value": float(p), "observed": obs, "null_mean": float(null.mean()),
        "null_sd": float(null.std(ddof=1)) if B > 1 else 0.0,
        "n_permutations": int(B), "exact": bool(exact), "p_floor": float(exact_p_floor(prep.n_a, prep.n_b) if exact else 1 / (B + 1)),
        "stratified_by_batch": prep.batch is not None, "n_tested_pairs": int(tested.sum()),
    }
    out["_set_p"] = {k: (set_ge[k] / B if exact else (1 + set_ge[k]) / (1 + B)) for k in set_pos}
    out["_set_obs"] = set_obs
    if keep_null:
        out["_null_abs_z"] = np.sort(pooled.ravel())
    return out


def permutation_fdr(abs_z: np.ndarray, null_sorted: np.ndarray, n_perm: int) -> np.ndarray:
    """SAM-style permutation q-values. For a threshold t, FDR(t) = (expected number of
    null pairs with |z| >= t, averaged over permutations) / (observed pairs with |z| >= t),
    with pi0 = 1 (conservative) and the null count floored at 1. q_i = min over
    thresholds that include pair i."""
    abs_z = np.asarray(abs_z, dtype=float)
    n = len(abs_z)
    if n == 0:
        return abs_z
    null_ge = len(null_sorted) - np.searchsorted(null_sorted, abs_z.astype(np.float32), side="left")
    obs_sorted = np.sort(abs_z)
    obs_ge = n - np.searchsorted(obs_sorted, abs_z, side="left")
    # at least one null exceedance is assumed, so q is never reported as exactly 0 (the
    # permutation resolution is 1 / n_perm)
    fdr = np.minimum(1.0, (np.maximum(null_ge, 1) / n_perm) / np.maximum(obs_ge, 1))
    # rejection region {|z| >= t'} contains pair i iff t' <= |z_i|; in ascending |z| order
    # that is every index up to i, so q is a running minimum from the smallest |z| upward
    order = np.argsort(abs_z, kind="stable")
    q = np.empty(n)
    q[order] = np.minimum.accumulate(fdr[order])
    return np.minimum(q, 1.0)


# ---------------------------------------------------------------------------
# 6. Hypothesis cards
# ---------------------------------------------------------------------------
_TAG_PATTERNS = [
    ("tca", r"\btca\b|citr(ate|ic acid) cycle|krebs|tricarboxylic"),
    ("tryptophan", r"tryptophan|kynuren"),
    ("glutamine", r"glutamin|glutamat|alanine,? aspartate|aspartate"),
    ("glycolysis", r"glycoly|gluconeo|pyruvate metab|warburg"),
    ("fao", r"fatty acid|carnitine|beta.?oxid|acylcarn"),
    ("purine", r"purine"),
    ("pyrimidine", r"pyrimidine"),
    ("onecarbon", r"one.?carbon|folate|methionine|serine|glycine"),
    ("redox", r"glutathione|redox|oxidative"),
    ("bile", r"bile"),
    ("lipid", r"sphingo|ceramide|glycerophospho|phospholipid|glycerolipid|lipid"),
    ("urea", r"urea|arginine|proline|polyamine"),
    ("ppp", r"pentose"),
]
_TRACERS = {
    "tca": "U-¹³C₆-glucose and U-¹³C₅-glutamine", "glutamine": "U-¹³C₅-glutamine",
    "glycolysis": "U-¹³C₆-glucose", "tryptophan": "¹³C₁₁-tryptophan",
    "fao": "U-¹³C₁₆-palmitate", "onecarbon": "U-¹³C₃-serine", "purine": "¹⁵N-amide-glutamine or ¹³C₂-glycine",
    "pyrimidine": "¹⁵N-amide-glutamine", "redox": "¹³C₂,¹⁵N-glycine (glutathione synthesis)",
    "urea": "¹⁵N₂-arginine", "ppp": "1,2-¹³C₂-glucose", "lipid": "U-¹³C₆-glucose (lipogenesis) or labelled fatty acids",
    "bile": "deuterated cholic acid",
}

# (sorted tag pair, class) -> (interpretation, experiment)
_INTERPRETATIONS = {
    (("glutamine", "tca"), "lost"): (
        "Glutamine/glutamate no longer track TCA intermediates. This is consistent with anaplerotic "
        "decoupling: glutamine carbon may be diverted to other fates (nucleotides, glutathione, secretion) "
        "while the cycle is fed from another source.",
        "U-¹³C₅-glutamine tracing; compare M+4 citrate and M+5 α-ketoglutarate enrichment between groups."),
    (("tryptophan", "tryptophan"), "flipped"): (
        "Kynurenine-branch metabolites move opposite to tryptophan, so their levels look set by catabolism "
        "rather than by tryptophan supply. That pattern is consistent with active IDO1/TDO2, an "
        "immunosuppressive axis reported in several tumor types.",
        "IDO1/TDO2 staining on the same specimens and ¹³C₁₁-tryptophan tracing with and without an IDO1 inhibitor."),
    (("tryptophan", "tryptophan"), "lost"): (
        "Branches of tryptophan metabolism stop co-varying, consistent with tryptophan being routed "
        "preferentially into one branch (for example kynurenine over serotonin).",
        "Measure TPH1 versus IDO1/TDO2 expression; ¹³C₁₁-tryptophan tracing reports the branch split directly."),
    (("glycolysis", "tryptophan"), "gained"): (
        "Glycolytic end-products and kynurenine-pathway metabolites co-vary only in the comparison group, "
        "suggesting a shared upstream driver such as hypoxia or inflammatory signalling.",
        "U-¹³C₆-glucose tracing with co-measured kynurenine; test whether lactate or hypoxia induces IDO1 in co-culture."),
    (("glycolysis", "glycolysis"), "lost"): (
        "Glycolytic intermediates partly decouple, meaning at least one of them (often lactate) is now "
        "controlled by something beyond glycolytic flux, such as export or a separate source.",
        "U-¹³C₆-glucose tracing; compare M+3 lactate with M+3 pyruvate fractions."),
    (("tca", "tca"), "lost"): (
        "TCA intermediates stop co-varying, which fits a truncated or split cycle (for example an SDH or FH "
        "block, or separate carbon supply to each half of the cycle).",
        "U-¹³C₆-glucose and U-¹³C₅-glutamine tracing; look for M+4 vs M+3 patterns around succinate/fumarate."),
    (("glycolysis", "tca"), "lost"): (
        "Glycolytic and TCA metabolites decouple, consistent with reduced pyruvate entry into mitochondria "
        "(PDH inhibition or MPC loss) and a shift towards lactate.",
        "U-¹³C₆-glucose tracing; compare M+2 citrate with M+3 lactate."),
    (("fao", "fao"), "lost"): (
        "Acylcarnitine species stop moving together, consistent with a change in which chain lengths are "
        "oxidised or exported rather than a uniform change in β-oxidation.",
        "U-¹³C₁₆-palmitate tracing into acetyl-carnitine and citrate; CPT1/CPT2 activity assay."),
    (("onecarbon", "redox"), "gained"): (
        "One-carbon and glutathione metabolites become coupled, consistent with serine/glycine being "
        "recruited to support glutathione synthesis under oxidative stress.",
        "U-¹³C₃-serine tracing into glycine and glutathione."),
}
_GENERIC = {
    "lost": ("A coordinated relationship present in the reference group breaks down. Common causes are a "
             "bypassed regulatory node, a new competing consumer, or loss of shared transcriptional control."),
    "gained": ("A coupling appears that the reference group lacks, which often means a shared upstream driver "
               "has switched on."),
    "flipped": ("The relationship reverses direction, typical of a substrate–product pair whose controlling "
                "step has changed."),
}
_FRAGILE_TEXT = ("Treat as a probable artifact. The signal does not survive resampling and/or a change of "
                 "correlation method, which usually means a few extreme samples are creating it.")
_OUTLIER_TEXT = ("Possibly driven by a few extreme samples: the effect largely disappears under the other "
                 "correlation method. Check the scatter plot before believing it.")
_IMPUTE_TEXT = ("Likely an imputation artifact: the coupling disappears when samples with imputed (not-detected) "
                "values are left out, so it mostly reflects which samples fell below the detection limit.")
_IMPUTE_NEXT = ("Check detection limits and missingness patterns for these metabolites; re-run with a stricter "
                "missingness filter or confirm with targeted quantitation.")
_OUTLIER_NEXT = ("Inspect the outlying samples in the scatter plot, check their QC (injection order, internal "
                 "standards, dilution), and re-run with Spearman correlation.")


def pathway_tag(name: str) -> Optional[str]:
    text = str(name).lower()
    for tag, pat in _TAG_PATTERNS:
        if re.search(pat, text):
            return tag
    return None


def _interpret(pw_a: str, pw_b: str, cls: str):
    ta, tb = pathway_tag(pw_a), pathway_tag(pw_b)
    if ta and tb:
        hit = _INTERPRETATIONS.get((tuple(sorted((ta, tb))), cls))
        if hit:
            return hit[0], hit[1], True
    tracers = [t for t in dict.fromkeys([_TRACERS.get(ta), _TRACERS.get(tb)]) if t]
    if tracers:
        exp = (f"Stable-isotope tracing with {' plus '.join(tracers)} to test whether flux through the connecting "
               "reactions differs; replicate the correlation change in an independent cohort.")
    else:
        exp = ("Replicate in an independent cohort with targeted quantitation, then trace the shared precursor "
               "with a stable-isotope label.")
    return _GENERIC[cls], exp, False


def _pair_stat(X, pi, pj, idx_a, idx_b, method, sign) -> float:
    """Mean over pairs of (r_B - r_A) * expected direction, for one resample."""
    def rs(idx):
        rows = np.unique(np.concatenate([pi, pj]))
        Z = _standardize(X[np.ix_(rows, idx)], method)
        pos = {r: k for k, r in enumerate(rows)}
        a = Z[[pos[r] for r in pi]]
        b = Z[[pos[r] for r in pj]]
        return np.clip((a * b).sum(axis=1) / (len(idx) - 1), -1, 1)
    return float(np.mean((rs(idx_b) - rs(idx_a)) * sign))


def _influence_check(prep, pi, pj, method, sign, observed, frac) -> dict:
    """Outlier-sample check. Within each group, a sample is outlying for this module if
    any module metabolite sits more than 4 robust SDs (median/MAD) from the group
    median. Those samples (at most `frac` of a group) are removed and the module
    statistic recomputed. A real module survives; one created by a handful of extreme
    samples collapses. (Leave-one-out influence is not used because several outliers
    mask each other.)"""
    rows = np.unique(np.concatenate([pi, pj]))
    drop = []
    for idx in (prep.idx_a, prep.idx_b):
        V = prep.X[np.ix_(rows, idx)]
        med = np.median(V, axis=1, keepdims=True)
        mad = 1.4826 * np.median(np.abs(V - med), axis=1, keepdims=True)
        mad[mad == 0] = np.inf
        rz = np.abs(V - med) / mad
        worst = rz.max(axis=0)
        cand = [(w, s) for w, s in zip(worst, idx) if w > 4.0]
        cand.sort(reverse=True)
        drop += [int(s) for _, s in cand[: max(1, int(frac * len(idx)))]]
    if not drop:
        return {"ok": True, "value": observed, "observed": observed, "dropped": []}
    keep_a = np.array([s for s in prep.idx_a if int(s) not in set(drop)])
    keep_b = np.array([s for s in prep.idx_b if int(s) not in set(drop)])
    value = _pair_stat(prep.X, pi, pj, keep_a, keep_b, method, sign)
    return {
        "ok": bool(value >= 0.5 * observed and value >= 0.2), "value": value, "observed": observed,
        "dropped": [prep.samples[s] for s in drop],
    }


def _complete_case_stat(prep, pi, pj, method, sign, min_n=5):
    """Module statistic using, for each pair, only samples where neither metabolite was
    imputed. Returns (value, n_pairs_used)."""
    imp = prep.imputed
    vals = []
    for a, b, sg in zip(pi, pj, sign):
        ok = ~(imp[a] | imp[b])
        ia = np.flatnonzero(ok & (prep.group == 0))
        ib = np.flatnonzero(ok & (prep.group == 1))
        if len(ia) < min_n or len(ib) < min_n:
            vals.append(0.0)       # can't be supported without imputed values
            continue
        za = _standardize(prep.X[np.ix_([a, b], ia)], method)
        zb = _standardize(prep.X[np.ix_([a, b], ib)], method)
        r_a = float(np.clip((za[0] * za[1]).sum() / (len(ia) - 1), -1, 1))
        r_b = float(np.clip((zb[0] * zb[1]).sum() / (len(ib) - 1), -1, 1))
        vals.append((r_b - r_a) * sg)
    return float(np.mean(vals)) if vals else 0.0


def _imputation_check(prep, pi, pj, method, sign, observed) -> dict:
    """Imputed values sit together at the bottom of the range, so two metabolites that
    are both 'not detected' in the same samples can look strongly correlated. Recompute
    the module from complete cases only."""
    if prep.imputed is None:
        return {"ok": True, "value": observed, "observed": observed, "n_imputed": 0}
    rows = np.unique(np.concatenate([pi, pj]))
    n_imp = int(prep.imputed[rows].sum())
    if n_imp == 0:
        return {"ok": True, "value": observed, "observed": observed, "n_imputed": 0}
    value = _complete_case_stat(prep, pi, pj, method, sign, min_n=min(5, prep.n_a, prep.n_b))
    return {"ok": bool(value >= 0.5 * observed and value >= 0.2), "value": value, "observed": observed,
            "n_imputed": n_imp}


SMALL_N_ROBUSTNESS = 5   # below this many samples in a group, resampling checks are not assessable


def _reaction_text(r) -> str:
    a, b = r["feature_a"], r["feature_b"]
    arrow = {"a>b": f"{a} → {b}", "b>a": f"{b} → {a}", "both": f"{a} ⇄ {b}"}.get(r["rxn_dir"], f"{a} – {b}")
    return arrow + (f" (via {r['rxn_via']})" if r["rxn_steps"] == 2 and r["rxn_via"] else "")


def build_hypotheses(prep: PreparedData, pairs: pd.DataFrame, feats: pd.DataFrame, cfg: RewiringConfig,
                     ptab: Optional[pd.DataFrame] = None) -> list:
    rew = pairs[pairs["class"].isin(REWIRED)]
    if rew.empty:
        return []
    ga, gb = prep.group_labels
    rng = np.random.default_rng(cfg.seed + 17)
    other = "spearman" if cfg.method == "pearson" else "pearson"
    fq = feats["q_value"].values
    small_n = min(prep.n_a, prep.n_b) < SMALL_N_ROBUSTNESS
    set_q = dict(zip(ptab["pathway"], ptab["q_value"])) if ptab is not None and not ptab.empty else {}
    if "rxn_steps" not in rew:
        rew = rew.assign(rxn_steps=0, rxn_dir="", rxn_via="")

    modules = []
    for (pa, pb, cls), grp in rew.assign(
        pw_lo=np.where(rew["pathway_a"] <= rew["pathway_b"], rew["pathway_a"], rew["pathway_b"]),
        pw_hi=np.where(rew["pathway_a"] <= rew["pathway_b"], rew["pathway_b"], rew["pathway_a"]),
        _rx=(rew["rxn_steps"] > 0).astype(int), _adr=rew["delta_r"].abs(),
    ).groupby(["pw_lo", "pw_hi", "class"], sort=False):
        # reaction-linked pairs first (they carry the mechanism), then by effect size
        grp = grp.sort_values(["_rx", "_adr"], ascending=[False, False]).head(cfg.max_pairs_per_module)
        n = len(grp)
        pre = float(grp["delta_r"].abs().mean()) * math.log2(1 + n)
        modules.append({"pw_a": pa, "pw_b": pb, "cls": cls, "pairs": grp, "pre": pre})
    modules.sort(key=lambda m: m["pre"], reverse=True)
    modules = modules[: cfg.max_modules_scored]

    cards = []
    for M in modules:
        grp = M["pairs"]
        pi, pj = grp["i"].values, grp["j"].values
        sign = np.sign(grp["delta_r"].values)
        observed = _pair_stat(prep.X, pi, pj, prep.idx_a, prep.idx_b, cfg.method, sign)

        if small_n:
            stability = None
            influence_check = {"ok": None, "value": observed, "observed": observed, "dropped": [],
                               "reason": f"Not assessable with fewer than {SMALL_N_ROBUSTNESS} samples per group."}
        else:
            keep = 0
            ia, ib = prep.idx_a, prep.idx_b
            for _ in range(cfg.n_bootstrap):
                ra = rng.choice(ia, size=len(ia), replace=True)
                rb = rng.choice(ib, size=len(ib), replace=True)
                if _pair_stat(prep.X, pi, pj, ra, rb, cfg.method, sign) >= cfg.robust_delta:
                    keep += 1
            stability = keep / cfg.n_bootstrap if cfg.n_bootstrap else None
            influence_check = _influence_check(prep, pi, pj, cfg.method, sign, observed, cfg.influence_frac)
        imputation_check = _imputation_check(prep, pi, pj, cfg.method, sign, observed)

        other_val = _pair_stat(prep.X, pi, pj, prep.idx_a, prep.idx_b, other, sign)
        method_ok = other_val >= 0.5 * observed
        method_check = {"method": other, "ok": bool(method_ok), "value": other_val, "observed": observed}

        batch_check = {"ok": None, "levels": [], "reason": "No batch column selected."}
        if prep.batch is not None:
            levels = []
            for lv in pd.unique(prep.batch):
                a = np.flatnonzero((prep.group == 0) & (prep.batch == lv))
                b = np.flatnonzero((prep.group == 1) & (prep.batch == lv))
                if len(a) >= cfg.batch_min_per_group and len(b) >= cfg.batch_min_per_group:
                    levels.append({"level": str(lv), "value": _pair_stat(prep.X, pi, pj, a, b, cfg.method, sign),
                                   "n_a": int(len(a)), "n_b": int(len(b))})
            if len(levels) >= 2:
                batch_check = {"ok": bool(all(l["value"] >= 0.2 for l in levels)), "levels": levels, "reason": ""}
            else:
                batch_check = {"ok": None, "levels": levels,
                               "reason": f"Fewer than two batches have ≥{cfg.batch_min_per_group} samples per group."}

        feat_idx = sorted(set(pi.tolist()) | set(pj.tolist()))
        hidden = [prep.features[k] for k in feat_idx if fq[k] >= 0.05]
        stab_fail = stability is not None and stability < cfg.stability_cutoff
        fails = (int(stab_fail) + int(not method_ok) + int(influence_check["ok"] is False)
                 + int(not imputation_check["ok"]) + int(batch_check["ok"] is False))
        if small_n:
            verdict = "fragile" if fails >= 2 else "exploratory"
        else:
            verdict = "robust" if fails == 0 else ("fragile" if fails >= 2 or stability < 0.4 else "mixed")
        mean_abs = float(grp["delta_r"].abs().mean())
        stab_w = stability if stability is not None else 0.5
        score = mean_abs * math.log2(1 + len(grp)) * (0.3 + stab_w) * (0.5 ** fails)

        rx = grp[grp["rxn_steps"] > 0]
        reaction_pairs = [{"text": _reaction_text(r), "steps": int(r["rxn_steps"]), "reactions": r.get("rxn_names", ""),
                           "ec": r.get("rxn_ec", "")} for _, r in rx.iterrows()]
        score *= 1 + 0.25 * min(1.0, len(rx) / max(1, len(grp)))   # mechanism-backed modules rank higher

        pa, pb, cls = M["pw_a"], M["pw_b"], M["cls"]
        same = pa == pb
        where = f"within {pa}" if same else f"between {pa} and {pb}"
        r_a, r_b = float(grp["r_a"].mean()), float(grp["r_b"].mean())
        n = len(grp)
        plural = "s" if n > 1 else ""
        if cls == "lost":
            claim = f"Coupling {where} breaks down in {gb}: mean r falls from {r_a:.2f} to {r_b:.2f} across {n} pair{plural}."
        elif cls == "gained":
            claim = f"{gb} gains a coupling {where} that {ga} lacks: mean r {r_a:.2f} → {r_b:.2f} across {n} pair{plural}."
        else:
            claim = f"The relationship {where} reverses sign: mean r {r_a:.2f} → {r_b:.2f} across {n} pair{plural}."
        interp, experiment, curated = _interpret(pa, pb, cls)
        if verdict not in ("robust", "exploratory") and not imputation_check["ok"]:
            interp, experiment = _IMPUTE_TEXT, _IMPUTE_NEXT
        elif verdict not in ("robust", "exploratory") and (not method_ok or influence_check["ok"] is False):
            interp, experiment = _OUTLIER_TEXT, _OUTLIER_NEXT
        if verdict == "fragile":
            interp, experiment = _FRAGILE_TEXT, _OUTLIER_NEXT

        cards.append({
            "key": f"{pa}|{pb}|{cls}", "pathway_a": pa, "pathway_b": pb, "class": cls,
            "n_pairs": n, "r_a": r_a, "r_b": r_b, "mean_abs_delta": mean_abs,
            "stability": stability, "method_check": method_check, "influence_check": influence_check,
            "imputation_check": imputation_check,
            "batch_check": batch_check,
            "features": [prep.features[k] for k in feat_idx], "feature_idx": feat_idx,
            "hidden_features": hidden,
            "pairs": [[int(a), int(b)] for a, b in zip(pi, pj)],
            "pair_names": [[prep.features[a], prep.features[b]] for a, b in zip(pi, pj)],
            "reaction_pairs": reaction_pairs,
            "pathway_q": {k: set_q[k] for k in dict.fromkeys([pa, pb]) if k in set_q},
            "claim": claim, "interpretation": interp, "next_experiment": experiment,
            "curated_interpretation": curated, "verdict": verdict, "score": score,
        })
    cards.sort(key=lambda c: c["score"], reverse=True)
    for k, c in enumerate(cards):
        c["rank"] = k + 1
    return cards


# ---------------------------------------------------------------------------
# Orchestration + serialization
# ---------------------------------------------------------------------------
@dataclass
class RewiringResult:
    prep: PreparedData
    cfg: RewiringConfig
    pairs: pd.DataFrame
    features: pd.DataFrame
    hypotheses: list
    permutation: Optional[dict]
    fdr_used: str = "bh"
    pathway_tests: pd.DataFrame = field(default_factory=pd.DataFrame)
    reactions: pd.DataFrame = field(default_factory=pd.DataFrame)

    def counts(self) -> dict:
        vc = self.pairs["class"].value_counts()
        return {c: int(vc.get(c, 0)) for c in CLASSES}

    def rewired_pairs(self) -> pd.DataFrame:
        cols = (["feature_a", "feature_b", "pathway_a", "pathway_b", "class", "r_a", "r_b", "delta_r", "z", "p_value",
                 "q_bh"] + [c for c in ("p_perm", "q_perm") if c in self.pairs] + ["q_value"])
        rx = [c for c in ("rxn_steps", "rxn_dir", "rxn_via", "rxn_names", "rxn_ec") if c in self.pairs]
        out = self.pairs[self.pairs["class"].isin(REWIRED)][cols + rx].copy()
        ga, gb = self.prep.group_labels
        out = out.rename(columns={"r_a": f"r_{ga}", "r_b": f"r_{gb}", "rxn_steps": "reaction_steps",
                                  "rxn_dir": "reaction_direction", "rxn_via": "via", "rxn_names": "reactions",
                                  "rxn_ec": "EC"})
        return out.sort_values("q_value").reset_index(drop=True)

    def hidden_feature_count(self) -> int:
        rew = self.pairs[self.pairs["class"].isin(REWIRED)]
        idx = set(rew["i"]).union(rew["j"])
        q = self.features["q_value"].values
        return int(sum(1 for k in idx if q[k] >= 0.05))

    def pathway_table_display(self) -> pd.DataFrame:
        if self.pathway_tests is None or self.pathway_tests.empty:
            return pd.DataFrame()
        return self.pathway_tests.drop(columns=["feature_idx"])

    def reaction_table_display(self) -> pd.DataFrame:
        if self.reactions is None or self.reactions.empty:
            return pd.DataFrame()
        ga, gb = self.prep.group_labels
        return self.reactions.drop(columns=["i", "j"]).rename(columns={
            "r_a": f"r_{ga}", "r_b": f"r_{gb}", "log2_ratio_a": f"log2(product/substrate) {ga}",
            "log2_ratio_b": f"log2(product/substrate) {gb}", "ratio_shift": f"ratio shift ({gb} - {ga})"})

    def summary(self) -> dict:
        P = self.prep
        rx = self.pairs["rxn_steps"] if "rxn_steps" in self.pairs else pd.Series(0, index=self.pairs.index)
        rew = self.pairs["class"].isin(REWIRED)
        pt = self.pathway_tests if self.pathway_tests is not None else pd.DataFrame()
        idm = P.id_map
        return {
            "n_features": len(P.features), "n_a": P.n_a, "n_b": P.n_b,
            "groups": list(P.group_labels), "method": self.cfg.method,
            "fdr": self.cfg.fdr, "cutoff_on": self.cfg.cutoff_on, "fdr_method": self.fdr_used,
            "min_delta": self.cfg.min_delta, "counts": self.counts(),
            "hidden_rewired_features": self.hidden_feature_count(),
            "permutation": self.permutation, "pathway_source": P.pathway_source,
            "pathway_library": P.pathway_library,
            "pair_scope": self.cfg.pair_scope, "n_tested_pairs": int(self.pairs["tested"].sum()) if "tested" in self.pairs else len(self.pairs),
            "n_pairs_total": len(self.pairs),
            "reaction_links": {"direct": int((rx == 1).sum()), "two_step": int((rx == 2).sum()),
                               "rewired_with_reaction": int((rew & (rx > 0)).sum())},
            "n_pathways_tested": int(len(pt)), "n_pathways_significant": int((pt["q_value"] < 0.05).sum()) if len(pt) else 0,
            "n_id_mapped": int((idm["Match_Method"] != "Unmapped").sum()) if idm is not None else None,
            "n_in_reaction_network": int((idm["GEM_ID"] != "").sum()) if idm is not None else None,
            "notes": list(P.notes),
        }

    def to_payload(self, max_intact: int = 400, max_cards: int = 6, max_reactions: int = 800) -> dict:
        """JSON-safe dict for the interactive view (no NaN/inf)."""
        P = self.prep
        pairs = self.pairs
        if "rxn_steps" not in pairs:
            pairs = pairs.assign(rxn_steps=0, rxn_dir="", rxn_via="", rxn_names="", rxn_ec="")
        cls = pairs[pairs["class"].notna()]
        intact = cls[cls["class"] == "intact"]
        if len(intact) > max_intact:
            keep = intact.assign(s=np.minimum(intact["r_a"].abs(), intact["r_b"].abs())).nlargest(max_intact, "s").index
            cls = pd.concat([cls[cls["class"] != "intact"], intact.loc[keep]])
        rew = cls[cls["class"] != "intact"]
        rxl = pairs[pairs["rxn_steps"] > 0]
        rxl = rxl.reindex(rxl["rxn_steps"].sort_values(kind="stable").index).head(max_reactions)
        ratio = {}
        if self.reactions is not None and not self.reactions.empty:
            for r in self.reactions.itertuples():
                ratio[(min(r.i, r.j), max(r.i, r.j))] = (r.i, r.j, r.ratio_shift, r.ratio_q, r.log2_ratio_a, r.log2_ratio_b)
        need = sorted(set(rew["i"]).union(rew["j"]).union(rxl["i"]).union(rxl["j"]))
        if len(P.features) * len(P.samples) <= 200_000:
            need = list(range(len(P.features)))
        cards = self.hypotheses[:max_cards]
        flagged = next((c for c in self.hypotheses if c["verdict"] in ("mixed", "fragile")), None)
        if flagged is not None and flagged not in cards and len(cards) == max_cards:
            cards = cards[: max_cards - 1] + [flagged]
        idm = P.id_map

        def pair_rec(r):
            d = {"i": int(r.i), "j": int(r.j), "rA": float(r.r_a), "rB": float(r.r_b), "dr": float(r.delta_r),
                 "z": float(r.z), "q": float(r.q_value), "p": float(getattr(r, "p_used", np.nan)),
                 "cls": r.cls, "rx": int(r.rxn_steps)}
            if r.rxn_steps:
                d.update({"rxd": r.rxn_dir, "rxv": r.rxn_via, "rxn": str(r.rxn_names)[:300], "ec": r.rxn_ec})
                rt = ratio.get((int(r.i), int(r.j)))
                if rt:
                    d.update({"sub": rt[0], "prod": rt[1], "rs": rt[2], "rsq": rt[3], "ra": rt[4], "rb": rt[5]})
            return d

        def rows(df):
            return [pair_rec(r) for r in df.rename(columns={"class": "cls"}).itertuples()]

        pt = self.pathway_tests if self.pathway_tests is not None else pd.DataFrame()
        top = pairs[pairs["tested"]] if "tested" in pairs else pairs
        top = top.dropna(subset=["z"])
        top = top.reindex(top["z"].abs().sort_values(ascending=False).index).head(40)
        payload = {
            "version": 2,
            "groups": list(P.group_labels), "method": self.cfg.method,
            "fdr": self.cfg.fdr, "cutoff_on": self.cfg.cutoff_on, "min_delta": self.cfg.min_delta,
            "scale_label": P.scale_label,
            "n": [P.n_a, P.n_b], "samples": P.samples, "group": P.group.tolist(),
            "features": [
                {"name": P.features[k], "pathway": P.pathways[k],
                 "lfc": float(self.features["log2fc"].iat[k]), "q": float(self.features["q_value"].iat[k]),
                 "gem": (idm["GEM_Name"].iat[k] if idm is not None else ""),
                 "hmdb": (idm["HMDB_ID"].iat[k] if idm is not None else ""),
                 "all": (idm["All_Pathways"].iat[k] if idm is not None and "All_Pathways" in idm else "")}
                for k in range(len(P.features))
            ],
            "pairs": rows(cls),
            "reactions": rows(rxl.assign(**{"class": rxl["class"]})),
            # largest correlation changes among tested pairs, for the exploratory list shown when no
            # pair passes the FDR threshold (typical with 3-5 samples per group)
            "top_pairs": rows(top),
            "pathways": [
                {"name": r.pathway, "n": int(r.n_measured), "np": int(r.n_pairs), "stat": float(r.statistic),
                 "dr": float(r.mean_abs_delta_r), "dir": r.direction, "p": float(r.p_value), "q": float(r.q_value),
                 "idx": list(r.feature_idx), "top": r.top_pairs, "nrx": int(r.n_reaction_pairs),
                 "nl": int(r.n_lost), "ng": int(r.n_gained), "nf": int(r.n_flipped)}
                for r in pt.head(80).itertuples()
            ] if len(pt) else [],
            "values": {str(k): [float(v) for v in P.X[k]] for k in need},
            "hypotheses": cards,
            "summary": self.summary(),
        }
        return _json_safe(payload)

    def evidence_json(self) -> str:
        """Everything behind the hypothesis cards, for audit or for handing to an LLM."""
        pt = self.pathway_table_display()
        rt = self.reaction_table_display()
        return json.dumps(_json_safe({
            "summary": self.summary(), "hypotheses": self.hypotheses,
            "pathway_level_tests": pt.head(50).to_dict("records") if len(pt) else [],
            "reaction_pairs": rt.head(100).to_dict("records") if len(rt) else [],
        }), indent=2, ensure_ascii=False, default=str)


def _json_safe(obj):
    if isinstance(obj, dict):
        return {str(k): _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating, float)):
        f = float(obj)
        return f if math.isfinite(f) else None
    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, np.ndarray):
        return _json_safe(obj.tolist())
    return obj


def pathway_table(prep: PreparedData, pairs: pd.DataFrame, perm: Optional[dict], cfg: RewiringConfig) -> pd.DataFrame:
    """Pathway-level (set) rewiring test results, one row per set with >= pathway_min_size
    measured members: mean z^2 of the pairs inside the set against label permutations,
    BH-adjusted across sets."""
    if not perm or "_set_p" not in perm:
        return pd.DataFrame()
    m = len(prep.features)
    pos = np.full((m, m), -1, dtype=np.int64)
    pos[pairs["i"].values, pairs["j"].values] = np.arange(len(pairs))
    rows = []
    for name, pval in perm["_set_p"].items():
        idx = np.asarray(sorted(prep.sets[name]))
        a, b = np.triu_indices(len(idx), k=1)
        sub = pairs.iloc[pos[idx[a], idx[b]]]
        ra, rb = sub["r_a"].abs(), sub["r_b"].abs()
        cls = sub["class"].value_counts()
        top = sub.reindex(sub["z"].abs().sort_values(ascending=False).index).head(3)
        rows.append({
            "pathway": name, "n_measured": int(len(idx)), "n_pairs": int(len(sub)),
            "members": "; ".join(prep.features[k] for k in idx),
            "statistic": perm["_set_obs"][name], "mean_abs_delta_r": float(sub["delta_r"].abs().mean()),
            "mean_abs_r_a": float(ra.mean()), "mean_abs_r_b": float(rb.mean()),
            "direction": "coupling lost" if rb.mean() < ra.mean() - 0.05 else
                         ("coupling gained" if rb.mean() > ra.mean() + 0.05 else "reorganised"),
            "n_lost": int(cls.get("lost", 0)), "n_gained": int(cls.get("gained", 0)),
            "n_flipped": int(cls.get("flipped", 0)),
            "n_reaction_pairs": int((sub["rxn_steps"] > 0).sum()) if "rxn_steps" in sub else 0,
            "top_pairs": "; ".join(f"{x}–{y}" for x, y in zip(top["feature_a"], top["feature_b"])),
            "p_value": float(pval), "feature_idx": [int(k) for k in idx],
        })
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df["q_value"] = multipletests(df["p_value"].values, method="fdr_bh")[1]
    return df.sort_values(["p_value", "statistic"], ascending=[True, False]).reset_index(drop=True)


def reaction_table(prep: PreparedData, pairs: pd.DataFrame) -> pd.DataFrame:
    """Every measured pair linked by a reaction (direct or via one unmeasured intermediate):
    correlation in each group plus the product/substrate log2 ratio per group, its shift
    (Welch t-test, BH across linked pairs). For reversible links the orientation is as listed."""
    L = pairs[pairs.get("rxn_steps", pd.Series(0, index=pairs.index)) > 0]
    if L.empty:
        return pd.DataFrame()
    rows = []
    for r, cls in zip(L.itertuples(), L["class"].values):
        sub, prod = (r.j, r.i) if r.rxn_dir == "b>a" else (r.i, r.j)
        ratio = prep.X[prod] - prep.X[sub]
        a, b = ratio[prep.idx_a], ratio[prep.idx_b]
        t, pv = stats.ttest_ind(b, a, equal_var=False)
        rows.append({
            "substrate": prep.features[sub], "product": prep.features[prod], "i": int(sub), "j": int(prod),
            "reversible": r.rxn_dir == "both", "steps": int(r.rxn_steps), "via": r.rxn_via,
            "reactions": r.rxn_names, "reaction_ids": r.rxn_ids, "subsystems": r.rxn_subsystems, "ec": r.rxn_ec,
            "r_a": float(r.r_a), "r_b": float(r.r_b), "delta_r": float(r.delta_r), "pair_q": float(r.q_value),
            "class": cls,
            "log2_ratio_a": float(a.mean()), "log2_ratio_b": float(b.mean()),
            "ratio_shift": float(b.mean() - a.mean()), "ratio_p": float(pv) if np.isfinite(pv) else np.nan,
        })
    df = pd.DataFrame(rows)
    df["ratio_q"] = _bh(df["ratio_p"].values)
    return df.sort_values(["steps", "ratio_p"], na_position="last").reset_index(drop=True)


def run_rewiring(prep: PreparedData, cfg: Optional[RewiringConfig] = None, *, run_permutation: bool = True) -> RewiringResult:
    """Full analysis. With permutations (default) pairs are called on the permutation FDR;
    without, on Benjamini-Hochberg over the analytic Fisher-z p-values (not available when a
    group has <= 3 samples)."""
    cfg = cfg or RewiringConfig()
    pairs = differential_correlation(prep, cfg.method)
    pairs = annotate_reactions(prep, pairs, cfg.max_steps)
    tested = tested_mask(pairs, cfg.pair_scope, cfg.max_steps)
    if cfg.pair_scope != "all" and not tested.any():
        raise RewiringInputError(
            "No measured metabolite pair is linked by a reaction in the Human-GEM library, so the "
            "reaction-linked scope has nothing to test. Choose 'All metabolite pairs'.")
    pairs["tested"] = tested
    pairs["q_bh"] = _bh(pairs["p_value"].values, tested)
    pairs = pairs.drop(columns=["q_value"])
    use_perm = run_permutation and cfg.n_permutations > 0
    if not use_perm and pair_scale(prep.n_a, prep.n_b, cfg.method) is None:
        raise RewiringInputError("With 3 samples in a group the analytic Fisher-z test is undefined; "
                                 "set label permutations above 0.")
    eligible = {k: v for k, v in prep.sets.items() if len(v) >= cfg.pathway_min_size}
    perm = None
    if use_perm:
        perm = permutation_test(prep, cfg.method, cfg.n_permutations, cfg.seed,
                                keep_null=cfg.fdr_method == "permutation" or pair_scale(prep.n_a, prep.n_b, cfg.method) is None,
                                tested=tested, sets=eligible)
        null = perm.pop("_null_abs_z", None)
        if null is not None:
            qp = np.full(len(pairs), np.nan)
            qp[tested] = permutation_fdr(np.abs(pairs["z"].values[tested]), null, perm["n_permutations"])
            pairs["q_perm"] = qp
            # per-pair permutation p-value against the pooled null of all tested pairs
            pp = np.full(len(pairs), np.nan)
            az = np.abs(pairs["z"].values[tested]).astype(np.float32)
            ge = len(null) - np.searchsorted(null, az, side="left")
            pp[tested] = (1 + ge) / (1 + len(null))
            pairs["p_perm"] = pp
    if "q_perm" in pairs and (cfg.fdr_method == "permutation" or pairs["q_bh"].isna().all()):
        pairs["q_value"] = pairs["q_perm"]
        pairs["p_used"] = pairs["p_perm"]
        fdr_used = "permutation"
    else:
        pairs["q_value"] = pairs["q_bh"]
        pairs["p_used"] = pairs["p_value"]
        fdr_used = "bh"
    pairs = classify_pairs(pairs, cfg)
    feats = feature_statistics(prep)
    ptab = pathway_table(prep, pairs, perm, cfg)
    if perm:
        perm.pop("_set_p", None)
        perm.pop("_set_obs", None)
    rtab = reaction_table(prep, pairs)
    hyps = build_hypotheses(prep, pairs, feats, cfg, ptab)
    res = RewiringResult(prep=prep, cfg=cfg, pairs=pairs, features=feats, hypotheses=hyps, permutation=perm,
                         pathway_tests=ptab, reactions=rtab)
    res.fdr_used = fdr_used
    return res


# ---------------------------------------------------------------------------
# Simulated demo with planted ground truth
# ---------------------------------------------------------------------------
_DEMO_SPECS = [
    ("Glucose", "gly", -0.4), ("Glucose-6-phosphate", "gly", 0.2), ("Fructose-1,6-bisphosphate", "gly", 0.5),
    ("3-Phosphoglycerate", "gly", 0.4), ("Phosphoenolpyruvate", "gly", 0.3), ("Pyruvate", "gly", 0.6), ("Lactate", "gly", 1.2),
    ("Citrate", "tca", 0.0), ("cis-Aconitate", "tca", 0.1), ("Isocitrate", "tca", 0.0), ("alpha-Ketoglutarate", "tca", 0.3),
    ("Succinate", "tca", 0.7), ("Fumarate", "tca", 0.5), ("Malate", "tca", 0.2),
    ("Glutamine", "gln", -0.8), ("Glutamate", "gln", 0.5), ("Aspartate", "gln", 0.0), ("Asparagine", "gln", -0.2), ("Proline", "gln", 0.3),
    ("Tryptophan", "trp", -0.7), ("Kynurenine", "trp", 1.3), ("Kynurenic acid", "trp", 0.9), ("3-Hydroxykynurenine", "trp", 0.6),
    ("Quinolinate", "trp", 1.0), ("Serotonin", "trp", -0.3),
    ("Acetylcarnitine", "ac", -0.2), ("Propionylcarnitine", "ac", -0.1), ("Butyrylcarnitine", "ac", -0.3),
    ("Octanoylcarnitine", "ac", -0.5), ("Palmitoylcarnitine", "ac", -1.0), ("Oleoylcarnitine", "ac", -0.9),
    ("Serine", "oc", 0.4), ("Glycine", "oc", 0.2), ("Methionine", "oc", -0.2), ("S-Adenosylmethionine", "oc", 0.1),
    ("S-Adenosylhomocysteine", "oc", 0.3),
    ("Hypoxanthine", "pur", 0.2), ("Xanthine", "pur", 0.1), ("Inosine", "pur", -0.1),
    ("Glutathione (GSH)", "red", -0.3), ("Glutathione disulfide (GSSG)", "red", 0.6), ("Cysteine", "red", 0.1), ("Taurine", "red", 0.0),
]
DEMO_PATHWAYS = {
    "gly": "Glycolysis", "tca": "TCA Cycle", "gln": "Glutamine / Amino Acid Metabolism",
    "trp": "Tryptophan–Kynurenine Metabolism", "ac": "Fatty Acid Oxidation (Acylcarnitines)",
    "oc": "One-Carbon Metabolism", "pur": "Purine Metabolism", "red": "Glutathione / Redox",
}


def simulate_rewiring_demo(n_per_group: int = 30, seed: int = 20260927, rewire: bool = True,
                           outliers: bool = True, batch_effect: bool = True) -> dict:
    """
    Simulated targeted panel (43 metabolites) for a Control vs Tumor comparison with
    planted rewiring. Returns raw peak areas in MetaboAI Pro's input format plus the
    ground truth.

    Planted in Tumor (when rewire=True):
      lost    : Glutamine/Glutamate/Aspartate x 7 TCA intermediates (21 pairs; TCA coupling
                replaced by an independent glutamine factor), Serotonin x Trp/kynurenine branch (5)
      flipped : Tryptophan x 4 kynurenine-branch metabolites (IDO1-driven catabolism)
      gained  : Lactate x 4 kynurenine-branch metabolites and Lactate x Tryptophan (5)
    Decoy (outliers=True): Taurine and Hypoxanthine both spiked in the same 3 tumor samples.
    Batch: two balanced batches; batch 2 shifts acylcarnitines and redox metabolites by +0.4
    log2 in both groups (a mean shift, so it must not create rewiring).
    Level changes: each metabolite has a tumor log2 fold change (some deliberately 0,
    e.g. Citrate and Aspartate, to show rewiring that fold-change analysis misses).
    """
    rng = np.random.default_rng(seed)
    n = 2 * n_per_group
    names = [s[0] for s in _DEMO_SPECS]
    base = 14 + rng.random(len(_DEMO_SPECS)) * 8
    X = np.zeros((len(_DEMO_SPECS), n))
    group = np.array([0] * n_per_group + [1] * n_per_group)
    batch = np.array(["B1" if s % 2 == 0 else "B2" for s in range(n)])
    outlier_samples = {n_per_group + 4, n_per_group + 17, n_per_group + 25}
    factors = ["gly", "tca", "gln", "aa", "trp", "ido", "sero", "ac", "oc", "pur", "red", "tau"]
    for s in range(n):
        g = group[s] if rewire else 0
        F = {k: rng.standard_normal() for k in factors}
        for i, (name, pw, fc) in enumerate(_DEMO_SPECS):
            if pw == "gly":
                lat = 0.3 * F["gly"] + 0.9 * F["ido"] if (name == "Lactate" and g == 1) else 0.8 * F["gly"]
            elif pw == "tca":
                lat = 0.8 * F["tca"]
            elif pw == "gln":
                lat = (0.8 * F["tca"] if g == 0 else 0.8 * F["gln"]) if name in ("Glutamine", "Glutamate", "Aspartate") else 0.7 * F["aa"]
            elif pw == "trp":
                if name == "Tryptophan":
                    lat = 0.7 * F["trp"] if g == 0 else -0.8 * F["ido"]
                elif name == "Serotonin":
                    lat = 0.7 * F["trp"] if g == 0 else 0.7 * F["sero"]
                else:
                    lat = 0.7 * F["trp"] if g == 0 else 0.8 * F["ido"]
            elif pw == "ac":
                lat = 0.8 * F["ac"]
            elif pw == "oc":
                lat = 0.7 * F["oc"]
            elif pw == "pur":
                lat = 0.7 * F["pur"]
            else:
                lat = 0.6 * F["tau"] if name == "Taurine" else 0.7 * F["red"]
            v = base[i] + lat + 0.35 * rng.standard_normal() + (fc if group[s] == 1 else 0.0)
            if batch_effect and batch[s] == "B2" and pw in ("ac", "red"):
                v += 0.4
            if outliers and group[s] == 1 and s in outlier_samples and name in ("Taurine", "Hypoxanthine"):
                v += 5.5
            X[i, s] = 2.0 ** v
    samples = [f"CTRL_{s + 1:02d}" if s < n_per_group else f"TUM_{s - n_per_group + 1:02d}" for s in range(n)]
    peak = pd.DataFrame(X, index=pd.Index(names, name="Metabolite"), columns=samples)
    meta = pd.DataFrame({
        "Sample": samples, "Group": ["Control" if g == 0 else "Tumor" for g in group],
        "IsQC": False, "Batch": batch,
    })
    annotations = pd.DataFrame({
        "Metabolite": names, "Method": "Targeted Metabolomics",
        "Pathway": [DEMO_PATHWAYS[s[1]] for s in _DEMO_SPECS],
    })

    def P(a, b):
        return frozenset((a, b))
    kyn = ["Kynurenine", "Kynurenic acid", "3-Hydroxykynurenine", "Quinolinate"]
    tca = [s[0] for s in _DEMO_SPECS if s[1] == "tca"]
    gly_up = [s[0] for s in _DEMO_SPECS if s[1] == "gly" and s[0] != "Lactate"]
    truth = {"lost": set(), "gained": set(), "flipped": set(), "ambiguous": set(), "decoy": set()}
    if rewire:
        truth["lost"] |= {P(a, b) for a in ("Glutamine", "Glutamate", "Aspartate") for b in tca}
        truth["lost"] |= {P("Serotonin", b) for b in kyn + ["Tryptophan"]}
        truth["flipped"] |= {P("Tryptophan", b) for b in kyn}
        truth["gained"] |= {P("Lactate", b) for b in kyn + ["Tryptophan"]}
        truth["ambiguous"] |= {P("Lactate", b) for b in gly_up}   # weakened coupling (r ~0.84 -> ~0.3)
    if outliers:
        truth["decoy"].add(P("Taurine", "Hypoxanthine"))
        # every other pair touching a spiked metabolite is outlier-contaminated in Tumor
        truth["ambiguous"] |= {P(a, b) for a in ("Taurine", "Hypoxanthine") for b in names if b != a} - truth["decoy"]
    return {"peak_df": peak, "meta": meta, "annotations": annotations, "truth": truth}
