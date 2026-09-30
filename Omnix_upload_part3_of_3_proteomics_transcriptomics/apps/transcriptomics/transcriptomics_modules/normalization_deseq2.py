"""
normalization_deseq2.py - Normalization for raw-count RNA-seq data via pyDESeq2's
median-of-ratios size factors, with a variance-stabilizing transform (VST) for
PCA/heatmap/boxplot visualization, and a CPM/log2-CPM path for the
'Pre-normalized logCPM' input mode (which skips this module's DESeq2 step
entirely -- there are no raw counts to model).

pyDESeq2 note: DeseqDataSet expects counts as SAMPLES (rows) x GENES (columns)
-- the transpose of the genes x samples convention used everywhere else in
this app. All transposing happens inside this module; every function here
still takes/returns genes x samples, matching the rest of the app.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use("Agg")

try:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference
    HAS_PYDESEQ2 = True
except Exception:  # pragma: no cover - exercised only when pydeseq2 isn't installed
    HAS_PYDESEQ2 = False


class NormalizationError(Exception):
    pass


def counts_to_log2cpm(counts_df: pd.DataFrame, pseudocount: float = 1.0) -> pd.DataFrame:
    """Simple log2(CPM + pseudocount) transform — used as a fast/robust fallback
    and for the 'Pre-normalized logCPM' input path's own optional re-derivation."""
    lib_size = counts_df.sum(axis=0).replace(0, np.nan)
    cpm = counts_df.div(lib_size, axis=1) * 1e6
    return np.log2(cpm + pseudocount)


def fit_deseq2(counts_df: pd.DataFrame, meta: pd.DataFrame, design_factors: list, n_cpus: int = 1,
               reference_level: str = None):
    """
    Fit a DeseqDataSet (size factors + dispersions + a full-model GLM) ONCE for
    a given design. Reused for every pairwise contrast drawn from the same
    design (fitting is the expensive step; each contrast's Wald test is cheap).

    counts_df    : genes x samples, raw integer counts (already gene-filtered).
    meta         : samples x covariates metadata, indexed by sample name,
                   containing every column named in design_factors.
    design_factors: list of metadata column names, e.g. ["Batch", "Diagnosis"]
                   -- pyDESeq2/DESeq2 convention: list the primary variable of
                   interest LAST, with any covariates to control for first.
    reference_level: if given, sets this level of the LAST design factor (the
                   primary variable of interest) as the model's reference/
                   baseline category, so its design-matrix coefficient names
                   read "<factor>[T.<other_level>]" for every other level --
                   required for apeGLM LFC shrinkage (DeseqStats.lfc_shrink),
                   which shrinks a named design coefficient, not an arbitrary
                   contrast. Leave None to use patsy's default (alphabetical)
                   reference.

    Returns a fitted DeseqDataSet.
    """
    if not HAS_PYDESEQ2:
        raise NormalizationError(
            "pydeseq2 is not installed. Run `pip install pydeseq2` (see requirements.txt) "
            "to enable DESeq2-based normalization and statistics for raw-count data."
        )
    counts_samples_x_genes = counts_df.T.copy()
    counts_samples_x_genes.index.name = None
    meta_aligned = meta.loc[counts_samples_x_genes.index, design_factors].copy()
    for col in design_factors:
        meta_aligned[col] = meta_aligned[col].astype(str)

    if reference_level is not None and design_factors:
        primary_col = design_factors[-1]
        levels = meta_aligned[primary_col].unique().tolist()
        if reference_level in levels:
            ordered = [reference_level] + [lvl for lvl in levels if lvl != reference_level]
            meta_aligned[primary_col] = pd.Categorical(meta_aligned[primary_col], categories=ordered)

    inference = DefaultInference(n_cpus=n_cpus)
    try:
        dds = DeseqDataSet(
            counts=counts_samples_x_genes.astype(int),
            metadata=meta_aligned,
            design="~ " + " + ".join(design_factors),
            refit_cooks=True,
            inference=inference,
        )
    except TypeError:
        # Older pyDESeq2 releases (<0.4) use design_factors=... instead of a
        # design="~ ..." formula string.
        dds = DeseqDataSet(
            counts=counts_samples_x_genes.astype(int),
            metadata=meta_aligned,
            design_factors=design_factors,
            refit_cooks=True,
            inference=inference,
        )
    dds.deseq2()
    return dds


def size_factors_table(dds) -> pd.DataFrame:
    """DESeq2 median-of-ratios size factor per sample -- values near 1.0 mean a
    sample's sequencing depth/composition is close to the dataset's geometric-mean
    reference; far from 1.0 flags samples that needed substantial rescaling.

    pyDESeq2 has stored this in different places across releases (`.obs["size_factors"]`
    in most current versions, `.obsm["size_factors"]` in some others) -- this checks
    both rather than assuming one.
    """
    if "size_factors" in dds.obs.columns:
        sf = pd.Series(dds.obs["size_factors"].values, index=dds.obs_names, name="Size_Factor")
    elif "size_factors" in dds.obsm:
        sf = pd.Series(np.asarray(dds.obsm["size_factors"]).ravel(), index=dds.obs_names, name="Size_Factor")
    else:
        raise NormalizationError(
            "Could not locate size factors on the fitted DeseqDataSet (checked "
            "obs['size_factors'] and obsm['size_factors']) — this pyDESeq2 version "
            "may use a different attribute name."
        )
    return sf.to_frame()


def normalized_counts(dds) -> pd.DataFrame:
    """DESeq2-normalized counts (raw counts / size factor), genes x samples."""
    if "normed_counts" in dds.layers:
        normed = pd.DataFrame(dds.layers["normed_counts"], index=dds.obs_names, columns=dds.var_names)
    else:
        # Fallback: compute directly from raw counts / size factors if this
        # pyDESeq2 version doesn't populate the 'normed_counts' layer.
        sf = size_factors_table(dds)["Size_Factor"].values
        raw = np.asarray(dds.X.todense()) if hasattr(dds.X, "todense") else np.asarray(dds.X)
        normed = pd.DataFrame(raw / sf[:, None], index=dds.obs_names, columns=dds.var_names)
    return normed.T


def vst_transform(dds) -> pd.DataFrame:
    """
    Variance-stabilizing transformation — DESeq2's recommended matrix for PCA,
    clustering, and heatmaps (log2-normalized counts alone still show
    higher variance for low-count genes; VST corrects that). Falls back to
    log2(normalized_counts + 1) if the installed pyDESeq2 version doesn't
    expose vst_fit/vst_transform (added in pydeseq2 >=0.4).
    Returns genes x samples.
    """
    try:
        dds.vst_fit()
        vst = dds.vst_transform()
        out = pd.DataFrame(vst, index=dds.obs_names, columns=dds.var_names)
        return out.T
    except (AttributeError, TypeError, KeyError, ValueError):
        normed = normalized_counts(dds)
        return np.log2(normed + 1)


def size_factor_plot(sf_table: pd.DataFrame):
    fig, ax = plt.subplots(figsize=(7, 4))
    order = sf_table["Size_Factor"].sort_values().index
    ax.bar(range(len(order)), sf_table.loc[order, "Size_Factor"], color="#4C72B0")
    ax.axhline(1.0, color="red", linestyle="--", linewidth=1, label="Size factor = 1.0")
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=90, fontsize=6)
    ax.set_ylabel("DESeq2 size factor")
    ax.set_title("Median-of-Ratios Size Factors per Sample")
    ax.legend()
    fig.tight_layout()
    return fig


def raw_counts_distribution_plot(raw_counts: pd.DataFrame):
    """'Before' half of the before/after distribution comparison: a single-panel
    histogram of raw count values (log10 scale) prior to DESeq2 normalization."""
    b = raw_counts.values.flatten()
    b = b[np.isfinite(b) & (b > 0)]
    fig, ax = plt.subplots(figsize=(6, 4))
    if b.size:
        ax.hist(np.log10(b), bins=60, color="#C44E52", alpha=0.85)
    ax.set_title("Before: Raw Counts (log10 scale)")
    ax.set_xlabel("log10(count + 1)")
    ax.set_ylabel("Frequency")
    fig.tight_layout()
    return fig


def vst_distribution_plot(vst_df: pd.DataFrame):
    """'After' half of the before/after distribution comparison: a single-panel
    histogram of DESeq2 variance-stabilized (VST) values."""
    a = vst_df.values.flatten()
    a = a[np.isfinite(a)]
    fig, ax = plt.subplots(figsize=(6, 4))
    if a.size:
        ax.hist(a, bins=60, color="#55A868", alpha=0.85)
    ax.set_title("After: DESeq2-normalized (VST)")
    ax.set_xlabel("VST value")
    ax.set_ylabel("Frequency")
    fig.tight_layout()
    return fig


def before_after_distribution_plot(raw_counts: pd.DataFrame, vst_df: pd.DataFrame):
    """Combined two-panel before/after figure, kept for any external callers that
    still want a single side-by-side comparison figure rather than the two
    separate panels used by the Normalization tab's figure grid."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    b = raw_counts.values.flatten()
    b = b[np.isfinite(b) & (b > 0)]
    if b.size:
        axes[0].hist(np.log10(b), bins=60, color="#C44E52", alpha=0.85)
    axes[0].set_title("Before: Raw Counts (log10 scale)")
    axes[0].set_xlabel("log10(count + 1)")
    axes[0].set_ylabel("Frequency")

    a = vst_df.values.flatten()
    a = a[np.isfinite(a)]
    if a.size:
        axes[1].hist(a, bins=60, color="#55A868", alpha=0.85)
    axes[1].set_title("After: DESeq2-normalized (VST)")
    axes[1].set_xlabel("VST value")
    fig.tight_layout()
    return fig
