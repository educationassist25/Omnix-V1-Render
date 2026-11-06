"""
filtering_rnaseq.py - "Data Cleaning" module for raw-count RNA-seq data.

Bulk RNA-seq raw counts are essentially never missing (an unmapped gene is a
true biological/technical zero, not "missing data" in the intensity-based sense),
so "cleaning" here means the standard low-expression gene filter used before
any DESeq2/edgeR/limma analysis -- keeping genes that are detectably expressed
in enough samples to support statistical testing, and dropping the long tail
of genes with near-zero counts everywhere (which contribute noise and inflate
multiple-testing burden without any realistic chance of significance). No
imputation is performed for either input mode (see app.py's Data Cleaning
tab) -- any stray non-numeric cell in a raw-count matrix is treated as zero,
and a supplied 'Pre-normalized logCPM' matrix is used as-is.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use("Agg")


def counts_to_cpm(counts_df: pd.DataFrame) -> pd.DataFrame:
    """Counts-per-million, computed per sample (column) against that sample's own library size."""
    lib_size = counts_df.sum(axis=0).replace(0, np.nan)
    return counts_df.div(lib_size, axis=1) * 1e6


def filter_low_expression_genes(counts_df: pd.DataFrame, min_cpm: float = 1.0,
                                 min_fraction_samples: float = 0.2) -> tuple:
    """
    Keep genes with CPM >= min_cpm in at least min_fraction_samples of samples
    -- the standard "filterByExpr"-style rule used before edgeR/limma/DESeq2
    (a fixed count threshold alone is unfair across samples with very
    different library sizes, which CPM normalizes away before thresholding).

    Returns (filtered_counts_df, filter_report_df) where filter_report_df has
    one row per gene with columns Mean_CPM, N_Samples_Detected, Kept.
    """
    cpm = counts_to_cpm(counts_df)
    n_samples = counts_df.shape[1]
    min_n = max(1, int(np.ceil(min_fraction_samples * n_samples)))
    n_pass = (cpm >= min_cpm).sum(axis=1)
    keep_mask = n_pass >= min_n

    report = pd.DataFrame({
        "Mean_CPM": cpm.mean(axis=1),
        "N_Samples_Detected": n_pass,
        "Kept": keep_mask,
    }, index=counts_df.index).sort_values("Mean_CPM", ascending=False)

    filtered = counts_df.loc[keep_mask]
    return filtered, report


def filter_low_variance_genes(log_df: pd.DataFrame, min_variance_percentile: float = 0.0) -> tuple:
    """
    For logCPM-input mode: optionally drop the least-variable genes (bottom
    `min_variance_percentile`% by variance across samples) -- these carry the
    least information for downstream PCA/clustering/statistics. Percentile 0
    (default) keeps everything.
    """
    var = log_df.var(axis=1)
    if min_variance_percentile <= 0:
        keep_mask = pd.Series(True, index=log_df.index)
    else:
        cutoff = np.percentile(var.dropna(), min_variance_percentile)
        keep_mask = var >= cutoff
    report = pd.DataFrame({"Variance": var, "Kept": keep_mask}).sort_values("Variance", ascending=False)
    return log_df.loc[keep_mask.reindex(log_df.index, fill_value=False)], report


def library_size_plot(lib_summary: pd.DataFrame):
    """Combined 2-subplot Library Size + Genes Detected figure.

    Kept for backward compatibility; the QC Validation tab in app.py now
    renders these as two independent single-axis figures (see
    `library_size_only_plot` and `genes_detected_plot` below) so each gets
    its own title/Format/DPI/Download controls and shares the same figsize
    as every other QC panel.
    """
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    order = lib_summary["Library_Size"].sort_values(ascending=False).index
    axes[0].bar(range(len(order)), lib_summary.loc[order, "Library_Size"] / 1e6, color="#4C72B0")
    axes[0].set_xticks(range(len(order)))
    axes[0].set_xticklabels(order, rotation=90, fontsize=6)
    axes[0].set_ylabel("Library size (million reads)")
    axes[0].set_title("Library Size per Sample")

    axes[1].bar(range(len(order)), lib_summary.loc[order, "Genes_Detected"], color="#55A868")
    axes[1].set_xticks(range(len(order)))
    axes[1].set_xticklabels(order, rotation=90, fontsize=6)
    axes[1].set_ylabel("Genes detected (count > 0)")
    axes[1].set_title("Genes Detected per Sample")
    fig.tight_layout()
    return fig


def library_size_only_plot(lib_summary: pd.DataFrame):
    """Library Size per Sample as its own single-axis figure (same figsize
    convention -- (6, 4) -- as the other QC Validation panels)."""
    fig, ax = plt.subplots(figsize=(6, 4))
    order = lib_summary["Library_Size"].sort_values(ascending=False).index
    ax.bar(range(len(order)), lib_summary.loc[order, "Library_Size"] / 1e6, color="#4C72B0")
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=90, fontsize=6)
    ax.set_ylabel("Library size (million reads)")
    ax.set_title("Library Size per Sample")
    fig.tight_layout()
    return fig


def genes_detected_plot(lib_summary: pd.DataFrame):
    """Genes Detected per Sample as its own single-axis figure (same figsize
    convention -- (6, 4) -- as the other QC Validation panels)."""
    fig, ax = plt.subplots(figsize=(6, 4))
    order = lib_summary["Library_Size"].sort_values(ascending=False).index
    ax.bar(range(len(order)), lib_summary.loc[order, "Genes_Detected"], color="#55A868")
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=90, fontsize=6)
    ax.set_ylabel("Genes detected (count > 0)")
    ax.set_title("Genes Detected per Sample")
    fig.tight_layout()
    return fig


def filtering_summary_plot(report: pd.DataFrame, min_cpm: float):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    counts = report["Kept"].value_counts()
    labels = ["Kept" if k else "Filtered out" for k in counts.index]
    axes[0].bar(labels, counts.values, color=["#55A868", "#C44E52"])
    for i, v in enumerate(counts.values):
        axes[0].text(i, v, str(v), ha="center", va="bottom")
    axes[0].set_ylabel("Number of genes")
    axes[0].set_title("Gene Filtering Summary")

    vals = report["Mean_CPM"].replace(0, np.nan).dropna()
    vals = vals[vals > 0]
    axes[1].hist(np.log10(vals), bins=60, color="#4C72B0")
    axes[1].axvline(np.log10(min_cpm), color="red", linestyle="--", label=f"{min_cpm} CPM threshold")
    axes[1].set_xlabel("log10(Mean CPM)")
    axes[1].set_ylabel("Number of genes")
    axes[1].set_title("Mean Expression Distribution")
    axes[1].legend()
    fig.tight_layout()
    return fig
