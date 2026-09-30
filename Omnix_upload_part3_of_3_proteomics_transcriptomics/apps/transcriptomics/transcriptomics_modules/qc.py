"""
qc.py - Quality Control validation module.

Implements:
  - Coefficient of Variation (CV) calculation across QC replicates
  - CV distribution / histogram plots
  - CV-based filtering table (<=20% acceptable, >20% variable)
  - QC sample correlation matrix (the only QC visualization retained per request —
    QC PCA, hierarchical clustering, and sample distance heatmap were removed)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use("Agg")


def calculate_cv(qc_df: pd.DataFrame) -> pd.DataFrame:
    """
    qc_df: features x QC-samples matrix (raw expression values).
    Returns a DataFrame with Mean, SD, CV(%) and Quality flag per feature.
    """
    mean = qc_df.mean(axis=1, skipna=True)
    sd = qc_df.std(axis=1, skipna=True, ddof=1)
    cv = (sd / mean.replace(0, np.nan)) * 100
    quality = np.where(cv <= 20, "Acceptable", "Variable")
    out = pd.DataFrame({"Mean": mean, "SD": sd, "CV(%)": cv, "Quality": quality}, index=qc_df.index)
    return out.sort_values("CV(%)")


def cv_distribution_plot(cv_table: pd.DataFrame):
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(cv_table["CV(%)"].dropna(), bins=30, color="#4C72B0", edgecolor="white")
    ax.axvline(20, color="red", linestyle="--", label="20% threshold")
    ax.set_xlabel("Coefficient of Variation (%)")
    ax.set_ylabel("Number of features")
    ax.set_title("QC CV Distribution")
    ax.legend()
    fig.tight_layout()
    return fig


def cv_histogram(cv_table: pd.DataFrame):
    fig, ax = plt.subplots(figsize=(6, 4))
    counts = cv_table["Quality"].value_counts()
    ax.bar(counts.index, counts.values, color=["#55A868", "#C44E52"])
    ax.set_ylabel("Number of features")
    ax.set_title("Feature Count by CV Quality")
    for i, v in enumerate(counts.values):
        ax.text(i, v, str(v), ha="center", va="bottom")
    fig.tight_layout()
    return fig


def sample_correlation_matrix(qc_log_data: pd.DataFrame):
    """QC-only Pearson correlation matrix. qc_log_data must contain ONLY QC replicates."""
    corr = qc_log_data.corr(method="pearson")
    fig, ax = plt.subplots(figsize=(6, 4))
    im = ax.imshow(corr, cmap="viridis", vmin=corr.values.min(), vmax=1, interpolation="nearest")
    ax.set_xticks(range(len(corr.columns)))
    ax.set_xticklabels(corr.columns, rotation=90, fontsize=6)
    ax.set_yticks(range(len(corr.columns)))
    ax.set_yticklabels(corr.columns, fontsize=6)
    ax.set_title("QC Sample Correlation Matrix (QC replicates only)")
    fig.colorbar(im, ax=ax, shrink=0.8, label="Pearson r")
    fig.tight_layout()
    return fig, corr


def expression_distribution_plot(log_df: pd.DataFrame, title: str = "Per-Sample Expression Distribution"):
    """Per-sample boxplot of log-scale expression values -- flags samples whose overall
    level or spread looks off relative to the rest of the cohort.

    Figsize is fixed at (6, 4) -- matching every other QC grid panel -- rather
    than growing with the sample count. render_figure_grid() renders every
    panel at the same DPI and stretches it to the same on-screen column
    width, so a wider-than-standard figure here would only get squeezed back
    down (shrinking its fonts further, not enlarging them). With the
    footprint now fixed, legibility for cohorts with many samples instead
    comes from a larger tick/label/title font plus thinning which sample
    labels are actually drawn (every data point is still boxplotted; only
    the x-tick LABELS are sparsified) -- standard practice for many-sample
    RNA-seq QC panels.
    """
    fig, ax = plt.subplots(figsize=(6, 4))
    order = list(log_df.columns)
    data = [log_df[c].dropna().values for c in order]
    ax.boxplot(data, tick_labels=order, showfliers=False, patch_artist=True,
               boxprops=dict(facecolor="#4C72B0", alpha=0.6))
    ax.set_xticks(range(1, len(order) + 1))
    n = len(order)
    # Show every sample's label up to ~20 samples; beyond that, thin the
    # labels (every Nth) so they stay readable at a larger font instead of
    # overlapping into an unreadable smear -- the boxplot itself still shows
    # every sample's data regardless of how many labels are drawn.
    step = max(1, -(-n // 20))  # ceil(n / 20)
    shown_ticks = list(range(1, n + 1, step))
    shown_labels = [order[i - 1] for i in shown_ticks]
    ax.set_xticks(shown_ticks)
    ax.set_xticklabels(shown_labels, rotation=90, fontsize=9)
    ax.set_ylabel("log2(CPM + 1) / VST value", fontsize=10)
    ax.set_title(title, fontsize=11)
    ax.tick_params(axis="y", labelsize=9)
    fig.tight_layout()
    return fig


def all_sample_correlation_matrix(log_df: pd.DataFrame,
                                   title: str = "Sample Correlation Matrix (all biological samples)"):
    """Pearson correlation across every biological sample (as opposed to
    sample_correlation_matrix above, which is restricted to QC replicates).

    Figsize is fixed at (6, 4) -- matching every other QC grid panel --
    rather than growing with the sample count. render_figure_grid() renders
    every panel at the same DPI and stretches it to the same on-screen
    column width, so a larger-than-standard figure here would only get
    squeezed back down on screen (shrinking its fonts further, not
    enlarging them). The full correlation matrix (every sample x every
    sample) is still computed and plotted -- only the tick LABELS are
    thinned for cohorts with many samples, so the heatmap stays legible at a
    larger font instead of overlapping into an unreadable smear.
    """
    corr = log_df.corr(method="pearson")
    n = len(corr)
    fig, ax = plt.subplots(figsize=(6, 4))
    im = ax.imshow(corr, cmap="viridis", vmin=float(corr.values.min()), vmax=1, interpolation="nearest")
    # Show every sample's label up to ~15 samples (fewer than the boxplot's
    # threshold since labels appear on both axes here); beyond that, thin to
    # every Nth label. Every row/column of the matrix is still drawn.
    step = max(1, -(-n // 15))  # ceil(n / 15)
    shown_ticks = list(range(0, n, step))
    shown_labels = [corr.columns[i] for i in shown_ticks]
    ax.set_xticks(shown_ticks)
    ax.set_xticklabels(shown_labels, rotation=90, fontsize=9)
    ax.set_yticks(shown_ticks)
    ax.set_yticklabels(shown_labels, fontsize=9)
    ax.set_title(title, fontsize=11)
    cbar = fig.colorbar(im, ax=ax, shrink=0.8, label="Pearson r")
    cbar.ax.tick_params(labelsize=9)
    cbar.set_label("Pearson r", fontsize=10)
    fig.tight_layout()
    return fig, corr


def detect_outlier_samples(log_df: pd.DataFrame, corr_z_thresh: float = 2.5,
                            pca_z_thresh: float = 3.0) -> pd.DataFrame:
    """
    Flags potential outlier samples by two independent criteria:

    1. Mean pairwise Pearson correlation to every other sample -- a sample whose
       mean correlation sits more than `corr_z_thresh` robust (MAD-scaled)
       deviations below the cohort median is flagged.
    2. PCA distance -- a sample farther than `pca_z_thresh` SD from the PC1/PC2
       centroid (computed on the same matrix passed in) is flagged.

    Returns a DataFrame indexed by sample: Mean_Correlation, Correlation_Z, PC1,
    PC2, PCA_Distance_Z, Flag_Reason, Is_Outlier.
    """
    corr = log_df.corr(method="pearson")
    n = corr.shape[1]
    mean_corr = (corr.sum(axis=1) - 1) / max(1, n - 1)
    med = mean_corr.median()
    mad = (mean_corr - med).abs().median() * 1.4826
    mad = mad if mad > 1e-9 else 1e-9
    corr_z = (mean_corr - med) / mad
    low_corr_flag = corr_z < -corr_z_thresh

    pc1 = pd.Series(np.nan, index=log_df.columns)
    pc2 = pd.Series(np.nan, index=log_df.columns)
    pca_z = pd.Series(0.0, index=log_df.columns)
    pca_flag = pd.Series(False, index=log_df.columns)
    X = log_df.T
    X = X.fillna(X.mean())
    if X.shape[0] >= 4 and X.shape[1] >= 2:
        from sklearn.decomposition import PCA
        pca = PCA(n_components=2)
        scores = pca.fit_transform(X.values)
        pc1 = pd.Series(scores[:, 0], index=X.index)
        pc2 = pd.Series(scores[:, 1], index=X.index)
        centroid = scores.mean(axis=0)
        dist = np.sqrt(((scores - centroid) ** 2).sum(axis=1))
        d_mean, d_sd = dist.mean(), (dist.std(ddof=1) or 1e-9)
        pca_z = pd.Series((dist - d_mean) / d_sd, index=X.index)
        pca_flag = pca_z > pca_z_thresh

    reasons = []
    for s in log_df.columns:
        r = []
        if bool(low_corr_flag.get(s, False)):
            r.append("Low correlation to other samples")
        if bool(pca_flag.get(s, False)):
            r.append("PCA outlier (distant from centroid)")
        reasons.append("; ".join(r))

    out = pd.DataFrame({
        "Mean_Correlation": mean_corr,
        "Correlation_Z": corr_z,
        "PC1": pc1,
        "PC2": pc2,
        "PCA_Distance_Z": pca_z,
        "Flag_Reason": reasons,
    }, index=log_df.columns)
    out["Is_Outlier"] = out["Flag_Reason"] != ""
    return out.sort_values(["Is_Outlier", "Correlation_Z"], ascending=[False, True])


def gene_expression_histogram(gene_df: pd.DataFrame, title: str = "Gene Expression Distribution"):
    """
    Histogram of all raw expression values in a features x samples matrix (every
    feature x every sample cell, flattened), on a log10 scale -- expression values
    span several orders of magnitude, so a linear-scale histogram would just show
    one tall bar near zero. Values <= 0 (already-imputed zeros/negatives shouldn't
    normally occur, but guard anyway) are dropped before logging.
    """
    values = gene_df.values.flatten()
    values = values[np.isfinite(values) & (values > 0)]
    log_values = np.log10(values)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(log_values, bins=50, color="#4C72B0", edgecolor="white")
    ax.set_xlabel("log10(Expression)")
    ax.set_ylabel("Count (feature × sample values)")
    ax.set_title(title)
    fig.tight_layout()
    return fig
