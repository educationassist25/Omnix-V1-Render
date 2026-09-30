"""
normalization.py - Reference-channel (TMT) normalization, Median-IQR (label-free)
normalization, and Log2 transformation.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use("Agg")


# ---------------------------------------------------------------------------
# 1. Reference-channel normalization (TMT)
# ---------------------------------------------------------------------------
def reference_channel_normalize(intensity_df: pd.DataFrame, reference_channel: str) -> pd.DataFrame:
    """
    Normalized Protein Intensity = Channel Protein Intensity / Reference Channel Protein Intensity.

    TMT (Tandem Mass Tag) plexes commonly include one channel that is a pooled/bridge
    reference sample (the same pooled material run in every plex, or a designated
    reference channel within the plex). Every protein's intensity in every other
    channel is expressed as a ratio to that same protein's intensity in the reference
    channel, which corrects for plex-to-plex and channel-to-channel loading differences
    the same way ISTD normalization does in targeted small-molecule assays -- just
    computed per-protein (row) against a reference SAMPLE (column) instead of against a
    single internal-standard feature.

    reference_channel must be a column (sample) present in intensity_df. The reference
    channel itself is dropped from the output (its own ratio is trivially 1 for every
    protein and carries no information).
    """
    if reference_channel not in intensity_df.columns:
        raise ValueError(f"Reference channel '{reference_channel}' not found among samples.")
    ref_col = intensity_df[reference_channel]
    normalized = intensity_df.drop(columns=reference_channel).div(ref_col.replace(0, np.nan), axis=0)
    return normalized


# Backward-compatible alias retained for any external callers.
istd_normalize = reference_channel_normalize


# ---------------------------------------------------------------------------
# 2. Median centering normalization
# ---------------------------------------------------------------------------
def median_center_normalize(df: pd.DataFrame, exclude_cols=None) -> pd.DataFrame:
    """
    Median centering normalization (per-sample):
        X_norm_ij = X_ij - Median_j(X) + Grand_Median

    For each sample (column), subtract that sample's own median (computed across
    proteins) so every sample is recentered to the same level — the standard fix for
    systematic sample-to-sample loading/injection offsets. Grand_Median (the median of
    all per-sample medians) is added back so the overall scale of the data is
    preserved rather than collapsed to zero.

    exclude_cols: sample column(s) to exclude entirely — both from the per-sample
    median/grand-median calculation AND from the output. Use this for a TMT
    reference/pooled channel, which isn't a directly-comparable study sample and
    would distort the grand median if included.
    """
    exclude_cols = list(exclude_cols) if exclude_cols else []
    calc_cols = [c for c in df.columns if c not in exclude_cols]
    sub = df[calc_cols]
    sample_medians = sub.median(axis=0)
    grand_median = sample_medians.median()
    return sub.sub(sample_medians, axis=1).add(grand_median)


# ---------------------------------------------------------------------------
# 3. Global MAD-based variance scaling
# ---------------------------------------------------------------------------
def mad_scale_normalize(df: pd.DataFrame, exclude_cols=None):
    """
    Global normalization by Median Absolute Deviation (MAD)-based variance scaling:
        X_norm = X / (MAD_global x 1.4826)

    A single global MAD is computed across every value in the dataset at once (not
    per-protein or per-sample), then used as one scaling factor applied uniformly —
    standardizing the overall variance/spread of the whole dataset in one step. The
    1.4826 constant rescales MAD to be comparable to a standard deviation for
    normally-distributed data (the standard MAD-to-SD consistency correction). Pure
    scaling (no centering), so — unlike Median Centering or Median-IQR — this never
    introduces negative values from originally-positive protein intensities.

    exclude_cols: sample column(s) to exclude entirely — both from the global MAD
    calculation AND from the output. Use this for a TMT reference/pooled channel,
    which isn't a directly-comparable study sample and would distort the global MAD
    if included.

    Returns (normalized_df, mad_scale_used).
    """
    exclude_cols = list(exclude_cols) if exclude_cols else []
    calc_cols = [c for c in df.columns if c not in exclude_cols]
    sub = df[calc_cols]
    vals = sub.values.astype(float).flatten()
    vals = vals[~np.isnan(vals)]
    global_median = np.median(vals) if vals.size else 0.0
    mad = np.median(np.abs(vals - global_median)) if vals.size else 0.0
    mad_scale = mad * 1.4826 if mad > 0 else 1.0
    return sub / mad_scale, mad_scale


# ---------------------------------------------------------------------------
# 4. IQR normalization
# ---------------------------------------------------------------------------
def iqr_normalize(
    df: pd.DataFrame,
    axis: str = "sample",
    batch_map: pd.Series = None,
    avoid_nan: bool = False,
    stats_source: pd.DataFrame = None,
) -> pd.DataFrame:
    """
    Median-IQR normalization (robust scaling):

        X_norm = (X - Median(X)) / IQR(X)

    Default (axis='sample'), indexed by protein i and sample j:

        normalized[i, j] = (X[i, j] - median_j) / IQR_j

    where median_j and IQR_j are computed from sample j's values across all
    proteins i -- correcting for total-intensity/loading differences between
    samples/channels. This is the recommended default for label-free
    proteomics, applied to the strict Log2-transformed intensities.

    Recommended usage: apply this AFTER log2 transformation, not before. Robust-scaled
    values are frequently negative (anything below the median), so taking log2 of the
    OUTPUT of this function will produce NaNs for roughly half the data. If you need
    both steps, always log2 first, then robust-scale the log2 values.

    Parameters
    ----------
    df : pd.DataFrame
        Numeric protein x sample dataframe.

    axis : str
        'sample' (default):
            Normalize each sample/channel across proteins:
            normalized[i, j] = (X[i, j] - median_j) / IQR_j.

        'feature':
            Normalize each protein (row) across samples:
            normalized[i, j] = (X[i, j] - median_i) / IQR_i.

        'batch':
            Normalize each protein independently within each batch.
            Requires batch_map.

    batch_map : pd.Series, optional
        Series indexed by sample name containing batch labels.

    avoid_nan : bool, default False
        When a feature/sample/batch-group has zero spread (IQR = 0, i.e.
        every value is identical), the default behavior (avoid_nan=False,
        unchanged from prior releases) divides by NaN, propagating NaN to
        every value in that group. When avoid_nan=True, a zero IQR is
        instead treated as 1 (no scaling), so a zero-spread group is left
        centered at 0 rather than converted to missing. This does not
        change any value for groups that already have nonzero spread; it
        only affects the degenerate zero-IQR edge case, and never
        introduces NaN that wasn't already present in the input.

    stats_source : pd.DataFrame, optional
        When given, the per-sample (axis='sample') or per-feature
        (axis='feature') median and IQR are computed from `stats_source`
        instead of from `df` itself, and then applied to `df`. `df` is
        only ever the dataframe that gets transformed and returned --
        `stats_source` never appears in the output. Useful when `df` is a
        reduced dataframe (e.g. after an optional CV-based feature filter)
        but the per-sample scaling statistics should still reflect the
        complete detected-feature background at that pipeline stage.
        Ignored when axis='batch'. Defaults to None, which reproduces the
        prior behavior of computing statistics from `df` itself.

    Returns
    -------
    pd.DataFrame
        IQR-normalized dataframe.
    """

    def _safe_iqr(q1, q3):
        iqr = q3 - q1
        return iqr.replace(0, 1) if avoid_nan else iqr.replace(0, np.nan)

    if axis == "feature":
        src = stats_source if stats_source is not None else df
        median = src.median(axis=1)
        q1 = src.quantile(0.25, axis=1)
        q3 = src.quantile(0.75, axis=1)
        iqr = _safe_iqr(q1, q3)
        return (
            df
            .sub(median.reindex(df.index), axis=0)
            .div(iqr.reindex(df.index), axis=0)
        )

    elif axis == "sample":
        src = stats_source if stats_source is not None else df
        median = src.median(axis=0)
        q1 = src.quantile(0.25, axis=0)
        q3 = src.quantile(0.75, axis=0)
        iqr = _safe_iqr(q1, q3)
        return (
            df
            .sub(median.reindex(df.columns), axis=1)
            .div(iqr.reindex(df.columns), axis=1)
        )

    elif axis == "batch":
        if batch_map is None:
            raise ValueError("batch_map (sample -> batch label) is required for batch-specific normalization.")
        out = df.copy()
        for batch in batch_map.dropna().unique():
            cols = batch_map.index[batch_map == batch].tolist()
            cols = [c for c in cols if c in df.columns]
            if not cols:
                continue
            sub = df[cols]
            median = sub.median(axis=1)
            q1 = sub.quantile(0.25, axis=1)
            q3 = sub.quantile(0.75, axis=1)
            iqr = _safe_iqr(q1, q3)
            out[cols] = sub.sub(median, axis=0).div(iqr, axis=0)
        return out

    else:
        raise ValueError("axis must be one of 'feature', 'sample', 'batch'")


# ---------------------------------------------------------------------------
# 5. Log2 transformation
# ---------------------------------------------------------------------------
def log2_transform(df: pd.DataFrame):
    """
    Strict log2(x) -- no pseudo-count, no constant added. Appropriate when
    missingness (zeros/NaNs) has already been fully resolved upstream, in the
    Data Cleaning & Imputation step, so every remaining value is strictly
    positive and a "handle zeros" mechanism inside the log transform itself is
    no longer needed.

    Raises ValueError -- rather than silently producing -inf or NaN -- if any
    zero, negative, or missing value is still present, since that means
    cleaning/imputation upstream hasn't fully resolved missingness yet (or an
    earlier normalization step produced negative values); either is worth
    surfacing explicitly rather than corrupting every statistic computed on
    this data downstream.

    Returns (transformed_df, 0.0) -- the second value is kept only so existing
    call sites that unpack a (transformed_df, constant) pair don't need to
    change; the constant is always 0 now since none is ever added.
    """
    n_nan = int(df.isna().sum().sum())
    n_nonpos = int((df <= 0).sum().sum())
    if n_nan or n_nonpos:
        problems = []
        if n_nonpos:
            problems.append(f"{n_nonpos} zero/negative value(s)")
        if n_nan:
            problems.append(f"{n_nan} missing (NaN) value(s)")
        raise ValueError(
            "Strict log2(x) requires every value to be strictly positive, but this data "
            f"still has {' and '.join(problems)}. Resolve these first — in Data Cleaning & "
            "Imputation for missing/zero biological values, or by checking the source data "
            "for QC replicates — before applying Log2."
        )
    return np.log2(df), 0.0


def shift_and_log2_transform(df: pd.DataFrame, shift: float = None):
    """
    LEGACY / retained for backward compatibility only -- not used by the current
    app workflow, and not affected by log2_transform()'s move to strict log2(x)
    (this function never relied on that behavior).

    Shift data to be strictly positive (if needed), then log2 transform.

    Use this — instead of log2_transform() — for data that has already passed through
    a CENTERING normalization such as the classic Median-IQR formula,
    (X - Median) / IQR, which produces negative values for any point below the median.
    log2_transform() requires strictly positive input and will raise rather than handle
    that; this function instead shifts the entire dataset by a constant just large
    enough to make the global minimum slightly positive, preserving every value's
    relative position, then logs the shifted data.

    If the data is already all-positive, shift defaults to 0 (equivalent to a plain
    log2 with no shift).

    Returns (transformed_df, shift_used).
    """
    finite_vals = df.values[np.isfinite(df.values)]
    min_val = finite_vals.min() if finite_vals.size else 0.0
    if shift is None:
        if min_val <= 0:
            shift = abs(min_val) + max(abs(min_val) * 0.01, 1e-3)
        else:
            shift = 0.0
    shifted = df + shift
    transformed = np.log2(shifted)
    return transformed, shift


def distribution_plots(before: pd.DataFrame, after: pd.DataFrame, sample_id: str = None):
    """
    Generate before/after density and box plots, each on ITS OWN appropriately-scaled
    panel. Raw protein intensities (before) typically span several orders of magnitude
    while log2-transformed values (after) span a much smaller range (~0-30); overlaying
    them on a single shared axis makes one distribution collapse to an invisible sliver.
    Using separate panels (with a log x-axis for the raw-scale density, since raw
    intensities are all positive and right-skewed) keeps both distributions legible.
    """
    b = before.values.flatten()
    b = b[~np.isnan(b) & (b > 0)]
    a = after.values.flatten()
    a = a[~np.isnan(a)]

    # Single row of 4 panels (Before density, Before box, After density, After
    # box) instead of a 2x2 block -- a compact 4x1 strip that reads at a glance
    # and takes up far less vertical space than the old 2x2 layout. Displayed
    # at full tab width (utils.render_figure uses width="stretch"), so the
    # figure is sized to fill that width well rather than looking small with
    # empty space on either side.
    fig, axes = plt.subplots(1, 4, figsize=(15, 3.4))
    fs = 9.5

    # Before: density (log-x, since raw intensities are positive and span orders of magnitude)
    if b.size:
        axes[0].hist(b, bins=60, color="#C44E52", alpha=0.8)
        axes[0].set_xscale("log")
    axes[0].set_title("Before: Density\n(raw, log scale)", fontsize=fs, fontweight="bold")
    axes[0].set_ylabel("Count", fontsize=fs)

    # Before: box plot (own y-axis, raw scale)
    if b.size:
        try:
            axes[1].boxplot([b], tick_labels=["Before"])
        except TypeError:
            axes[1].boxplot([b], labels=["Before"])
        axes[1].set_yscale("log")
    axes[1].set_title("Before: Box Plot\n(raw, log scale)", fontsize=fs, fontweight="bold")

    # After: density (linear x, already log2 scale)
    if a.size:
        axes[2].hist(a, bins=60, color="#55A868", alpha=0.8)
    axes[2].set_title("After: Density\n(log2-transformed)", fontsize=fs, fontweight="bold")
    axes[2].set_xlabel("log2(intensity)", fontsize=fs)
    axes[2].set_ylabel("Count", fontsize=fs)

    # After: box plot (own y-axis, log2 scale)
    if a.size:
        try:
            axes[3].boxplot([a], tick_labels=["After"])
        except TypeError:
            axes[3].boxplot([a], labels=["After"])
    axes[3].set_title("After: Box Plot\n(log2-transformed)", fontsize=fs, fontweight="bold")

    for ax in axes:
        ax.tick_params(labelsize=fs - 1)

    fig.tight_layout()
    return fig
