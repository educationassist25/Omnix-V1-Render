"""
normalization.py - ISTD normalization, IQR normalization, and strict Log2 transformation.

Log2 transformation is strictly:

    log2(x)

No pseudocount, constant addition, zero replacement, shifting, or arbitrary offset
is applied. Values <= 0 cannot be log2 transformed and are converted to NaN.
"""

import re

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# 1. ISTD normalization
# ---------------------------------------------------------------------------
def istd_normalize(peak_df: pd.DataFrame, istd_name: str) -> pd.DataFrame:
    """
    Internal Standard (ISTD) normalization.

    Normalized Peak Area =
        Endogenous Metabolite Peak Area / Internal Standard Peak Area

    The ISTD feature itself is removed from the returned dataframe.

    If an ISTD value is zero, division is undefined and the corresponding
    normalized values become NaN.
    """
    if istd_name not in peak_df.index:
        raise ValueError(
            f"Internal standard '{istd_name}' not found among features."
        )

    istd_row = peak_df.loc[istd_name]

    # Avoid division by zero.
    denominator = istd_row.replace(0, np.nan)

    normalized = (
        peak_df
        .drop(index=istd_name)
        .div(denominator, axis=1)
    )

    return normalized


# ---------------------------------------------------------------------------
# ISTD <-> metabolite mapping (targeted metabolomics / lipidomics)
# ---------------------------------------------------------------------------
# Row-annotation column that names each analyte's internal standard. The first
# column whose name matches one of these (case/space/underscore-insensitive) is used.
ISTD_COLUMN_NAMES = ("istd", "internalstandard", "is", "istdname", "matchedistd", "assignedistd")
POOLED_ISTD = "Pooled ISTD (geometric mean)"
MATCH_COLUMN = "ISTD column (row annotations)"
MATCH_CLASS = "Class match"
MATCH_SINGLE = "Only ISTD selected"
MATCH_POOLED = "Pooled fallback (no assigned ISTD)"

_ISTD_NAME_RE = re.compile(r"istd|internal[\s_-]*standard|^is[_\s-]", re.IGNORECASE)


def _norm_key(x) -> str:
    return re.sub(r"[\s_\-]+", "", str(x)).lower()


def find_istd_column(row_annotations):
    """Name of the row-annotation column that holds each feature's internal standard, or None."""
    if row_annotations is None or not hasattr(row_annotations, "columns"):
        return None
    for col in row_annotations.columns:
        if _norm_key(col) in ISTD_COLUMN_NAMES:
            return col
    return None


def _blank(v) -> bool:
    return v is None or (isinstance(v, float) and np.isnan(v)) or str(v).strip() in ("", "nan", "None", "NA")


def annotation_values(features, row_annotations, column, label=None, method=None) -> dict:
    """
    {feature: value of `column`} from the row annotations, for every feature that has one.

    Looks a feature up as "<label>::<feature>" first (multi-dataset annotation files may use
    the combined, prefixed names), then as the plain name. When a plain name appears more
    than once (e.g. the same compound measured by two methods), the row whose 'Method'
    equals `method` is preferred.
    """
    out = {}
    if row_annotations is None or column is None or column not in getattr(row_annotations, "columns", []):
        return out
    ann = row_annotations
    idx = pd.Index([str(i) for i in ann.index])
    has_method = method is not None and "Method" in ann.columns
    for feat in features:
        keys = ([f"{label}::{feat}"] if label else []) + [str(feat)]
        for key in keys:
            pos = np.flatnonzero(idx == key)
            if len(pos) == 0:
                continue
            if len(pos) > 1 and has_method:
                same = [p for p in pos if str(ann["Method"].iloc[p]) == str(method)]
                if same:
                    pos = same
            val = ann[column].iloc[pos[0]]
            if not _blank(val):
                out[feat] = str(val).strip()
                break
    return out


def detect_istds(features, row_annotations=None, label=None, method=None) -> list:
    """
    Internal-standard features among `features`, in their original order:
      * any feature named in the row annotations' ISTD column (an ISTD row may name itself), or
      * any feature whose name marks it as a standard ("ISTD", "Internal Standard", "IS_...").
    """
    features = list(features)
    col = find_istd_column(row_annotations)
    referenced = set()
    if col is not None:
        by_key = {_norm_key(f): f for f in features}
        for v in row_annotations[col].dropna().astype(str):
            f = by_key.get(_norm_key(v))
            if f is not None:
                referenced.add(f)
    return [f for f in features if f in referenced or _ISTD_NAME_RE.search(str(f))]


def default_istd(features, row_annotations=None, label=None, method=None):
    """
    The dataset's internal standard: the ISTD named most often in the row annotations' ISTD
    column for these features, otherwise the first detected ISTD row; None if there is none.
    """
    features = list(features)
    istds = detect_istds(features, row_annotations, label, method)
    if not istds:
        return None
    col = find_istd_column(row_annotations)
    if col is not None:
        by_key = {_norm_key(i): i for i in istds}
        named = [by_key.get(_norm_key(v)) for f, v in
                 annotation_values(features, row_annotations, col, label, method).items() if f not in istds]
        named = [n for n in named if n is not None]
        if named:
            return max(dict.fromkeys(named), key=named.count)
    return istds[0]


def build_istd_mapping(features, istd_names, row_annotations=None, label=None, method=None) -> pd.DataFrame:
    """
    Decide which internal standard normalizes each (non-ISTD) feature.

    Priority for every feature:
      1. the ISTD named for it in the row annotations' ISTD column, if that ISTD is selected;
      2. otherwise a selected ISTD whose row-annotation 'Class' equals the feature's Class;
      3. otherwise, when exactly one ISTD is selected, that ISTD;
      4. otherwise the pooled ISTD: the per-sample geometric mean of all selected ISTDs.

    Returns a DataFrame indexed by feature with columns 'ISTD', 'Matched by', 'Class',
    'Annotated ISTD' (what the annotation asked for, even if it was not selected).
    """
    istd_names = list(istd_names)
    if not istd_names:
        raise ValueError("At least one internal standard feature must be selected.")
    istd_set = set(istd_names)
    by_key = {_norm_key(i): i for i in istd_names}
    features = [f for f in features if f not in istd_set]
    col = find_istd_column(row_annotations)
    assigned = annotation_values(features, row_annotations, col, label, method) if col else {}
    feat_cls = annotation_values(features, row_annotations, "Class", label, method)
    istd_cls = annotation_values(istd_names, row_annotations, "Class", label, method)
    cls_to_istd = {}
    for istd in istd_names:                                # first selected ISTD of a class wins
        c = istd_cls.get(istd)
        if c is not None:
            cls_to_istd.setdefault(_norm_key(c), istd)

    rows = []
    for feat in features:
        wanted = assigned.get(feat)
        cls = feat_cls.get(feat)
        if wanted is not None and _norm_key(wanted) in by_key:
            istd, how = by_key[_norm_key(wanted)], MATCH_COLUMN
        elif cls is not None and _norm_key(cls) in cls_to_istd:
            istd, how = cls_to_istd[_norm_key(cls)], MATCH_CLASS
        elif len(istd_names) == 1:
            istd, how = istd_names[0], MATCH_SINGLE
        else:
            istd, how = POOLED_ISTD, MATCH_POOLED
        rows.append({"Feature": feat, "ISTD": istd, "Matched by": how,
                     "Class": cls if cls is not None else "", "Annotated ISTD": wanted or ""})
    return pd.DataFrame(rows, columns=["Feature", "ISTD", "Matched by", "Class", "Annotated ISTD"]
                        ).set_index("Feature")


def istd_normalize_mapped(peak_df: pd.DataFrame, istd_names, row_annotations: pd.DataFrame = None,
                          label=None, method=None, drop_features=()):
    """
    ISTD normalization with an explicit ISTD <-> metabolite mapping (see build_istd_mapping):

        normalized[i, j] = peak_area[i, j] / peak_area[ISTD(i), j]

    Every selected ISTD row, plus any other feature listed in `drop_features` (e.g. internal
    standards that were detected but not selected), is removed from the result. A zero or
    missing ISTD value gives NaN for that sample (division is undefined).

    Returns (normalized_df, mapping_df).
    """
    istd_names = list(istd_names)
    if not istd_names:
        raise ValueError("At least one internal standard feature must be selected.")
    missing = [n for n in istd_names if n not in peak_df.index]
    if missing:
        raise ValueError(f"Internal standard(s) not found among features: {', '.join(map(str, missing))}")
    istd_rows = peak_df.loc[istd_names].apply(pd.to_numeric, errors="coerce").replace(0, np.nan)
    drop = [f for f in list(istd_names) + list(drop_features) if f in peak_df.index]
    remaining = peak_df.drop(index=list(dict.fromkeys(drop)))
    mapping = build_istd_mapping(remaining.index, istd_names, row_annotations, label, method)

    denominators = {n: istd_rows.loc[n] for n in istd_names}
    if (mapping["ISTD"] == POOLED_ISTD).any():
        with np.errstate(divide="ignore", invalid="ignore"):
            denominators[POOLED_ISTD] = np.exp(np.log(istd_rows).mean(axis=0, skipna=False))
    denom = pd.DataFrame([denominators[mapping.at[f, "ISTD"]].values for f in remaining.index],
                         index=remaining.index, columns=remaining.columns)
    normalized = remaining.apply(pd.to_numeric, errors="coerce").div(denom)
    return normalized, mapping


def istd_normalize_multi(peak_df: pd.DataFrame, istd_names, row_annotations: pd.DataFrame = None):
    """
    Backwards-compatible wrapper around istd_normalize_mapped().
    Returns (normalized_df, unmatched_features): the features that fell back to the pooled ISTD.
    """
    normalized, mapping = istd_normalize_mapped(peak_df, istd_names, row_annotations)
    return normalized, mapping.index[mapping["ISTD"] == POOLED_ISTD].tolist()


def istd_summary(mapping: pd.DataFrame, istd_names, cv_table: pd.DataFrame = None,
                 sample_df: pd.DataFrame = None) -> pd.DataFrame:
    """
    One row per selected ISTD: how many features it normalizes, its QC CV (%) when a QC CV
    table is available, and how many study samples have a zero/missing value for it.
    """
    counts = mapping["ISTD"].value_counts() if mapping is not None and len(mapping) else pd.Series(dtype=int)
    rows = []
    for istd in list(istd_names) + ([POOLED_ISTD] if POOLED_ISTD in counts.index else []):
        row = {"ISTD": istd, "Features normalized": int(counts.get(istd, 0))}
        if cv_table is not None and "CV(%)" in cv_table.columns:
            row["QC CV (%)"] = (round(float(cv_table.at[istd, "CV(%)"]), 1)
                                if istd in cv_table.index and pd.notna(cv_table.at[istd, "CV(%)"]) else np.nan)
        if sample_df is not None:
            if istd in sample_df.index:
                v = pd.to_numeric(sample_df.loc[istd], errors="coerce")
                row["Zero/missing samples"] = int((v.isna() | (v <= 0)).sum())
            else:
                row["Zero/missing samples"] = np.nan
        rows.append(row)
    return pd.DataFrame(rows).set_index("ISTD")


def keep_istds(filtered_df: pd.DataFrame, source_df: pd.DataFrame, istds) -> pd.DataFrame:
    """
    Put internal-standard rows that a feature filter (missingness, QC CV) removed back into
    `filtered_df`, in their original order: internal standards are needed for normalization
    and are removed only at the ISTD normalization step itself.
    """
    lost = [i for i in istds if i in source_df.index and i not in filtered_df.index]
    if not lost:
        return filtered_df
    keep = set(filtered_df.index) | set(lost)
    return source_df.loc[[f for f in source_df.index if f in keep]]


def drop_istds(df: pd.DataFrame, istds) -> pd.DataFrame:
    """Remove internal-standard rows (they are not biological features)."""
    istds = [i for i in istds if i in df.index]
    return df.drop(index=istds) if istds else df


# ---------------------------------------------------------------------------
# 2. IQR normalization
# ---------------------------------------------------------------------------
def iqr_normalize(
    df: pd.DataFrame,
    axis: str = "sample",
    batch_map: pd.Series = None,
    avoid_nan: bool = False,
    stats_source: pd.DataFrame = None
) -> pd.DataFrame:
    """
    Median-IQR normalization (robust scaling):

        X_norm = (X - Median(X)) / IQR(X)

    Default (axis='sample'), indexed by metabolite i and sample j:

        normalized[i, j] = (X[i, j] - median_j) / IQR_j

    where median_j and IQR_j are computed from sample j's values across all
    metabolites i -- correcting for total-intensity/loading differences
    between samples. This is the recommended default for the app's
    Untargeted Metabolomics/Lipidomics pipeline, applied to the strict
    Log2-transformed data (see stats_analysis.strict_log2 / strict_log2 is
    applied before this function, not by it): normalized[i, j] =
    (log2(raw[i, j]) - median_j) / IQR_j.

    Parameters
    ----------
    df : pd.DataFrame
        Numeric feature x sample dataframe.

    axis : str
        'sample' (default):
            Normalize each sample across features:
            normalized[i, j] = (X[i, j] - median_j) / IQR_j.

        'feature':
            Normalize each feature across samples:
            normalized[i, j] = (X[i, j] - median_i) / IQR_i.

        'batch':
            Normalize each feature independently within each batch.
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
        `stats_source` never appears in the output.

        This is for the case where `df` is a reduced dataframe (e.g. after
        an optional CV-based feature filter removed some rows), but the
        per-sample scaling statistics should still reflect the complete
        detected-feature background at that pipeline stage (every feature
        present before the optional filter), matching the recommended
        practice of computing median_j/IQR_j from the full feature table
        rather than a smaller, curated/exported subset. Concretely:

            normalized[i, j] = (df[i, j] - median_j(stats_source)) / IQR_j(stats_source)

        For axis='sample', `stats_source` must share the same columns
        (samples) as `df`; its rows (features) may differ or be a
        superset. For axis='feature', it must share the same index
        (features); its columns (samples) may differ or be a superset.
        Ignored when axis='batch'. Defaults to None, which reproduces the
        prior behavior of computing statistics from `df` itself.

    Returns
    -------
    pd.DataFrame
        IQR-normalized dataframe.

    Notes
    -----
    This is a centering/scaling transformation and can produce negative
    values. Therefore, the result should NOT be passed directly to the
    strict log2 transformation unless all values are positive.
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

        # Statistics are computed from `src` (features x its own samples) but
        # applied to `df`, aligned by feature (row) index -- df's rows are the
        # ones actually transformed/returned.
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

        # Statistics are computed from `src` (all features detected for each
        # sample column, e.g. before an optional CV-based feature filter) but
        # applied to `df` (the, possibly filtered, feature set actually
        # normalized/returned), aligned by sample (column).
        return (
            df
            .sub(median.reindex(df.columns), axis=1)
            .div(iqr.reindex(df.columns), axis=1)
        )

    elif axis == "batch":

        if batch_map is None:
            raise ValueError(
                "batch_map (sample -> batch label) is required "
                "for batch-specific normalization."
            )

        out = df.copy()

        for batch in batch_map.dropna().unique():

            cols = batch_map.index[
                batch_map == batch
            ].tolist()

            # Keep only columns that actually exist in the dataframe.
            cols = [c for c in cols if c in df.columns]

            if not cols:
                continue

            sub = df[cols]

            median = sub.median(axis=1)
            q1 = sub.quantile(0.25, axis=1)
            q3 = sub.quantile(0.75, axis=1)

            iqr = _safe_iqr(q1, q3)

            out[cols] = (
                sub
                .sub(median, axis=0)
                .div(iqr, axis=0)
            )

        return out

    else:
        raise ValueError(
            "axis must be one of 'feature', 'sample', or 'batch'"
        )


# ---------------------------------------------------------------------------
# 3. STRICT Log2 transformation
# ---------------------------------------------------------------------------
def log2_transform(df: pd.DataFrame) -> pd.DataFrame:
    """
    Strict Log2 transformation.

        transformed = log2(x)

    IMPORTANT
    ---------
    No pseudocount is added.

    No constant is added.

    No zero replacement is performed.

    No shifting is performed.

    No arbitrary offset is applied.

    Values <= 0 cannot be log2 transformed and are converted to NaN.

    NaN values remain NaN.

    Parameters
    ----------
    df : pd.DataFrame
        Numeric dataframe.

    Returns
    -------
    pd.DataFrame
        Strict log2-transformed dataframe.
    """

    work = df.copy()

    # Convert values <= 0 to NaN.
    work = work.where(work > 0, np.nan)

    # Strict mathematical log2(x).
    transformed = np.log2(work)

    return transformed


# ---------------------------------------------------------------------------
# 4. Distribution plots
# ---------------------------------------------------------------------------
def distribution_plots(
    before: pd.DataFrame,
    after: pd.DataFrame,
    sample_id: str = None
):
    """
    Generate before/after distribution and box plots.

    Before:
        Raw peak-area values are displayed using a logarithmic x/y scale
        where appropriate.

    After:
        Strict log2-transformed values are displayed on a linear scale.

    Values that are <= 0 before transformation are excluded from the
    transformed distribution because they become NaN under strict log2.
    """

    # ------------------------------------------------------------------
    # Before data
    # ------------------------------------------------------------------
    b = before.values.flatten()

    b = b[
        np.isfinite(b) &
        (b > 0)
    ]

    # ------------------------------------------------------------------
    # After data
    # ------------------------------------------------------------------
    a = after.values.flatten()

    a = a[
        np.isfinite(a)
    ]

    # Single row of 4 panels (Before density, Before box, After density, After
    # box) instead of a 2x2 block -- a compact 4x1 strip that reads at a glance
    # and takes up far less vertical space than the old 2x2 layout. Displayed
    # at full tab width (see render_single_figure's width_ratio=None), so the
    # figure is sized to fill that width well rather than looking small with
    # empty space on either side.
    fig, axes = plt.subplots(
        1,
        4,
        figsize=(15, 3.4)
    )
    fs = 9.5

    # ------------------------------------------------------------------
    # Before: Distribution
    # ------------------------------------------------------------------
    if b.size:
        axes[0].hist(
            b,
            bins=60,
            alpha=0.8
        )

        axes[0].set_xscale("log")

    axes[0].set_title(
        "Before: Raw Peak Area\nDistribution",
        fontsize=fs, fontweight="bold", loc="center"
    )

    axes[0].set_xlabel(
        "Peak Area (log scale)", fontsize=fs
    )

    axes[0].set_ylabel(
        "Count", fontsize=fs
    )

    # ------------------------------------------------------------------
    # Before: Box plot
    # ------------------------------------------------------------------
    if b.size:

        try:
            axes[1].boxplot(
                [b],
                tick_labels=["Before"]
            )
        except TypeError:
            # Compatibility with older matplotlib versions.
            axes[1].boxplot(
                [b],
                labels=["Before"]
            )

        axes[1].set_yscale("log")

    axes[1].set_title(
        "Before: Box Plot",
        fontsize=fs, fontweight="bold", loc="center"
    )

    axes[1].set_ylabel(
        "Peak Area (log scale)", fontsize=fs
    )

    # ------------------------------------------------------------------
    # After: Strict Log2 distribution
    # ------------------------------------------------------------------
    if a.size:

        axes[2].hist(
            a,
            bins=60,
            alpha=0.8
        )

    axes[2].set_title(
        "After: Strict Log2\nDistribution",
        fontsize=fs, fontweight="bold", loc="center"
    )

    axes[2].set_xlabel(
        "log2(Peak Area)", fontsize=fs
    )

    axes[2].set_ylabel(
        "Count", fontsize=fs
    )

    # ------------------------------------------------------------------
    # After: Box plot
    # ------------------------------------------------------------------
    if a.size:

        try:
            axes[3].boxplot(
                [a],
                tick_labels=["After"]
            )
        except TypeError:
            # Compatibility with older matplotlib versions.
            axes[3].boxplot(
                [a],
                labels=["After"]
            )

    axes[3].set_title(
        "After: Strict Log2\nBox Plot",
        fontsize=fs, fontweight="bold", loc="center"
    )

    axes[3].set_ylabel(
        "log2(Peak Area)", fontsize=fs
    )

    for ax in axes:
        ax.tick_params(labelsize=fs - 1)

    # ------------------------------------------------------------------
    # Optional sample identifier
    # ------------------------------------------------------------------
    if sample_id:
        fig.suptitle(
            f"Normalization / Log2 Transformation: {sample_id}",
            fontsize=13, fontweight="bold"
        )

    fig.tight_layout()

    return fig
