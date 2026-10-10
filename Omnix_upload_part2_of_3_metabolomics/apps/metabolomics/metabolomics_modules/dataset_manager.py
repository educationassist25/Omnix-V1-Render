"""
dataset_manager.py - Multi-dataset management for combined LC-MS profiling workflows.

Real-world untargeted/targeted metabolomics studies routinely run the SAME biological
samples through multiple analytical methods (e.g., HILIC vs RP chromatography) and/or
ionization modes (positive/negative), since no single method captures the full
metabolome. Each method/mode combination is cleaned, QC'd, and normalized
INDEPENDENTLY (their raw intensity scales and technical characteristics differ), and
only the final normalized (log2) data is combined into one unified feature matrix for
downstream statistics/visualization.

This module provides the combination logic and cross-dataset QC comparison. Per-dataset
cleaning/QC/normalization reuses the existing single-dataset modules (imputation_module,
qc, normalization) — app.py loops over datasets and calls those as before.
"""

import io
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use("Agg")

from metabolomics_modules import qc as _qc
from metabolomics_modules import normalization as _normalization
from metabolomics_modules import stats_analysis as _stats_analysis
from metabolomics_modules import imputation_module as _imputation_module


def make_dataset_entry(label: str, data_type: str) -> dict:
    """Create a fresh, empty dataset entry for the multi-dataset registry."""
    return {
        "label": label,
        "data_type": data_type,
        "raw_df": None,
        "qc_cols": [],
        "sample_cols": [],
        "cleaned_data": None,
        "imputed_data": None,
        "raw_peak_df_qc": None,
        "raw_peak_df_qc_full": None,
        "log2_data": None,
        "log2_constant": None,
        "istd_normalized": None,
        "istd_mapping": None,
        "iqr_normalized": None,
        "log2_only": None,
        "log2_full": None,
        "fig_dist_log2": None,
        "fig_dist_iqr": None,
        "cv_table": None,
        "qc_log_data": None,
        "processing_notes": [],
    }


def fig_to_png_bytes(fig, dpi: int = 300) -> bytes:
    """Serialize a matplotlib figure to PNG bytes, for st.download_button."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
    buf.seek(0)
    return buf.read()


def fig_to_bytes(fig, fmt: str = "png", dpi: int = 300) -> bytes:
    """
    Serialize a matplotlib figure to bytes in the given format, for
    st.download_button. Supports 'png'/'tiff'/'jpg' (raster, dpi matters) and
    'svg'/'pdf' (vector, dpi is ignored by matplotlib but harmless to pass).
    JPEG has no alpha channel, so a transparent figure/axes background is
    flattened to white rather than left to matplotlib's default (black).
    """
    buf = io.BytesIO()
    save_fmt = "jpeg" if fmt == "jpg" else fmt
    save_kwargs = {"format": save_fmt, "dpi": dpi, "bbox_inches": "tight"}
    if save_fmt == "jpeg":
        save_kwargs["facecolor"] = "white"
    fig.savefig(buf, **save_kwargs)
    buf.seek(0)
    return buf.read()


FIGURE_EXPORT_FORMATS = {
    "PNG (raster)": ("png", "image/png"),
    "JPEG (raster)": ("jpg", "image/jpeg"),
    "TIFF (raster, publication)": ("tiff", "image/tiff"),
    "SVG (vector)": ("svg", "image/svg+xml"),
    "PDF (vector)": ("pdf", "application/pdf"),
}


def render_figure_download(st_module, fig, base_filename: str, key_prefix: str,
                            dpi_options=(150, 300, 600), default_dpi_index: int = 1):
    """
    Shared 'download this figure' widget: format selector (PNG/TIFF/SVG/PDF) + DPI
    selector (for raster formats) + a download button, all in one row, all
    vertically aligned to the same baseline. Used everywhere a plot needs a
    high-resolution, multi-format download (QC plots, Combined QC plots, PCA
    plots) so the control looks and behaves identically across tabs. `key_prefix`
    must be unique per call site (dataset id + plot name) to avoid Streamlit
    widget-key collisions when the same plot type is rendered once per dataset.
    """
    c1, c2, c3 = st_module.columns([1.8, 1, 2])
    fmt_label = c1.selectbox(
        "Format", list(FIGURE_EXPORT_FORMATS.keys()), key=f"{key_prefix}_fmt", index=0
    )
    fmt, mime = FIGURE_EXPORT_FORMATS[fmt_label]
    is_raster = fmt in ("png", "tiff", "jpg")
    if is_raster:
        dpi = c2.selectbox("DPI", list(dpi_options), key=f"{key_prefix}_dpi",
                            index=default_dpi_index)
    else:
        dpi = 300
        c2.write("")
        c2.write("")
        c2.caption("(vector — resolution-independent)")
    file_name = f"{base_filename}.{fmt}"
    # Selectboxes render their own label line above the input; the download
    # button has no such label, so a matching blank spacer keeps its top edge
    # level with the Format/DPI inputs instead of sitting visibly higher.
    c3.write("")
    c3.write("")
    c3.download_button(
        f"Download {fmt.upper()}", fig_to_bytes(fig, fmt=fmt, dpi=dpi),
        file_name=file_name, mime=mime, icon="⬇️", key=f"{key_prefix}_dl"
    )


def validate_shared_samples(datasets: dict, meta: pd.DataFrame):
    """
    Check that every dataset's columns are present in the shared metadata index.
    Returns a list of warning strings (empty if everything lines up).
    """
    warnings = []
    if not datasets:
        return warnings
    for ds_id, ds in datasets.items():
        if ds["raw_df"] is None:
            continue
        missing = [c for c in ds["raw_df"].columns if c not in meta.index]
        if missing:
            warnings.append(
                f"Dataset '{ds['label']}': {len(missing)} sample column(s) not found in the "
                f"shared metadata (e.g. {missing[:3]}). Every dataset must use the same sample "
                f"names as the shared metadata table."
            )
    return warnings


def _combine_by_label(datasets: dict, sample_cols: list, getter) -> tuple:
    """
    Shared concatenation helper: for each dataset, pull a per-dataset DataFrame via
    `getter(ds)`, align it to `sample_cols`, prefix its feature (row) names with the
    dataset label (e.g. "Targeted Metabolomics::Glucose") so the same compound name across
    methods/modes is never silently merged, and concatenate across datasets.

    Datasets for which `getter(ds)` returns None are skipped. Returns
    (combined_df_or_None, per_dataset_feature_counts dict).
    """
    pieces = []
    counts = {}
    for ds_id, ds in datasets.items():
        df = getter(ds)
        if df is None:
            continue
        df = df.copy().reindex(columns=sample_cols)  # align columns, NaN-fill any missing sample
        df.index = [f"{ds['label']}::{feat}" for feat in df.index]
        pieces.append(df)
        counts[ds["label"]] = df.shape[0]
    if not pieces:
        return None, counts
    combined = pd.concat(pieces, axis=0)
    return combined, counts


def combine_imputed_datasets(datasets: dict, sample_cols: list) -> tuple:
    """
    Concatenate each dataset's post-cleaning/imputation peak area matrix (raw
    intensity scale, NOT yet normalized) into one unified feature matrix, prefixed by
    dataset label. Falls back to 'cleaned_data' (if imputation wasn't run because there
    was nothing to impute) or the raw biological samples (if Tab 2 was skipped
    entirely for that dataset), so every loaded dataset can still contribute.

    Returns (combined_df_or_None, per_dataset_feature_counts dict).
    """
    def _getter(ds):
        if ds.get("imputed_data") is not None:
            return ds["imputed_data"]
        if ds.get("cleaned_data") is not None:
            return ds["cleaned_data"]
        if ds.get("raw_df") is not None:
            cols = [c for c in ds.get("sample_cols", []) if c in ds["raw_df"].columns]
            return ds["raw_df"][cols] if cols else None
        return None
    return _combine_by_label(datasets, sample_cols, _getter)


def imputation_summary_table(datasets: dict) -> pd.DataFrame:
    """Summary of each dataset's cleaning/imputation status, for the Combine step in Tab 2."""
    rows = []
    for ds_id, ds in datasets.items():
        if ds.get("imputed_data") is not None:
            status, n_features = "Imputed", ds["imputed_data"].shape[0]
        elif ds.get("cleaned_data") is not None:
            status, n_features = "Cleaned (no imputation needed)", ds["cleaned_data"].shape[0]
        elif ds.get("raw_df") is not None:
            status, n_features = "Not cleaned (using raw)", ds["raw_df"].shape[0]
        else:
            status, n_features = "No data", 0
        rows.append({"Dataset": ds["label"], "Analysis Type": ds["data_type"],
                     "Status": status, "Features": n_features})
    return pd.DataFrame(rows)


def combine_normalized_datasets(datasets: dict, sample_cols: list, use_log2: bool) -> tuple:
    """
    Concatenate each dataset's normalized matrix into one unified feature matrix.

    use_log2=True  -> uses each dataset's 'log2_data' (final matrix for that dataset's
                       pipeline -- already fully normalized AND log2-transformed,
                       regardless of which order those two steps ran in).
    use_log2=False -> uses each dataset's matrix from immediately before its final
                       pipeline step: 'istd_normalized' for targeted datasets,
                       'iqr_normalized' for untargeted datasets whose pipeline runs
                       IQR-normalization before log2, or 'log2_only' for untargeted
                       datasets whose pipeline runs log2 before IQR-normalization
                       (whichever is present).

    Metabolite names are prefixed with the dataset label (e.g. "Targeted Metabolomics::Glucose")
    to keep them distinct across methods/modes -- the same compound name can appear in
    multiple methods/ionization modes as a genuinely different measurement (different
    adduct, different chromatographic behavior), so collapsing them by name would
    silently merge unrelated measurements.

    Only datasets with the relevant data present are included. All included datasets
    are aligned to the same `sample_cols` (biological samples only); any dataset
    missing a sample gets NaN for that sample's values in the combined matrix.

    Returns (combined_df_or_None, per_dataset_feature_counts dict).
    """
    def _getter(ds):
        if use_log2:
            return ds.get("log2_data")
        if ds.get("istd_normalized") is not None:
            return ds["istd_normalized"]
        if ds.get("iqr_normalized") is not None:
            return ds["iqr_normalized"]
        return ds.get("log2_only")
    return _combine_by_label(datasets, sample_cols, _getter)


def combine_log2_datasets(datasets: dict, sample_cols: list) -> tuple:
    """Backward-compatible alias for combine_normalized_datasets(..., use_log2=True)."""
    return combine_normalized_datasets(datasets, sample_cols, use_log2=True)


def qc_comparison_table(datasets: dict) -> pd.DataFrame:
    """Side-by-side QC summary across all datasets that have a computed CV table."""
    rows = []
    for ds_id, ds in datasets.items():
        cv_table = ds.get("cv_table")
        if cv_table is None:
            continue
        n_total = len(cv_table)
        n_acceptable = (cv_table["Quality"] == "Acceptable").sum()
        rows.append({
            "Dataset": ds["label"],
            "Total Features": n_total,
            "Acceptable (CV≤20%)": n_acceptable,
            "% Acceptable": round(100 * n_acceptable / n_total, 1) if n_total else 0,
            "Median CV%": round(cv_table["CV(%)"].median(), 1),
        })
    return pd.DataFrame(rows)


def qc_comparison_plot(datasets: dict):
    """Side-by-side bar chart comparing % acceptable-CV features across datasets."""
    table = qc_comparison_table(datasets)
    if table.empty:
        return None
    fig, ax = plt.subplots(figsize=(max(5, 1.0 * len(table)), 3.6))
    colors = ["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B2", "#937860", "#DA8BC3", "#8C8C8C",
              "#CCB974", "#64B5CD", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02",
              "#A6761D", "#666666", "#8DD3C7", "#FB8072"]
    bars = ax.bar(table["Dataset"], table["% Acceptable"], color=[colors[i % len(colors)] for i in range(len(table))])
    ax.set_ylabel("% Features with\nAcceptable CV (≤20%)")
    ax.set_title("QC Comparison Across Datasets", fontsize=11, fontweight="bold", loc="center")
    ax.set_ylim(0, 112)  # headroom so a 100% value label stays inside the plot
    for bar, val in zip(bars, table["% Acceptable"]):
        ax.annotate(f"{val}%", (bar.get_x() + bar.get_width() / 2, val), xytext=(0, 3),
                    textcoords="offset points", ha="center", va="bottom", fontsize=9)
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right")
    fig.tight_layout()
    return fig


def combine_cv_tables(datasets: dict) -> pd.DataFrame:
    """
    Concatenate each dataset's per-feature CV table (Mean, SD, CV(%), Quality) into one
    combined table, with feature names prefixed by dataset label. Feeds the Combined QC
    section's CV Distribution and Feature-Counts-by-Quality plots (reusing qc.py's
    existing plot functions directly on this combined table).
    """
    pieces = []
    for ds_id, ds in datasets.items():
        cv_table = ds.get("cv_table")
        if cv_table is None:
            continue
        t = cv_table.copy()
        t.index = [f"{ds['label']}::{feat}" for feat in t.index]
        pieces.append(t)
    if not pieces:
        return None
    return pd.concat(pieces, axis=0).sort_values("CV(%)")


def combine_qc_log_data(datasets: dict, qc_cols: list) -> pd.DataFrame:
    """
    Concatenate each dataset's QC-only log2 matrix (features x QC replicates) into one
    combined matrix, prefixed by dataset label. Feeds the Combined QC section's
    Sample Correlation Matrix (reusing qc.sample_correlation_matrix directly).
    """
    def _getter(ds):
        return ds.get("qc_log_data")
    combined, _ = _combine_by_label(datasets, qc_cols, _getter)
    return combined


def normalization_summary_table(datasets: dict) -> pd.DataFrame:
    """Summary of each dataset's normalization status, for the Combine step."""
    rows = []
    for ds_id, ds in datasets.items():
        status = "Complete" if ds.get("log2_data") is not None else "Not yet normalized"
        n_features = ds["log2_data"].shape[0] if ds.get("log2_data") is not None else 0
        rows.append({
            "Dataset": ds["label"], "Analysis Type": ds["data_type"],
            "Status": status, "Features": n_features,
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Pipeline status dashboard — shown above the dataset dropdown in Tabs 2-4 so
# it's always clear, at a glance, which datasets still need attention.
# ---------------------------------------------------------------------------
def dataset_status(ds: dict) -> dict:
    """One-row-per-dataset status snapshot across the Cleaning/QC/Normalization stages."""
    if ds.get("imputed_data") is not None:
        clean_status = "✅ Done"
    elif ds.get("cleaned_data") is not None:
        clean_status = "🟡 Filtered only"
    else:
        clean_status = "⚪ Not started"

    if ds.get("raw_peak_df_qc") is not None:
        qc_status = "✅ Confirmed"
    elif ds.get("cv_table") is not None:
        qc_status = "🟡 Reviewed, not confirmed"
    else:
        qc_status = "⚪ Not started"

    if ds.get("log2_data") is not None:
        norm_status = "✅ Done"
    else:
        norm_status = "⚪ Not started"

    return {
        "Dataset": ds["label"], "Analysis Type": ds["data_type"],
        "1. Cleaning": clean_status, "2. QC": qc_status, "3. Normalization": norm_status,
    }


def dataset_status_table(datasets: dict) -> pd.DataFrame:
    """Status dashboard across all loaded datasets, for display above the dataset dropdown."""
    return pd.DataFrame([dataset_status(ds) for ds in datasets.values()])


def all_datasets_ready(datasets: dict) -> bool:
    """True once every loaded dataset has completed normalization (log2_data present)."""
    return bool(datasets) and all(ds.get("log2_data") is not None for ds in datasets.values())


# ---------------------------------------------------------------------------
# One-click "recommended defaults" pipeline — lets a demo (or any dataset) reach
# downstream analysis (PCA/Statistics/...) without manually stepping through
# Tabs 2-4 for every dataset. Uses the same underlying functions the manual UI
# calls, just with sensible defaults pre-selected.
# ---------------------------------------------------------------------------
def auto_process_dataset(ds: dict, max_pct_missing: float = 30, cv_filter: bool = True,
                          imputation_method: str = "Half-Minimum (LOD/2)", row_annotations=None) -> list:
    """
    Run the recommended default pipeline on ONE dataset entry, mutating it in place:
    missingness filter -> imputation (if needed) -> QC CV computation/filtering ->
    ISTD (targeted) or Median-IQR (untargeted) normalization -> Log2 transform.
    Returns a list of human-readable note strings describing what was done.

    cv_filter (CV>20% feature exclusion) is only actually applied for TARGETED
    Metabolomics/Lipidomics. For Untargeted Metabolomics/Lipidomics, this
    parameter is ignored and no feature is ever excluded from Statistics based
    on CV -- matching the manual Normalization tab's own default (its CV
    filter checkbox defaults to unchecked) and external reference tools, since
    a smaller/different feature set changes the Benjamini-Hochberg FDR
    denominator (n) and produces different FDR values even when every
    individual p-value is unchanged. CV is still computed and shown either
    way; it's just not used to drop untargeted rows automatically.

    Internal standards (ISTD rows: named "ISTD..." or listed in the row annotations' ISTD
    column) are never removed by the missingness or CV filters. Each targeted dataset is
    normalized with its ONE internal standard (normalization.default_istd: the ISTD named in
    the row annotations' ISTD column, else the first ISTD row); untargeted datasets simply
    drop ISTD rows, since an internal standard is not a biological feature.
    """
    notes = []
    sample_cols = ds["sample_cols"]
    qc_cols = ds["qc_cols"]
    raw_df = ds["raw_df"]
    is_targeted = "Targeted" in ds["data_type"]
    istds = _normalization.detect_istds(raw_df.index, row_annotations, ds.get("label"), ds["data_type"])

    # 1. Cleaning + imputation
    source_df = raw_df[sample_cols]
    miss_table = _imputation_module.compute_missingness(source_df)
    cleaned = _imputation_module.filter_by_missingness(source_df, miss_table, max_pct_missing)
    cleaned = _normalization.keep_istds(cleaned, source_df, istds)
    ds["cleaned_data"] = cleaned
    n_missing = _imputation_module.count_missing(cleaned)
    if n_missing > 0:
        imputed = _imputation_module.impute(cleaned, imputation_method)
        ds["imputation_method_used"] = imputation_method
    else:
        imputed = cleaned
    ds["imputed_data"] = imputed
    notes.append(
        f"{ds['label']}: cleaned ({cleaned.shape[0]}/{source_df.shape[0]} features kept at "
        f"{max_pct_missing}% missingness threshold); {n_missing} missing values imputed via "
        f"{imputation_method}."
    )

    # 2. QC
    #    cv_filter is only actually applied for TARGETED assays here (unchanged
    #    default behavior). For UNTARGETED Metabolomics/Lipidomics, Quick Start
    #    now matches the manual Normalization tab's own default (its "Filter
    #    out 'Variable' features (CV>20%)" checkbox defaults to unchecked) and
    #    external reference statistics tools: no feature is excluded from
    #    Statistics/FDR based on CV, since a smaller/different feature set
    #    changes the Benjamini-Hochberg denominator (n) and produces different
    #    FDR values even when every individual p-value is identical. CV is
    #    still computed and shown (cv_table), just not used to drop rows.
    effective_cv_filter = cv_filter if is_targeted else False
    if len(qc_cols) >= 2:
        cv_table = _qc.calculate_cv(raw_df[qc_cols])
        ds["cv_table"] = cv_table
        qc_log = _stats_analysis.strict_log2(raw_df[qc_cols])
        ds["qc_log_data"] = qc_log
        working_df = imputed.copy()
        # Keep the full (pre-CV-filter) feature table too: even when cv_filter
        # removes "Variable" features from the table used downstream, the
        # per-sample Median-IQR statistics below are computed from the
        # complete detected-feature background (this full table), not just
        # the reduced/exported subset -- matching the recommended practice
        # of basing median_j/IQR_j on every feature detected in that sample.
        ds["raw_peak_df_qc_full"] = working_df
        if effective_cv_filter:
            keep = cv_table.index[cv_table["Quality"] == "Acceptable"]
            working_df = _normalization.keep_istds(working_df.loc[working_df.index.intersection(keep)],
                                                   imputed, istds)
        ds["raw_peak_df_qc"] = working_df
        n_accept = int((cv_table["Quality"] == "Acceptable").sum())
        notes.append(
            f"{ds['label']}: QC — {n_accept}/{len(cv_table)} features Acceptable (CV≤20%)"
            + (", variable features removed" if effective_cv_filter
               else ", no CV-based filtering applied (untargeted default: all detected "
                    "features proceed to Statistics)" if not is_targeted else "")
            + "."
        )
    else:
        ds["raw_peak_df_qc"] = imputed
        ds["raw_peak_df_qc_full"] = imputed
        notes.append(f"{ds['label']}: fewer than 2 QC replicates available — QC filtering skipped.")

    # 3. Normalization + Log2
    #    Targeted (Metabolomics/Lipidomics): ISTD normalization -> strict Log2.
    #    Untargeted (Metabolomics/Lipidomics): strict Log2 -> Median-IQR, matching the
    #    manual pipeline in the Normalization tab. Running Log2 before the IQR scaling
    #    step (rather than after) avoids the large number of missing values that
    #    Log2-of-negative-numbers would otherwise create, since Median-IQR scaling
    #    centers values around zero; avoid_nan=True additionally keeps a zero-spread
    #    row/column (IQR = 0) centered at 0 instead of becoming missing.
    #    Default is SAMPLE-based (axis="sample"):
    #        normalized[i, j] = (log2(raw[i, j]) - median_j) / IQR_j
    #    where median_j and IQR_j are the median and interquartile range of sample j's
    #    log2 values across all metabolites i -- correcting for total-intensity/loading
    #    differences between samples, matching the Normalization tab's default.
    #    median_j/IQR_j are computed from the FULL detected-feature table at this stage
    #    (raw_peak_df_qc_full, i.e. before the optional CV filter removed any "Variable"
    #    features), not just the reduced table that cv_filter may have exported -- so
    #    turning cv_filter on never changes the per-sample scaling statistics, only which
    #    rows are kept in the final normalized/exported table.
    base_df = ds["raw_peak_df_qc"]
    full_df = ds.get("raw_peak_df_qc_full")
    if full_df is None:
        full_df = base_df
    if is_targeted:
        istd = _normalization.default_istd(base_df.index, row_annotations, ds.get("label"), ds["data_type"])
        if istd is None:
            raise ValueError(
                f"{ds['label']}: no internal standard found. Targeted datasets need an ISTD row (named "
                "'ISTD_...' or listed in an 'ISTD' column of the row annotations); normalize this "
                "dataset on the Normalization page instead.")
        norm, mapping = _normalization.istd_normalize_mapped(
            base_df, [istd], row_annotations, label=ds.get("label"), method=ds["data_type"],
            drop_features=[i for i in istds if i != istd])
        ds["istd_normalized"] = norm
        ds["istd_mapping"] = mapping
        log2_df = _stats_analysis.strict_log2(norm)
        notes.append(f"{ds['label']}: ISTD-normalized using '{istd}' (all {len(norm)} features divided by it, "
                     "sample by sample), then Log2-transformed.")
    else:
        if istds:
            base_df = _normalization.drop_istds(base_df, istds)
            full_df = _normalization.drop_istds(full_df, istds)
            notes.append(f"{ds['label']}: {len(istds)} internal standard row(s) excluded (not biological "
                         "features).")
        log2_only = _stats_analysis.strict_log2(base_df)
        ds["log2_only"] = log2_only
        log2_full = _stats_analysis.strict_log2(full_df)
        ds["log2_full"] = log2_full
        log2_df = _normalization.iqr_normalize(
            log2_only, axis="sample", avoid_nan=True, stats_source=log2_full
        )
        notes.append(f"{ds['label']}: Log2-transformed, then Median-IQR normalized (sample-based: "
                     f"normalized[i,j] = (log2(raw[i,j]) - median_j) / IQR_j, median_j/IQR_j from "
                     f"all {log2_full.shape[0]} detected features per sample).")
    ds["log2_data"] = log2_df
    ds["log2_constant"] = 0.0

    return notes


def auto_process_all_datasets(datasets: dict, max_pct_missing: float = 30, cv_filter: bool = True,
                               imputation_method: str = "Half-Minimum (LOD/2)", row_annotations=None) -> list:
    """Run auto_process_dataset() on every loaded dataset. Returns the combined note list."""
    notes = []
    for ds in datasets.values():
        notes.extend(auto_process_dataset(
            ds, max_pct_missing=max_pct_missing, cv_filter=cv_filter,
            imputation_method=imputation_method, row_annotations=row_annotations
        ))
    return notes
