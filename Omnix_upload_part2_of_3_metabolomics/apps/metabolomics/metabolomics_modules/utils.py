"""
utils.py - Core data loading, validation, and helper utilities for MetaboAI Pro.

Expected input formats
-----------------------
Peak area matrix (CSV/XLSX):
    Metabolite, Sample1, Sample2, QC1, QC2, ...
    (rows = metabolites/features, columns = samples)

Metadata table (CSV/XLSX):
    Sample, Group, IsQC, Batch
    - Sample   : must match a column name in the peak area matrix
    - Group    : experimental group / condition label
    - IsQC     : True/False (or 1/0) flag for QC samples
    - Batch    : optional batch identifier (used for batch-specific normalization)
"""

import io
import numpy as np
import pandas as pd


REQUIRED_METADATA_COLS = ["Sample", "Group"]


def display_feature_name(name) -> str:
    """
    Strip a multi-dataset 'Dataset Label::Metabolite' combine-step prefix (see
    dataset_manager._combine_by_label) down to just the metabolite name, for
    anything shown to the user (tables, plot labels, downloads). Names with no
    '::' (single-dataset mode) are returned unchanged.
    """
    text = str(name)
    return text.split("::", 1)[1] if "::" in text else text


def display_index(index) -> list:
    """Vectorized form of display_feature_name() for a pandas Index/list of names."""
    return [display_feature_name(v) for v in index]


def with_display_feature_names(df: pd.DataFrame) -> pd.DataFrame:
    """
    Return a copy of df (features x samples) whose row index shows only the
    metabolite name -- never the 'Dataset Label::' combine-step prefix -- for
    every downstream analysis tab (Statistics, PCA, Volcano, Biomarker Discovery,
    Heatmap, Boxplot). If stripping the prefix makes two different features share
    the same name (the same metabolite name appearing in more than one dataset),
    a small ' (2)', ' (3)', ... suffix is added to the repeats -- in that order of
    first appearance -- so every row stays individually selectable; the dataset
    name itself is never used for this. Deterministic: the same input index
    always produces the same output index, so results stay matched up across
    tabs (e.g. Statistics results vs. the Heatmap's own copy of this data).
    """
    if df is None:
        return df
    stripped = display_index(df.index)
    seen = {}
    deduped = []
    for name in stripped:
        if name not in seen:
            seen[name] = 1
            deduped.append(name)
        else:
            seen[name] += 1
            deduped.append(f"{name} ({seen[name]})")
    out = df.copy()
    out.index = deduped
    return out


class DataValidationError(Exception):
    pass


def load_table(file_or_buffer, filename_hint: str = "") -> pd.DataFrame:
    """
    Load a CSV or Excel file into a DataFrame, auto-detecting format and, for CSV,
    auto-detecting text encoding. Many lab instrument/software exports (Excel "CSV"
    saves, older Windows tools) use Windows-1252/Latin-1, not UTF-8 -- e.g. '±' (as
    in 'mean ± SD'), 'µ' (micro), or curly quotes are classic culprits that raise
    UnicodeDecodeError under pandas' UTF-8 default. Falls back through common
    encodings in order; Latin-1 can decode any byte sequence, so this never raises
    UnicodeDecodeError itself (a genuinely corrupt/binary file will instead fail
    validate_peak_matrix/validate_metadata with a clearer structural error).
    """
    name = filename_hint.lower() if filename_hint else getattr(file_or_buffer, "name", "").lower()
    if name.endswith(".xlsx") or name.endswith(".xls"):
        return pd.read_excel(file_or_buffer)

    for encoding in ("utf-8-sig", "utf-8", "cp1252", "latin-1"):
        try:
            if hasattr(file_or_buffer, "seek"):
                file_or_buffer.seek(0)
            return pd.read_csv(file_or_buffer, encoding=encoding)
        except UnicodeDecodeError:
            continue
    # Unreachable in practice (latin-1 accepts every byte value 0-255), but keeps
    # the function's contract honest if that ever changes.
    if hasattr(file_or_buffer, "seek"):
        file_or_buffer.seek(0)
    return pd.read_csv(file_or_buffer)


def validate_peak_matrix(df: pd.DataFrame) -> pd.DataFrame:
    """Validate raw peak area matrix. First column must be metabolite identifiers."""
    if df.shape[1] < 2:
        raise DataValidationError("Peak area matrix must have a metabolite column plus at least one sample column.")
    first_col = df.columns[0]
    df = df.rename(columns={first_col: "Metabolite"})
    df["Metabolite"] = df["Metabolite"].astype(str)
    if df["Metabolite"].duplicated().any():
        dupes = df["Metabolite"][df["Metabolite"].duplicated()].unique().tolist()
        raise DataValidationError(f"Duplicate metabolite names found: {dupes[:5]}...")
    sample_cols = df.columns[1:]
    non_numeric = [c for c in sample_cols if not pd.api.types.is_numeric_dtype(pd.to_numeric(df[c], errors="coerce"))]
    for c in sample_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df.set_index("Metabolite")


def validate_metadata(meta: pd.DataFrame, sample_cols) -> pd.DataFrame:
    """Validate metadata table and align it to the sample columns present in the peak matrix."""
    missing = [c for c in REQUIRED_METADATA_COLS if c not in meta.columns]
    if missing:
        raise DataValidationError(f"Metadata missing required column(s): {missing}")
    if "IsQC" not in meta.columns:
        meta["IsQC"] = False
    else:
        meta["IsQC"] = meta["IsQC"].astype(str).str.lower().isin(["true", "1", "yes", "qc"])
    if "Batch" not in meta.columns:
        meta["Batch"] = "1"
    meta = meta.set_index("Sample")
    missing_samples = [s for s in sample_cols if s not in meta.index]
    if missing_samples:
        raise DataValidationError(f"Samples present in peak matrix but missing from metadata: {missing_samples}")
    return meta.loc[list(sample_cols)]


def split_qc_and_samples(peak_df: pd.DataFrame, meta: pd.DataFrame):
    """Return (qc_columns, sample_columns) based on metadata IsQC flag."""
    qc_cols = meta.index[meta["IsQC"]].tolist()
    sample_cols = meta.index[~meta["IsQC"]].tolist()
    return qc_cols, sample_cols


def get_categorical_metadata_columns(meta: pd.DataFrame, sample_cols=None, max_unique: int = 15) -> list:
    """
    Which metadata columns are usable as a grouping/coloring variable (for PCA,
    Statistics, Heatmap column annotation, Boxplot) -- not just the hardcoded
    'Group' column. A column qualifies if, among biological samples only (QC rows'
    placeholder values like 'QC'/NaN are excluded from this check), it's non-numeric
    (object/category dtype -- covers text labels like Diagnosis, Gender, Treatment,
    Ethnicity) or numeric with few enough distinct values to plausibly be a
    group/category rather than a continuous measurement (excludes Age, Body Weight,
    and similar continuous covariates, which aren't meaningful t-test/ANOVA/boxplot
    groups). 'Group' is always included first if present, for a stable default.
    """
    if sample_cols is None:
        sample_cols = meta.index.tolist()
    bio_meta = meta.loc[meta.index.intersection(sample_cols)]
    candidates = []
    for col in meta.columns:
        if col in ("IsQC",):
            continue
        series = bio_meta[col].dropna()
        if series.empty:
            continue
        n_unique = series.nunique()
        if n_unique < 2 or n_unique > len(series):
            continue
        if pd.api.types.is_numeric_dtype(series):
            if n_unique <= max_unique:
                candidates.append(col)
        else:
            candidates.append(col)
    # stable, predictable ordering: Group first (if present), then the rest as-authored
    ordered = [c for c in ("Group",) if c in candidates]
    ordered += [c for c in candidates if c not in ordered]
    return ordered


def zero_replacement(df: pd.DataFrame, method: str = "min_fraction", fraction: float = 0.5) -> pd.DataFrame:
    """
    Replace zeros / missing values prior to log transform.
    method:
      - 'min_fraction': replace with fraction * (smallest non-zero value in that feature's row)
      - 'global_min': replace with fraction * (smallest non-zero value across whole matrix)
    """
    out = df.copy()
    if method == "min_fraction":
        for idx in out.index:
            row = out.loc[idx]
            nonzero = row[(row > 0) & row.notna()]
            fill_val = nonzero.min() * fraction if len(nonzero) else np.nan
            out.loc[idx] = row.where((row > 0) & row.notna(), fill_val)
    else:  # global_min
        nonzero = out.values[(out.values > 0) & (~np.isnan(out.values))]
        fill_val = nonzero.min() * fraction if nonzero.size else np.nan
        out = out.where((out > 0) & out.notna(), fill_val)
    return out


def _format_float_str(v, decimals: int = 6) -> str:
    """
    Format a single float for on-screen display as a string: values that need
    more than `decimals` (default 6) decimal places to show meaningful
    precision -- in practice, any non-zero value smaller in magnitude than
    10^-decimals, where fixed-point notation would otherwise round it to all
    zeros or hide its real precision -- are shown in scientific notation (e.g.
    "1.234568e-08"). Everything else (0 and ordinary-magnitude numbers) is
    shown as a plain fixed-point string with up to `decimals` decimal places.
    Missing/invalid values become "" (empty). Always returns a string, so a
    formatted column is never a mix of str and float -- mixed dtypes are what
    trip up Streamlit/PyArrow's table serialization.
    """
    if v is None or pd.isna(v):
        return ""
    if np.isinf(v):
        return "inf" if v > 0 else "-inf"
    if v == 0:
        return f"{0:.{decimals}f}"
    if abs(v) < 10 ** (-decimals):
        return f"{v:.{decimals}e}"
    return f"{v:.{decimals}f}"


def format_df_for_display(df: pd.DataFrame, decimals: int = 6):
    """
    Format every float column of a DataFrame for on-screen display, for use
    with st.dataframe(...) wherever the app shows a data table. Values with
    more than `decimals` decimal places of real precision (i.e. very small
    magnitudes that would otherwise display as 0.000000 or lose precision
    under fixed-point rounding) render in scientific notation; everything else
    renders as a plain number with up to `decimals` decimal places. Integer,
    boolean, string, and other non-float columns pass through untouched. The
    returned object is a new DataFrame (the original is never mutated), so
    this only changes what is shown on screen -- it never affects downstream
    computation or CSV/XLSX export, which always use full, unrounded values.
    """
    if df is None or not isinstance(df, pd.DataFrame):
        return df
    out = df.copy()
    for col in out.columns:
        if pd.api.types.is_float_dtype(out[col]):
            out[col] = out[col].map(lambda v: _format_float_str(v, decimals))
    return out


def to_download_bytes_csv(df: pd.DataFrame) -> bytes:
    """
    Serialize a DataFrame to CSV bytes for st.download_button, with a UTF-8 BOM
    (utf-8-sig). Metabolite/lipid names routinely contain Greek letters (α, β, Δ),
    symbols (±, µ, °), or accented Latin characters -- plain 'utf-8' without a BOM
    is valid and round-trips fine in Python/pandas, but Excel (still the most common
    tool users re-open these exports in) assumes the system locale encoding for a
    plain UTF-8 CSV and renders those characters as mojibake unless a BOM marks it
    explicitly as UTF-8. The BOM is a no-op for pandas/any UTF-8-aware reader.
    """
    return df.to_csv().encode("utf-8-sig")


def to_download_bytes_xlsx(sheets: dict) -> bytes:
    """sheets: dict of {sheet_name: DataFrame}"""
    buf = io.BytesIO()
    with pd.ExcelWriter(buf, engine="openpyxl") as writer:
        for name, d in sheets.items():
            safe_name = name[:31]
            d.to_excel(writer, sheet_name=safe_name)
    buf.seek(0)
    return buf.read()
