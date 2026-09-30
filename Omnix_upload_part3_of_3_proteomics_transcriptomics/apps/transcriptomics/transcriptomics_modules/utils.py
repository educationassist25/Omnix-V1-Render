"""
utils.py - Core data loading, validation, and helper utilities for the Bulk
RNA-seq Analysis App.

Expected input formats
-----------------------
Gene expression / counts matrix (CSV/XLSX):
    Gene, Sample1, Sample2, QC1, QC2, ...
    (rows = genes/features, columns = samples)

Metadata table (CSV/XLSX):
    Sample, Group, IsQC, Batch
    - Sample   : must match a column name in the gene expression matrix
    - Group    : experimental group / condition label
    - IsQC     : True/False (or 1/0) flag for QC samples
    - Batch    : optional batch identifier (used for batch-specific normalization)
"""

import io
import numpy as np
import pandas as pd
import matplotlib
from matplotlib.colors import LinearSegmentedColormap


REQUIRED_METADATA_COLS = ["Sample", "Group"]
GENE_COLUMN_ALIASES = {"gene", "gene name", "gene_name", "genesymbol", "gene symbol",
                        "gene names", "genes"}


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
    validate_gene_matrix/validate_metadata with a clearer structural error).
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


def validate_gene_matrix(df: pd.DataFrame) -> pd.DataFrame:
    """
    Validate a raw gene expression / counts matrix. First column must be gene
    identifiers. An optional SECOND column named 'Gene Symbol' (or a close
    variant, see GENE_COLUMN_ALIASES) is recognized as an identifier -> gene
    symbol mapping and stashed in the returned DataFrame's `.attrs['gene_map']`
    (a {Gene: Symbol} dict) rather than kept as a data column -- every
    downstream numeric operation on the matrix (cleaning, QC, normalization,
    stats, PCA, heatmap, ...) sees only the sample columns.
    """
    if df.shape[1] < 2:
        raise DataValidationError("Gene expression matrix must have a gene column plus at least one sample column.")
    first_col = df.columns[0]
    df = df.rename(columns={first_col: "Gene"})
    df["Gene"] = df["Gene"].astype(str)
    if df["Gene"].duplicated().any():
        dupes = df["Gene"][df["Gene"].duplicated()].unique().tolist()
        raise DataValidationError(f"Duplicate gene identifiers found: {dupes[:5]}...")

    gene_map = None
    if df.shape[1] >= 3 and str(df.columns[1]).strip().lower() in GENE_COLUMN_ALIASES:
        gene_col = df.columns[1]
        gene_map = dict(zip(df["Gene"], df[gene_col].astype(str)))
        df = df.drop(columns=[gene_col])

    sample_cols = df.columns[1:]
    for c in sample_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    out = df.set_index("Gene")
    if gene_map is not None:
        out.attrs["gene_map"] = gene_map
    return out


def validate_metadata(meta: pd.DataFrame, sample_cols) -> pd.DataFrame:
    """Validate metadata table and align it to the sample columns present in the gene matrix."""
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
        raise DataValidationError(f"Samples present in gene matrix but missing from metadata: {missing_samples}")
    return meta.loc[list(sample_cols)]


def split_qc_and_samples(gene_df: pd.DataFrame, meta: pd.DataFrame):
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


def format_for_display(df: pd.DataFrame, sci_threshold: float = 1e-6, decimals: int = 6):
    """
    Format a DataFrame for on-screen display: any numeric value that would need
    more than `decimals` (default 6) decimal places to show a significant digit
    -- i.e. abs(value) < 1e-6 -- is rendered in scientific notation instead of
    a long string of zeros. Everything else is shown with up to 6 decimals,
    trailing zeros trimmed.

    Returns a pandas Styler, so Streamlit still renders a sortable, interactive
    table and the UNDERLYING VALUES ARE UNCHANGED -- this only affects display.
    CSV downloads keep full precision, since they use the raw DataFrame.

    P-values and FDRs in RNA-seq results routinely reach 1e-14 or smaller, where
    a fixed-decimal view would show nothing but "0.000000" and silently destroy
    the ranking information the user is looking at.
    """
    numeric_cols = df.select_dtypes(include=[np.number]).columns

    def _fmt(v):
        if pd.isna(v):
            return ""
        if isinstance(v, (bool, np.bool_)):
            return str(bool(v))
        try:
            fv = float(v)
        except (TypeError, ValueError):
            return str(v)
        if fv == 0:
            return "0"
        if abs(fv) < sci_threshold or abs(fv) >= 1e7:
            return f"{fv:.3e}"
        if float(fv).is_integer() and abs(fv) < 1e7:
            return f"{int(fv)}"
        return f"{fv:.{decimals}f}".rstrip("0").rstrip(".")

    try:
        return df.style.format({c: _fmt for c in numeric_cols})
    except Exception:
        # If styling fails for any reason, fall back to the plain frame rather
        # than breaking the page.
        return df


def significance_flag(pval, fdr, threshold: float):
    """
    Shared "Significant" boolean rule used everywhere a p-value/FDR pair is
    turned into a Significant column (two-group Welch t-test, DESeq2/logCPM
    Wald tests, and ANOVA): FDR < threshold OR raw p-value < threshold,
    against the SAME user-adjustable threshold on both sides.

    This is deliberately an OR, not an AND: the earlier "AND" version
    silently required a gene to clear both cutoffs at once, which is not
    what "use either FDR<0.05 or p<0.05" means and could make the flag
    stricter than either cutoff alone. NaN on either side (e.g. a gene with
    too few non-missing values for a test) never counts as satisfying that
    side -- pandas' `<` already evaluates to False against NaN, but this is
    made explicit with `.notna()` so the behavior doesn't depend on that
    implicit pandas quirk.

    pval, fdr : array-like / pd.Series (same index/order)
    threshold : float, the single shared cutoff applied to both p-value and FDR
    """
    pval = pd.Series(pval) if not isinstance(pval, pd.Series) else pval
    fdr = pd.Series(fdr) if not isinstance(fdr, pd.Series) else fdr
    fdr_hit = fdr.notna() & (fdr < threshold)
    pval_hit = pval.notna() & (pval < threshold)
    fdr_hit = fdr_hit.to_numpy() if hasattr(fdr_hit, "to_numpy") else fdr_hit
    pval_hit = pval_hit.to_numpy() if hasattr(pval_hit, "to_numpy") else pval_hit
    return fdr_hit | pval_hit


def apply_significance_filter(df: pd.DataFrame, metric: str, threshold: float, apply_filter: bool,
                               pval_col: str = "p-value", fdr_col: str = "FDR") -> pd.DataFrame:
    """
    Shared "Significance Threshold" filter/flag helper for the Statistics section's
    Two-Group Comparison and ANOVA (>=3 groups) results tables.

    Replaces the earlier FDR-OR-p-value "Significant" flag with a single, explicit
    user choice: pick ONE metric (FDR or raw p-value) and ONE threshold, then either
    just re-flag "Significant" against that single rule, or additionally drop every
    row that doesn't meet it.

    metric : "FDR" or "P-value" -- which column drives both the Significant flag
      and (when apply_filter) the row filter. Any other value is treated as "P-value".
    threshold : the single cutoff compared against the selected metric.
    apply_filter : if True, only rows with Significant == True are returned (i.e.
      selected_metric <= threshold); if False, every row is returned (with the
      Significant column still recomputed against the current metric/threshold).
    pval_col, fdr_col : the actual column names holding the raw p-value and FDR in
      `df` -- callers pass "ANOVA p-value" for the ANOVA table, the default
      "p-value" everywhere else.

    Note the comparison is `<=` (not `<`), per this filter's spec, which is a
    deliberate difference from the strict `<` used elsewhere. Returns a copy;
    never mutates the passed-in DataFrame.
    """
    df = df.copy()
    metric_col = fdr_col if metric == "FDR" else pval_col
    if metric_col not in df.columns:
        raise KeyError(f"apply_significance_filter: column '{metric_col}' not found in DataFrame")
    values = pd.to_numeric(df[metric_col], errors="coerce")
    sig = values.notna() & (values <= threshold)
    df["Significant"] = sig.to_numpy() if hasattr(sig, "to_numpy") else sig
    if apply_filter:
        df = df[df["Significant"]]
    return df


TWO_GROUP_DISPLAY_COLUMNS = ["Gene", "p-value", "FDR", "Log2FC", "Linear_FC", "Significant"]


def two_group_display_view(result: pd.DataFrame) -> pd.DataFrame:
    """
    Build the standardized, display/export-only 6-column view for ANY two-group
    (single-pair) comparison result table -- pyDESeq2 Wald contrasts, the
    logCPM Wald test, or the Welch's t-test fallback -- regardless of which
    module produced it: Gene, p-value, FDR, Log2FC, Linear_FC, Significant, in
    that exact order.

    This is always a fresh copy built from a copy of `result`; it never
    mutates or truncates the passed-in DataFrame, so callers can keep using
    the full result (baseMean, lfcSE, stat, Mean_A/Mean_B, CI bounds,
    Group_A/Group_B, Higher_In, etc.) for Volcano/Biomarker/Heatmap/other
    downstream consumers.
    """
    df = result.copy()
    if df.index.name == "Gene" or "Gene" not in df.columns:
        df = df.reset_index().rename(columns={df.index.name or "index": "Gene"})
    missing = [c for c in TWO_GROUP_DISPLAY_COLUMNS if c not in df.columns]
    if missing:
        raise KeyError(f"two_group_display_view: result is missing required column(s) {missing}")
    return df[TWO_GROUP_DISPLAY_COLUMNS]


def to_download_bytes_csv(df: pd.DataFrame) -> bytes:
    """
    Serialize a DataFrame to CSV bytes for st.download_button, with a UTF-8 BOM
    (utf-8-sig). Gene names routinely contain Greek letters (α, β, Δ),
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


def build_colormap(base_cmap: str = "RdBu_r", custom_colors=None, reverse: bool = False):
    """
    Build a matplotlib colormap either from a predefined name or a custom list of hex
    colors picked via a 2D color picker (creates a smooth gradient through them via
    LinearSegmentedColormap). Shared by the Heatmap and P-P Interaction tabs so both
    offer the same "predefined palette OR pick your own colors" experience.
    """
    if custom_colors:
        cmap = LinearSegmentedColormap.from_list("custom_gradient", custom_colors, N=256)
    else:
        cmap = matplotlib.colormaps.get_cmap(base_cmap)
    if reverse:
        cmap = cmap.reversed()
    return cmap


# ---------------------------------------------------------------------------
# Figure export helpers (shared by every tab that renders a downloadable plot)
# ---------------------------------------------------------------------------
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


# On-screen rendering DPI for every matplotlib figure shown via render_figure()/
# render_figure_grid() below -- independent of download resolution (which the user
# picks separately in render_figure_download()). Streamlit's st.pyplot() defaults to
# width='stretch', so a figure placed in a narrower st.columns() cell is shrunk to
# fit; without a high enough source DPI, that downscaling (combined with a HiDPI/
# retina screen rendering the shrunk image back up) is exactly what produces the
# blurry/hazy look this DPI raises to prevent.
ON_SCREEN_FIGURE_DPI = 160


def render_figure(st_module, fig, dpi: int = ON_SCREEN_FIGURE_DPI):
    """
    Display one matplotlib figure at a crisp, consistent on-screen resolution.
    Sets the figure's own DPI (matplotlib's default is 100) before handing it to
    st.pyplot(), so text, axis labels, legends, and annotations stay sharp even
    when the figure is displayed at reduced width (e.g. inside a column) or on a
    high-density display. Purely a display-quality change -- it does not touch
    the figure's data, size in inches, or any computed values.
    """
    try:
        fig.set_dpi(dpi)
    except Exception:
        pass
    st_module.pyplot(fig, width="stretch")


def _grid_dims(n: int):
    """
    Choose (nrows, ncols) for a compact, visually balanced grid of n figures.

    Fixed cases the app relies on:
      1 -> (1, 1)            single full-width panel
      2 -> (1, 2)            one row
      3 -> (1, 3)            one row
      4 -> (2, 2)            exactly 2x2
      5 -> (2, 3)            3-then-2 (last row has one empty slot)
      6 -> (3, 2)            exactly 3 rows of 2 columns
      7 -> (3, 3)            4-then-3 (last row has two empty slots... see below)
      8 -> (3, 3)            3+3+2
      9 -> (3, 3)            exactly 3x3

    General rule: start from ncols = ceil(sqrt(n)), nrows = ceil(n / ncols), then
    check the "transposed" factorization (swap roles) and keep whichever produces
    fewer empty cells (ties keep the wider layout, i.e. more columns than rows,
    since figures usually read better wide than tall). n<=3 is always one row;
    n==6 is special-cased to (3, 2) per the compact "3 rows x 2 columns" spec.
    """
    if n <= 0:
        return (0, 0)
    if n <= 3:
        return (1, n)
    if n == 4:
        return (2, 2)
    if n == 6:
        return (3, 2)
    if n == 9:
        return (3, 3)

    import math
    ncols_a = math.ceil(math.sqrt(n))
    nrows_a = math.ceil(n / ncols_a)
    empty_a = ncols_a * nrows_a - n

    nrows_b = math.ceil(math.sqrt(n))
    ncols_b = math.ceil(n / nrows_b)
    empty_b = ncols_b * nrows_b - n

    if empty_b < empty_a:
        return (nrows_b, ncols_b)
    return (nrows_a, ncols_a)


def render_figure_grid(st_module, panels: list, dpi: int = ON_SCREEN_FIGURE_DPI, row_layout=None):
    """
    Lay out a set of related figures in a compact, balanced grid that
    automatically reorganizes based on how many panels are given, instead of
    a single long, sparsely-filled row or an ever-taller single column:
      1 figure   -> full width, single panel
      2-3 figures -> all in one row (2 or 3 columns)
      4 figures  -> exactly 2 rows x 2 columns
      6 figures  -> exactly 3 rows x 2 columns
      9 figures  -> exactly 3 rows x 3 columns
      other N    -> the most compact, near-square arrangement with the fewest
                    empty cells (see _grid_dims)

    `row_layout`: optional explicit list of per-row panel counts, e.g. [4, 3]
    for "4 figures in row 1, then 3 in row 2". This is an irregular/pyramid
    layout (each row gets its own `st.columns(row_len)`, so row 1 and row 2
    are independently full-width) rather than a single rectangular ncols x
    nrows grid. Use this when a section needs a specific, hand-picked row
    split (e.g. QC Validation's 7 figures as 4-then-3, or Normalization's 5
    figures as 3-then-2) instead of the automatic near-square layout above.
    `sum(row_layout)` must equal `len(panels)`; if it doesn't, this falls back
    to the automatic layout so a stale row_layout can never silently drop or
    misplace panels.

    Each entry in `panels` is a dict describing ONE figure and everything that
    goes with it, so nothing is dropped when the layout changes:
      fig            (required)  the matplotlib Figure
      title          (optional)  markdown header shown above the figure
      caption        (optional)  st.caption text shown above the figure
      download_name  (optional)  base filename for the download row
      download_key   (required if download_name is set) unique widget key prefix

    All figures are rendered at the same on-screen DPI (see render_figure), so
    the grid looks uniform regardless of how many panels are shown per row.
    """
    def _draw(container, spec):
        if spec.get("title"):
            container.markdown(spec["title"])
        if spec.get("caption"):
            container.caption(spec["caption"])
        render_figure(container, spec["fig"], dpi=dpi)
        if spec.get("download_name"):
            render_figure_download(container, spec["fig"], spec["download_name"],
                                    key_prefix=spec["download_key"])

    n = len(panels)
    if n == 0:
        return
    if n == 1:
        _draw(st_module, panels[0])
        return

    if row_layout and sum(row_layout) == n:
        idx = 0
        for row_len in row_layout:
            if row_len <= 0:
                continue
            cols = st_module.columns(row_len)
            for i in range(row_len):
                _draw(cols[i], panels[idx])
                idx += 1
        return

    nrows, ncols = _grid_dims(n)
    idx = 0
    for _ in range(nrows):
        remaining = n - idx
        row_len = min(ncols, remaining)
        if row_len <= 0:
            break
        cols = st_module.columns(ncols)
        for i in range(row_len):
            _draw(cols[i], panels[idx])
            idx += 1


def align_download_button(container):
    """
    st.selectbox / st.number_input render a visible label line above their input
    box; st.download_button has no such label line, so when placed in a column
    next to a selectbox it renders noticeably higher -- the button and the
    selectbox's actual input boxes end up on different baselines. This inserts
    an INVISIBLE label-sized spacer (same font-size/line-height/margin as a real
    Streamlit widget label, just hidden) immediately before the button, so its
    top edge lines up with the selectbox input box beside it rather than the
    selectbox's label. Used everywhere a Format/DPI selector sits next to a
    download button, so every export row in the app aligns the same way.
    """
    container.markdown(
        "<div style='visibility:hidden; font-size:14px; line-height:1.6; "
        "margin-bottom:0.25rem;'>&nbsp;</div>",
        unsafe_allow_html=True,
    )


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
