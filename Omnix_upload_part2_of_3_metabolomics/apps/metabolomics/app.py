"""
Omnix Metabolomics — Automated LC-MS Metabolomics Statistical Analysis and Reporting Platform
Streamlit application entry point.
"""

import os
import re
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as _plt
import streamlit as st

sys.path.insert(0, os.path.dirname(__file__))
from metabolomics_modules import utils, qc, normalization, imputation_module, dataset_manager, stats_analysis, pca_module, volcano, biomarker, heatmap_module, boxplot_module
from metabolomics_modules import pathway_analysis
from metabolomics_modules import rewiring, rewiring_ui

# Render every matplotlib figure at a high enough DPI for crisp, non-blurry display on
# HiDPI/retina screens (figure.dpi controls on-screen sharpness; savefig.dpi is the
# fallback for any export path that doesn't pass its own explicit dpi -- most already
# do, at 300+). This only affects rendering resolution, never the underlying data,
# statistics, or figure layout/size in inches.
_plt.rcParams["figure.dpi"] = 150
_plt.rcParams["savefig.dpi"] = 300
_plt.rcParams["savefig.bbox"] = "tight"

# ---------------------------------------------------------------------------
# Which app is running. metabolomics_app.py and lipidomics_app.py run this same file with
# APP_FAMILY set, so each app offers only its own two analysis types (Untargeted / Targeted);
# every page and analysis is otherwise identical. Running app.py directly keeps all four.
# ---------------------------------------------------------------------------
APP_FAMILY = globals().get("APP_FAMILY", "")          # "Metabolomics", "Lipidomics" or ""
ALL_DATA_TYPES = ["Untargeted Metabolomics", "Targeted Metabolomics", "Untargeted Lipidomics", "Targeted Lipidomics"]
DATA_TYPES = [t for t in ALL_DATA_TYPES if not APP_FAMILY or t.endswith(APP_FAMILY)]
APP_NAME = APP_FAMILY or "Metabolomics"

_GROUPS_DESC = "24 biological samples across 4 groups (Control/Mild/Moderate/Severe), 6 QC replicates"
# Standard demo for each analysis type: (matrix CSV, row-annotation CSV, Method filter, description)
STANDARD_DEMOS = {
    "Untargeted Metabolomics": (
        "sample_peak_area_matrix_untargeted.csv", "metabolite_row_annotations_untargeted.csv", None,
        f"300 metabolites (real names with HMDB IDs; no internal standards), {_GROUPS_DESC}"),
    "Targeted Metabolomics": (
        "sample_peak_area_matrix_targeted.csv", "metabolite_row_annotations_targeted.csv", None,
        f"59 metabolites + 1 internal standard (ISTD_D4-Alanine; every metabolite is linked to it in the "
        f"row annotations' ISTD column), {_GROUPS_DESC}"),
    "Untargeted Lipidomics": (
        "sample_peak_area_matrix_multimethod_untargeted_lipidomics.csv",
        "metabolite_row_annotations_multimethod.csv", "Untargeted Lipidomics",
        f"90 lipids (no internal standards), {_GROUPS_DESC}"),
    "Targeted Lipidomics": (
        "sample_peak_area_matrix_multimethod_targeted_lipidomics.csv",
        "metabolite_row_annotations_multimethod.csv", "Targeted Lipidomics",
        f"39 lipids + 1 lipid internal standard (ISTD_PC(15:0/18:1-d7); every lipid is linked to it in the "
        f"row annotations' ISTD column), {_GROUPS_DESC}"),
}

# Inside the Omnix portal (omnix_app.py) the portal sets the page title and header; run on its own,
# the app sets its own.
OMNIX = globals().get("OMNIX", False)
if not OMNIX:
    st.set_page_config(page_title=f"Omnix · {APP_NAME}", layout="wide", page_icon="🫧" if APP_FAMILY == "Lipidomics" else "🧪")

st.markdown(
    """
    <style>
    /* Sidebar section ordering: title/caption above the group box, no gap between them. */
    [data-testid="stSidebarContent"] { display: flex; flex-direction: column; }
    [data-testid="stSidebarHeader"] { order: 0; }
    [data-testid="stSidebarUserContent"] { order: 1; padding-bottom: 0px !important; }

    /* ---- Sidebar: group-name-only box, one group per row, professional nav-list look ---- */
    div[data-testid="stVerticalBlock"][class*="st-key-nav_group_box"] {
        background: #e6f4f2;
        border: 1px solid #bfe3dc;
        border-radius: 10px;
        padding: 8px !important;
        gap: 4px !important;
    }
    [class*="st-key-nav_group_box"] [data-testid="stButton"] button {
        background: transparent !important;
        border: 1px solid transparent !important;
        border-left: 3px solid #bfe3dc !important;
        box-shadow: none !important;
        border-radius: 6px !important;
        justify-content: flex-start !important;
        text-align: left !important;
        white-space: nowrap !important;
        overflow: visible !important;
        padding: 10px 12px !important;
        width: 100%;
        transition: background-color 0.12s ease-in-out, color 0.12s ease-in-out,
                    border-left-color 0.12s ease-in-out;
    }
    /* The button's visible label is a nested <p> with its own font-size from
       Streamlit's base styles -- it doesn't inherit from the button, so every
       text property has to be set here directly, not on the button itself. */
    [class*="st-key-nav_group_box"] [data-testid="stButton"] button p {
        color: #00695c !important;
        font-weight: 700 !important;
        font-size: 18px !important;
        line-height: 1.15 !important;
        letter-spacing: 0.1px;
        white-space: nowrap !important;
        word-break: keep-all !important;
        overflow-wrap: normal !important;
        hyphens: none !important;
    }
    [class*="st-key-nav_group_box"] [data-testid="stButton"] button:hover {
        background: rgba(0, 105, 92, 0.10) !important;
        border-left-color: #00897b !important;
    }
    [class*="st-key-nav_group_box"] [data-testid="stButton"] button:hover p {
        color: #00695c !important;
    }
    [class*="st-key-nav_group_box"] [data-testid="stButton"] button:focus:not(:active) p {
        color: #00695c !important;
    }

    /* ---- Main page: row of page names belonging to the selected group ---- */
    [class*="st-key-nav_page_row"] [data-testid="stHorizontalBlock"] {
        display: flex;
        flex-wrap: wrap;
        row-gap: 8px;
        column-gap: 10px;
    }
    [class*="st-key-nav_page_row"] [data-testid="stColumn"] {
        width: auto !important;
        flex: 0 0 auto !important;
        min-width: 0 !important;
    }
    [class*="st-key-nav_page_row"] [data-testid="stButton"] button {
        background: #f8fafc !important;
        border: 1px solid #e2e8f0 !important;
        box-shadow: none !important;
        border-radius: 8px !important;
        white-space: nowrap !important;
        overflow: visible !important;
        padding: 10px 22px !important;
        transition: background-color 0.12s ease-in-out, border-color 0.12s ease-in-out;
    }
    [class*="st-key-nav_page_row"] [data-testid="stButton"] button p {
        color: #334155 !important;
        font-weight: 600 !important;
        font-size: 22px !important;
        line-height: 1.15 !important;
        white-space: nowrap !important;
        word-break: keep-all !important;
        overflow-wrap: normal !important;
        hyphens: none !important;
    }
    [class*="st-key-nav_page_row"] [data-testid="stButton"] button:hover {
        background: #eef2f7 !important;
        border-color: #cbd5e1 !important;
        text-decoration: none !important;
    }
    [class*="st-key-nav_page_row"] [data-testid="stButton"] button:focus:not(:active) p {
        color: #334155 !important;
    }

    /* ---- Sidebar panel: subtle gradient so it reads as a distinct "app shell" rail,
       the way most professional dashboards (Notion, Linear, Stripe) separate nav from
       content instead of using the same flat white everywhere. ---- */
    [data-testid="stSidebar"] {
        background: linear-gradient(180deg, #f8fafc 0%, #e6f4f2 100%) !important;
        border-right: 1px solid #e2e8f0;
    }
    [data-testid="stSidebar"] [data-testid="stSidebarHeader"] h1,
    [data-testid="stSidebar"] .stSidebar h1 {
        color: #00695c !important;
    }

    /* ---- Global action buttons: every st.button in the main content area gets a
       bold, colorful, "clearly clickable" gradient pill with a hover lift — instead
       of blending into the page like a plain outlined control. Higher-specificity
       nav rules above still win inside the nav rails, so this only reaches ordinary
       page buttons (Run, Apply, Confirm, Generate, Discover, etc.). ---- */
    [data-testid="stButton"] button {
        background: linear-gradient(135deg, #00897b 0%, #00695c 100%) !important;
        border: none !important;
        border-radius: 8px !important;
        padding: 0.5rem 1.1rem !important;
        box-shadow: 0 2px 6px rgba(0, 105, 92, 0.35) !important;
        transition: transform 0.12s ease-in-out, box-shadow 0.12s ease-in-out,
                    background 0.12s ease-in-out !important;
        height: auto !important;
        min-height: 2.5rem !important;
        white-space: normal !important;
    }
    [data-testid="stButton"] button p {
        color: #ffffff !important;
        font-weight: 600 !important;
        white-space: normal !important;
        overflow-wrap: break-word !important;
        word-break: break-word !important;
        line-height: 1.25 !important;
    }
    [data-testid="stButton"] button:hover {
        background: linear-gradient(135deg, #00695c 0%, #004d40 100%) !important;
        box-shadow: 0 4px 12px rgba(0, 105, 92, 0.45) !important;
        transform: translateY(-1px) !important;
    }
    [data-testid="stButton"] button:hover p { color: #ffffff !important; }
    [data-testid="stButton"] button:active { transform: translateY(0) !important; }
    [data-testid="stButton"] button:focus:not(:active) p { color: #ffffff !important; }
    [data-testid="stButton"] button:disabled {
        background: #cbd5e1 !important;
        box-shadow: none !important;
        transform: none !important;
    }
    [data-testid="stButton"] button:disabled p { color: #64748b !important; }

    /* ---- Download buttons: a distinct teal/emerald accent so "get a file" reads as
       visually different from "run/apply" at a glance, across every tab. ---- */
    [data-testid="stDownloadButton"] button {
        background: linear-gradient(135deg, #f57c00 0%, #e65100 100%) !important;
        border: none !important;
        border-radius: 8px !important;
        padding: 0.5rem 1.1rem !important;
        box-shadow: 0 2px 6px rgba(230, 81, 0, 0.35) !important;
        transition: transform 0.12s ease-in-out, box-shadow 0.12s ease-in-out,
                    background 0.12s ease-in-out !important;
        height: auto !important;
        min-height: 2.5rem !important;
        white-space: normal !important;
    }
    [data-testid="stDownloadButton"] button p {
        color: #ffffff !important;
        font-weight: 600 !important;
        white-space: normal !important;
        overflow-wrap: break-word !important;
        word-break: break-word !important;
        line-height: 1.25 !important;
    }
    [data-testid="stDownloadButton"] button:hover {
        background: linear-gradient(135deg, #e65100 0%, #bf360c 100%) !important;
        box-shadow: 0 4px 12px rgba(230, 81, 0, 0.45) !important;
        transform: translateY(-1px) !important;
    }
    [data-testid="stDownloadButton"] button:hover p { color: #ffffff !important; }
    [data-testid="stDownloadButton"] button:active { transform: translateY(0) !important; }
    </style>
    """,
    unsafe_allow_html=True,
)

# ---------------------------------------------------------------------------
# Session state initialization
# ---------------------------------------------------------------------------
for key, default in [
    ("raw_peak_df", None), ("meta", None), ("data_type", DATA_TYPES[0]),
    ("row_annotations", None),
    ("heatmap_annotation_colors_applied", None), ("heatmap_colors_reset_pending", False),
    ("qc_cols", []), ("sample_cols", []), ("istd_normalized", None), ("istd_mapping", None),
    ("iqr_normalized", None), ("log2_data", None), ("log2_only", None), ("log2_constant", 0.0),
    ("stats_result", None), ("stats_result_groups", None),
    ("anova_result", None), ("posthoc_result", None), ("anova_groups_used", None),
    ("stats_result_group_col", None),
    ("processing_notes", []), ("qc_figs", {}), ("viz_figs", {}), ("heatmap_fig", None), ("heatmap_zscore", None),
    ("volcano_fig", None), ("volcano_annotated", None), ("volcano_settings", None),
    ("boxplot_fig", None), ("boxplot_stats", None),
    ("missingness_table", None), ("cleaned_data", None), ("imputed_data", None), ("imputation_method_used", None),
    ("data_mode", "single"), ("datasets", {}), ("combined_done", False),
    ("raw_peak_df_qc", None), ("raw_peak_df_qc_full", None), ("log2_full", None),
    ("cv_table", None), ("qc_log_data", None),
    ("combined_normalized_prelog2", None),
    ("fig_dist_log2", None), ("fig_dist_iqr", None),
]:
    if key not in st.session_state:
        st.session_state[key] = default


class _SingleDatasetState(dict):
    """Adapter so the shared render_cleaning_ui/render_qc_ui/render_normalization_ui
    functions can read/write the flat single-dataset session_state keys through the
    same dict-like interface used for multi-dataset mode (where each dataset is its
    own plain dict)."""
    _MAP = {
        "cleaned_data": "cleaned_data", "imputed_data": "imputed_data",
        "missingness_table": "missingness_table", "imputation_method_used": "imputation_method_used",
        "raw_peak_df_qc": "raw_peak_df_qc", "raw_peak_df_qc_full": "raw_peak_df_qc_full",
        "cv_table": "cv_table", "qc_figs": "qc_figs",
        "qc_log_data": "qc_log_data",
        "istd_normalized": "istd_normalized", "istd_mapping": "istd_mapping",
        "iqr_normalized": "iqr_normalized",
        "log2_only": "log2_only", "log2_full": "log2_full",
        "log2_data": "log2_data", "log2_constant": "log2_constant",
        "fig_dist_log2": "fig_dist_log2", "fig_dist_iqr": "fig_dist_iqr",
    }

    def __getitem__(self, key):
        if key == "processing_notes":
            return st.session_state.processing_notes
        return st.session_state[self._MAP[key]]

    def __setitem__(self, key, value):
        if key == "processing_notes":
            st.session_state.processing_notes = value
        else:
            st.session_state[self._MAP[key]] = value

    def get(self, key, default=None):
        try:
            return self[key]
        except KeyError:
            return default


def reset_downstream_analysis_state():
    """
    Clear every downstream-analysis result (Statistics, ANOVA, post-hoc, Biomarkers)
    whenever new data is loaded. Without this, a stale stats_result from a previous
    dataset/mode can persist across a reload, and tabs that read it (Heatmap, Volcano,
    Biomarker Discovery) can then reference a log2 matrix that no longer matches --
    e.g. crashing on mismatched columns -- since those tabs assume stats_result is
    None until freshly (re)computed against the currently-loaded data.
    """
    st.session_state.stats_result = None
    st.session_state.stats_result_groups = None
    st.session_state.stats_result_group_col = None
    st.session_state.anova_result = None
    st.session_state.posthoc_result = None
    st.session_state.anova_groups_used = None
    st.session_state.biomarkers = None


# ---------------------------------------------------------------------------
# FDR column/display label
#
# The multiple-testing correction is always Benjamini-Hochberg
# (statsmodels multipletests(method="fdr_bh"), equivalent to R's
# p.adjust(method="BH")) regardless of assay type -- this never changes.
# Only the displayed/exported column name changes: for Untargeted
# Metabolomics and Untargeted Lipidomics ("unbiased" assays), it's shown
# as "BH P Value" instead of "FDR". Targeted Metabolomics/Lipidomics (and,
# in multi-dataset mode, any combination that mixes in a targeted dataset)
# keep the "FDR" label.
# ---------------------------------------------------------------------------
UNBIASED_DATA_TYPES = {"Untargeted Metabolomics", "Untargeted Lipidomics"}


def current_fdr_label() -> str:
    if st.session_state.data_mode == "multi":
        types = {ds["data_type"] for ds in st.session_state.datasets.values()} if st.session_state.datasets else set()
        if types and types.issubset(UNBIASED_DATA_TYPES):
            return "BH P Value"
        return "FDR"
    return "BH P Value" if st.session_state.data_type in UNBIASED_DATA_TYPES else "FDR"


def stats_fdr_label(stats_df) -> str:
    """
    For tabs that consume an already-computed stats table (Volcano, Heatmap,
    Biomarker Discovery) rather than computing one fresh: read the FDR column
    name directly off that table, so it always matches whatever label it was
    actually built with (current_fdr_label() is for computing a NEW table).
    """
    return "BH P Value" if "BH P Value" in stats_df.columns else "FDR"


def _safe_tag(s):
    """
    Turn any label (a group name, comparison string, dataset name, ...) into a
    filesystem-safe fragment for building descriptive download filenames, e.g.
    "Healthy vs. Diabetic (T2D)" -> "Healthy_vs._Diabetic_T2D". Mirrors
    pathway_analysis.py's `_safe()` so every download across the app names
    files the same way.
    """
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(s)).strip("_") or "Data"


# ---------------------------------------------------------------------------
# Shared "figure grid" layout: lays out a section's figures side-by-side
# instead of stacked vertically, so each section fits on one screen with
# minimal scrolling -- 2 figures -> 2 columns, 3 -> 3 columns, 4 -> a 2x2
# grid, more than 4 -> wraps into additional rows of up to 3. Each figure
# keeps its own heading/caption and its existing download button; only the
# on-screen arrangement changes, nothing about the plot content itself.
# ---------------------------------------------------------------------------
def render_figure_row(items):
    items = [it for it in items if it.get("fig") is not None]
    n = len(items)
    if n == 0:
        return
    if n == 4:
        rows = [items[0:2], items[2:4]]
    elif n <= 3:
        rows = [items]
    else:
        rows = [items[i:i + 3] for i in range(0, n, 3)]

    for row in rows:
        cols = st.columns(len(row))
        for col, item in zip(cols, row):
            with col:
                if item.get("title"):
                    st.markdown(f"**{item['title']}**")
                if item.get("caption"):
                    st.caption(item["caption"])
                st.pyplot(item["fig"], width="content")
                if item.get("download_name"):
                    dataset_manager.render_figure_download(
                        st, item["fig"], item["download_name"],
                        key_prefix=item.get("key_prefix", item["download_name"])
                    )


def render_single_figure(fig, download_name=None, key_prefix=None, title=None, caption=None,
                          width_ratio=(1, 3, 1)):
    """
    Show one standalone figure constrained to a centered fraction of the page
    width (width_ratio, default 60%) instead of letting it sit small and
    left-aligned inside the full 'wide' page container -- the figure's own
    size (in inches) governs how large it renders up to that width, so this
    keeps single-figure sections compact and centered rather than surrounded
    by uneven empty space. Pass a wider ratio (e.g. (1, 6, 1) for 75%) for
    figures that are already sized bigger/denser (heatmap, volcano, boxplot).
    Pass width_ratio=None for a wide, multi-panel strip (e.g. a 4-panel
    distribution row) that should stretch across the full tab width instead
    of being squeezed into a narrow centered column.
    """
    if fig is None:
        return
    if title:
        st.markdown(f"**{title}**")
    if caption:
        st.caption(caption)
    if width_ratio is None:
        st.pyplot(fig, width="stretch")
        if download_name:
            dataset_manager.render_figure_download(st, fig, download_name,
                                                     key_prefix=key_prefix or download_name)
        return
    cols = st.columns(list(width_ratio))
    with cols[1]:
        st.pyplot(fig, width="content")
        if download_name:
            dataset_manager.render_figure_download(st, fig, download_name,
                                                     key_prefix=key_prefix or download_name)


def render_dataset_dropdown(datasets: dict, stage_key: str, stage_col: str) -> str:
    """
    Shared multi-dataset selector for Tabs 2-4: a status dashboard (one row per
    dataset, showing Cleaning/QC/Normalization progress) followed by a dropdown to
    pick which single dataset to configure below — replaces the old nested-tabs UI
    (which meant up to 5 tab levels deep with 4+ datasets) with one flat, always-visible
    overview plus a single selection control. Returns the selected dataset's key.
    `stage_col` (e.g. "2. QC") names which stage this tab is currently working on.
    """
    status_df = dataset_manager.dataset_status_table(datasets)
    st.dataframe(utils.format_df_for_display(status_df), width='stretch', hide_index=True)
    st.caption(f"This page configures the **{stage_col}** column above, one dataset at a time.")

    labels = {ds_id: f"{ds['label']} ({ds['data_type']})" for ds_id, ds in datasets.items()}
    ids = list(labels.keys())
    state_key = f"multi_select_{stage_key}"
    if st.session_state.get(state_key) not in ids:
        st.session_state[state_key] = ids[0]
    selected_id = st.selectbox(
        "Select a dataset/method to configure", options=ids,
        format_func=lambda k: labels[k], key=state_key
    )
    return selected_id


# ===========================================================================
# TAB 1 — DATA UPLOAD
# ===========================================================================
def page_data_upload():
    st.header("Data Upload & Study Configuration")

    st.session_state.data_mode = st.radio(
        "How many datasets are you uploading?",
        ["Single dataset", "Multiple datasets (combine methods/modes)"],
        index=0 if st.session_state.data_mode == "single" else 1,
        horizontal=True,
        help="Real LC-MS studies often run the same samples through multiple chromatography "
             "methods and/or ionization modes (e.g. Method1-Positive, Method1-Negative, "
             "Method2-Positive, Method2-Negative) since no single method captures the full "
             "metabolome. Choose 'Multiple datasets' to clean, QC, and normalize each one "
             "independently, then combine them into one unified matrix for downstream analysis."
    )
    st.session_state.data_mode = "single" if st.session_state.data_mode == "Single dataset" else "multi"

    # =======================================================================
    # SINGLE-DATASET MODE (original workflow, unchanged)
    # =======================================================================
    if st.session_state.data_mode == "single":
        st.session_state.data_type = st.selectbox(
            "Analysis type",
            DATA_TYPES,
        )
        demo_peak_csv, demo_annot_csv, demo_method, demo_desc = STANDARD_DEMOS[st.session_state.data_type]

        col1, col2 = st.columns(2)
        with col1:
            st.subheader("Peak Area Matrix")
            st.caption("First column = Metabolite/Lipid ID, remaining columns = sample peak areas (incl. QC).")
            peak_file = st.file_uploader("Upload peak area matrix (CSV/XLSX)", type=["csv", "xlsx"], key="peak_upload")
            use_demo = st.checkbox("Use built-in demo dataset instead", value=(peak_file is None))
            demo_choice = "Standard"
            if use_demo or peak_file is None:
                demo_choice = st.radio(
                    "Demo dataset", ["Standard", "Rich Clinical Demo"], horizontal=True,
                    key="single_demo_choice"
                )
                if demo_choice == "Standard":
                    st.caption(f"**{st.session_state.data_type}**: {demo_desc}.")
                else:
                    st.caption(
                        "200 metabolites + 1 internal standard (ISTD_D4-Alanine), 36 biological samples "
                        "across 4 diagnosis groups (9 each) + 6 QC replicates = 42 total. Metadata: "
                        "Diagnosis, Age, Gender, Treatment, Ethnicity, Body Weight. Row annotations: Method, "
                        "Pathway, HMDB ID and ISTD (targeted types normalize with the ISTD; untargeted "
                        "types leave the ISTD row out)."
                    )
        with col2:
            st.subheader("Sample Metadata")
            st.caption("Columns: Sample, Group, IsQC (True/False), Batch (optional) — any additional "
                       "columns (Diagnosis, Age, Gender, Treatment, Ethnicity, ...) are kept and become "
                       "selectable grouping/coloring variables in PCA, Statistics, Heatmap, and Boxplot.")
            meta_file = st.file_uploader("Upload metadata table (CSV/XLSX)", type=["csv", "xlsx"], key="meta_upload")

        st.subheader("Metabolite Row Annotations (optional)")
        st.caption("Columns: Metabolite, then any annotation columns (e.g. Method, Pathway) — enables "
                   "the Heatmap page's row annotation tracks. Not required for PCA/Statistics/Volcano/"
                   "Biomarker/Boxplot.")
        row_annot_file = st.file_uploader("Upload row annotation table (CSV/XLSX)",
                                           type=["csv", "xlsx"], key="row_annot_upload")

        if st.button("Load & Validate Data", type="primary", icon="✅"):
            try:
                row_annotations = None
                if use_demo or peak_file is None:
                    if demo_choice == "Rich Clinical Demo":
                        peak_path = os.path.join(os.path.dirname(__file__), "sample_peak_area_matrix_richdemo.csv")
                        meta_path = os.path.join(os.path.dirname(__file__), "sample_metadata_richdemo.csv")
                        annot_path = os.path.join(os.path.dirname(__file__), "metabolite_row_annotations_richdemo.csv")
                        peak_raw = pd.read_csv(peak_path)
                        meta_raw = pd.read_csv(meta_path)
                        row_annotations = pd.read_csv(annot_path).set_index("Metabolite")
                        st.info("Using the built-in **Rich Clinical Demo** dataset (200 metabolites + 1 "
                                "internal standard, 42 samples, 6 metadata variables, row annotations "
                                "with Method, Pathway, HMDB ID and ISTD).")
                    else:
                        here = os.path.dirname(__file__)
                        peak_raw = pd.read_csv(os.path.join(here, demo_peak_csv))
                        meta_raw = pd.read_csv(os.path.join(here, "sample_metadata_targeted.csv"))
                        annot_path = os.path.join(here, demo_annot_csv)
                        if os.path.exists(annot_path):
                            ann = pd.read_csv(annot_path)
                            if demo_method:
                                ann = ann[ann["Method"] == demo_method]
                            row_annotations = ann.set_index("Metabolite")
                        st.info(f"Using the built-in **Standard** demo dataset for "
                                f"'{st.session_state.data_type}': {demo_desc}.")
                else:
                    peak_raw = utils.load_table(peak_file)
                    meta_raw = utils.load_table(meta_file) if meta_file is not None else None
                    if meta_raw is None:
                        st.error("Please upload a metadata table (or check 'use demo dataset').")
                        st.stop()

                peak_df = utils.validate_peak_matrix(peak_raw)
                meta = utils.validate_metadata(meta_raw, peak_df.columns)
                qc_cols, sample_cols = utils.split_qc_and_samples(peak_df, meta)

                if row_annot_file is not None and row_annotations is None:
                    row_annot_raw = utils.load_table(row_annot_file)
                    if "Metabolite" not in row_annot_raw.columns:
                        st.warning("Row annotation file has no 'Metabolite' column — ignoring it.")
                    else:
                        row_annotations = row_annot_raw.set_index("Metabolite")

                st.session_state.raw_peak_df = peak_df
                st.session_state.meta = meta
                st.session_state.qc_cols = qc_cols
                st.session_state.sample_cols = sample_cols
                st.session_state.row_annotations = row_annotations
                st.session_state.processing_notes = [f"Analysis type: {st.session_state.data_type}"]
                st.session_state.log2_data = None
                st.session_state.log2_only = None
                st.session_state.log2_full = None
                st.session_state.fig_dist_log2 = None
                st.session_state.fig_dist_iqr = None
                st.session_state.raw_peak_df_qc = None
                st.session_state.raw_peak_df_qc_full = None
                # results of the previous dataset's cleaning/QC/normalization must not carry over
                for _k in ("cleaned_data", "imputed_data", "missingness_table", "imputation_method_used",
                           "cv_table", "qc_log_data", "istd_normalized", "istd_mapping", "iqr_normalized"):
                    st.session_state[_k] = None
                reset_downstream_analysis_state()

                st.success(f"Loaded {peak_df.shape[0]} features × {peak_df.shape[1]} samples "
                           f"({len(qc_cols)} QC, {len(sample_cols)} study samples)."
                           + (f" Row annotations: {', '.join(row_annotations.columns)}."
                              if row_annotations is not None else ""))
            except utils.DataValidationError as e:
                st.error(f"Validation error: {e}")
            except Exception as e:
                st.error(f"Unexpected error loading data: {e}")

        if st.session_state.raw_peak_df is not None:
            st.subheader("Preview")
            st.dataframe(utils.format_df_for_display(st.session_state.raw_peak_df.head(10)), width='stretch')
            st.write("**Metadata:**")
            st.dataframe(utils.format_df_for_display(st.session_state.meta), width='stretch')
            if st.session_state.get("row_annotations") is not None:
                st.write("**Row Annotations:**")
                st.dataframe(utils.format_df_for_display(st.session_state.row_annotations.head(10)), width='stretch')

    # =======================================================================
    # MULTI-DATASET MODE (multiple methods/modes, combined downstream)
    # =======================================================================
    else:
        st.caption(
            "Each dataset is cleaned, QC'd, and normalized **independently**. Once all are "
            "normalized, click **🔗 Generate Combined Normalized Data** in the Normalization "
            "page — the combined table then feeds every downstream page."
        )

        data_source = st.radio(
            "Data source", ["Demo Data", "Real Data"], horizontal=True, key="multi_data_source"
        )

        st.subheader("Metabolite Row Annotations (optional)")
        st.caption("Columns: Metabolite, then any annotation columns (e.g. Method, Pathway) — enables "
                   "the Heatmap page's row annotation tracks. Matches either the full prefixed name "
                   "(e.g. `Targeted Metabolomics::Glucose`) or the plain name.")
        multi_row_annot_file = st.file_uploader("Upload row annotation table (CSV/XLSX)",
                                                 type=["csv", "xlsx"], key="multi_row_annot_upload")

        # ===================================================================
        # MULTIPLE DATASETS — DEMO DATA
        # ===================================================================
        if data_source == "Demo Data":
            st.caption(
                (f"Loads {len(DATA_TYPES)} demo datasets ({', '.join(DATA_TYPES)}) sharing 24 biological "
                 f"samples — real {'lipid' if APP_FAMILY == 'Lipidomics' else 'metabolite'} names "
                 "(row annotations included).")
                if APP_FAMILY else
                "Loads 4 demo datasets (Untargeted/Targeted Metabolomics, Untargeted/Targeted "
                "Lipidomics) sharing 24 biological samples — real metabolite and lipid names with "
                "HMDB IDs (row annotations included)."
            )
            if st.button("Load Demo Datasets", type="primary", icon="🧪"):
                meta_path = os.path.join(os.path.dirname(__file__), "sample_metadata_multimethod.csv")
                meta_raw = pd.read_csv(meta_path)
                demo_types = list(DATA_TYPES)
                datasets = {}
                first_peak_df = None
                for dtype in demo_types:
                    prefix = dtype.lower().replace(" ", "_")
                    peak_path = os.path.join(
                        os.path.dirname(__file__),
                        f"sample_peak_area_matrix_multimethod_{prefix}.csv"
                    )
                    peak_raw = pd.read_csv(peak_path)
                    peak_df = utils.validate_peak_matrix(peak_raw)
                    if first_peak_df is None:
                        first_peak_df = peak_df
                    entry = dataset_manager.make_dataset_entry(dtype, dtype)
                    entry["raw_df"] = peak_df
                    datasets[dtype] = entry
                meta = utils.validate_metadata(meta_raw, first_peak_df.columns)
                qc_cols, sample_cols = utils.split_qc_and_samples(first_peak_df, meta)
                for entry in datasets.values():
                    entry["qc_cols"] = qc_cols
                    entry["sample_cols"] = sample_cols
                st.session_state.datasets = datasets
                st.session_state.meta = meta
                st.session_state.qc_cols = qc_cols
                st.session_state.sample_cols = sample_cols
                st.session_state.combined_done = False
                st.session_state.log2_data = None
                st.session_state.log2_only = None
                st.session_state.fig_dist_log2 = None
                st.session_state.fig_dist_iqr = None
                row_annotations = None
                default_annot_path = os.path.join(os.path.dirname(__file__),
                                                    "metabolite_row_annotations_multimethod.csv")
                if os.path.exists(default_annot_path):
                    row_annotations = pd.read_csv(default_annot_path).set_index("Metabolite")
                if multi_row_annot_file is not None:
                    row_annot_raw = utils.load_table(multi_row_annot_file)
                    if "Metabolite" not in row_annot_raw.columns:
                        st.warning("Row annotation file has no 'Metabolite' column — ignoring it.")
                    else:
                        row_annotations = row_annot_raw.set_index("Metabolite")
                st.session_state.row_annotations = row_annotations
                reset_downstream_analysis_state()
                st.session_state.processing_notes = [
                    f"Multi-dataset mode (demo): {len(datasets)} datasets loaded "
                    f"({', '.join(datasets.keys())}), sharing {len(sample_cols)} biological samples."
                ]
                st.success(f"Loaded {len(datasets)} demo datasets sharing {len(sample_cols)} biological "
                           f"samples across {meta['Group'].nunique() - 1} groups."
                           + (f" Row annotations: {', '.join(row_annotations.columns)}."
                              if row_annotations is not None else ""))

        # ===================================================================
        # MULTIPLE DATASETS — REAL DATA
        # ===================================================================
        else:
            st.subheader("Upload Your Datasets")
            n_datasets = st.number_input("Number of datasets", min_value=1, max_value=20, value=4, step=1)
            st.caption("Typical: 2 methods × 2 ionization modes = 4. Add more for additional methods "
                       "— e.g. targeted lipidomics run as several class-specific batches (up to 20).")

            upload_entries = []
            for i in range(int(n_datasets)):
                with st.expander(f"Dataset {i+1}", expanded=(i < 2)):
                    c1, c2 = st.columns(2)
                    label = c1.text_input(f"Dataset {i+1} label", value=f"Method{i+1}", key=f"ds_label_{i}")
                    dtype = c2.selectbox(f"Analysis type", DATA_TYPES,
                                          key=f"ds_type_{i}")
                    ds_file = st.file_uploader(f"Peak area matrix for {label}", type=["csv", "xlsx"], key=f"ds_file_{i}")
                    upload_entries.append((label, dtype, ds_file))

            st.subheader("Sample Metadata")
            st.caption(
                "One metadata table applies to ALL datasets above — the same biological samples "
                "(by name) must appear in every peak area matrix you uploaded. Columns: Sample, "
                "Group, IsQC (True/False), Batch (optional)."
            )
            shared_meta_file = st.file_uploader("Upload shared metadata table (CSV/XLSX)",
                                                 type=["csv", "xlsx"], key="shared_meta_upload")

            if st.button("Load All Datasets", type="primary", icon="📥"):
                if shared_meta_file is None:
                    st.error("Please upload the shared metadata table above — it's required to "
                             "identify QC vs. biological samples and group assignments for every "
                             "dataset.")
                    st.stop()
                try:
                    meta_raw = utils.load_table(shared_meta_file)
                    datasets = {}
                    first_peak_df = None
                    for label, dtype, ds_file in upload_entries:
                        if ds_file is None:
                            st.warning(f"Skipping '{label}' — no file uploaded.")
                            continue
                        peak_raw = utils.load_table(ds_file)
                        peak_df = utils.validate_peak_matrix(peak_raw)
                        if first_peak_df is None:
                            first_peak_df = peak_df
                        entry = dataset_manager.make_dataset_entry(label, dtype)
                        entry["raw_df"] = peak_df
                        datasets[label] = entry

                    if not datasets:
                        st.error("No datasets were uploaded.")
                        st.stop()

                    meta = utils.validate_metadata(meta_raw, first_peak_df.columns)
                    warnings = dataset_manager.validate_shared_samples(datasets, meta)
                    for w in warnings:
                        st.warning(w)

                    qc_cols, sample_cols = utils.split_qc_and_samples(first_peak_df, meta)
                    for entry in datasets.values():
                        entry["qc_cols"] = [c for c in qc_cols if c in entry["raw_df"].columns]
                        entry["sample_cols"] = [c for c in sample_cols if c in entry["raw_df"].columns]

                    st.session_state.datasets = datasets
                    st.session_state.meta = meta
                    st.session_state.qc_cols = qc_cols
                    st.session_state.sample_cols = sample_cols
                    st.session_state.combined_done = False
                    st.session_state.log2_data = None
                    st.session_state.log2_only = None
                    st.session_state.fig_dist_log2 = None
                    st.session_state.fig_dist_iqr = None
                    row_annotations = None
                    if multi_row_annot_file is not None:
                        row_annot_raw = utils.load_table(multi_row_annot_file)
                        if "Metabolite" not in row_annot_raw.columns:
                            st.warning("Row annotation file has no 'Metabolite' column — ignoring it.")
                        else:
                            row_annotations = row_annot_raw.set_index("Metabolite")
                    st.session_state.row_annotations = row_annotations
                    reset_downstream_analysis_state()
                    st.session_state.processing_notes = [
                        f"Multi-dataset mode (real data): {len(datasets)} datasets loaded "
                        f"({', '.join(datasets.keys())}), sharing {len(sample_cols)} biological samples."
                    ]
                    st.success(f"Loaded {len(datasets)} datasets."
                               + (f" Row annotations: {', '.join(row_annotations.columns)}."
                                  if row_annotations is not None else ""))
                except utils.DataValidationError as e:
                    st.error(f"Validation error: {e}")
                except Exception as e:
                    st.error(f"Unexpected error loading data: {e}")

        if st.session_state.datasets:
            st.subheader("Loaded Datasets")
            summary_rows = [{"Label": ds["label"], "Analysis Type": ds["data_type"],
                             "Features": ds["raw_df"].shape[0] if ds["raw_df"] is not None else 0,
                             "Samples": ds["raw_df"].shape[1] if ds["raw_df"] is not None else 0}
                            for ds in st.session_state.datasets.values()]
            st.dataframe(utils.format_df_for_display(pd.DataFrame(summary_rows)), width='stretch')
            st.write("**Shared Metadata:**")
            st.dataframe(utils.format_df_for_display(st.session_state.meta), width='stretch')

            st.divider()
            st.subheader("🚀 Quick Start")
            already_processed = dataset_manager.all_datasets_ready(st.session_state.datasets)
            if already_processed:
                st.success(
                    "All datasets are cleaned, QC'd, and normalized — go to **Combined Normalized "
                    "Data → 🔗 Generate Combined Normalized Data** to combine them. That combined "
                    "table is required before any downstream page (PCA, Statistics, Volcano, "
                    "Biomarker Discovery, Heatmap, Boxplot) can be used."
                )
            st.caption(
                "Applies recommended defaults to every dataset: missingness filter (30%) → "
                "Half-Minimum imputation → QC CV filter (CV > 20% removed for targeted datasets only; "
                "untargeted datasets keep every feature, with CV still reported) → normalization → Log2. "
                "Each targeted dataset is normalized with its internal standard (the ISTD named in the row "
                "annotations' ISTD column, or the ISTD_… row). "
                "Does not combine them — that's a separate step (Combined Normalization Data "
                "(Multi-Dataset) page)."
            )
            if st.button("🚀 Auto-Process All Datasets (recommended defaults)", type="primary"):
                with st.spinner("Cleaning, QC'ing, and normalizing every dataset independently..."):
                    try:
                        notes = dataset_manager.auto_process_all_datasets(
                            st.session_state.datasets, row_annotations=st.session_state.get("row_annotations"))
                    except ValueError as e:
                        st.error(str(e))
                        st.stop()
                st.session_state.processing_notes = (
                    [f"Multi-dataset mode: {len(st.session_state.datasets)} datasets loaded "
                     f"({', '.join(ds['label'] for ds in st.session_state.datasets.values())}), "
                     f"sharing {len(st.session_state.sample_cols)} biological samples."]
                    + notes
                )
                st.success(
                    f"All {len(st.session_state.datasets)} datasets cleaned, QC'd, and normalized "
                    "independently. Go to **Combined Normalized Data → 🔗 Generate Combined "
                    "Normalized Data** to combine them — required before PCA, Statistics, or any "
                    "other downstream page."
                )
                st.dataframe(utils.format_df_for_display(dataset_manager.dataset_status_table(st.session_state.datasets)),
                             width='stretch', hide_index=True)


def _state_istds(features, state):
    """Internal-standard features of a dataset (named 'ISTD...' or listed in the row annotations'
    ISTD column). Feature filters never remove these: they are needed for ISTD normalization."""
    label = state.get("label")
    data_type = state.get("data_type") or st.session_state.get("data_type")
    return normalization.detect_istds(features, st.session_state.get("row_annotations"), label, data_type)


# ---------------------------------------------------------------------------
# Reusable QC Validation UI, shared between single-dataset mode and each
# dataset in multi-dataset mode.
# ---------------------------------------------------------------------------
def render_qc_ui(peak_df, qc_cols, sample_cols, key_prefix, state):
    if len(qc_cols) < 2:
        st.warning("At least 2 QC samples are required for CV calculation. Skipping QC module.")
        return

    cv_table = qc.calculate_cv(peak_df[qc_cols])
    state["cv_table"] = cv_table
    n_accept = (cv_table["Quality"] == "Acceptable").sum()
    n_var = (cv_table["Quality"] == "Variable").sum()

    m1, m2, m3 = st.columns(3)
    m1.metric("Total Features", len(cv_table))
    m2.metric("Acceptable (CV≤20%)", n_accept)
    m3.metric("Variable (CV>20%)", n_var)

    st.markdown("**QC Metrics & Visualization**")
    st.caption("The sample correlation panel uses **QC replicates only**.")
    fig_cv_dist = qc.cv_distribution_plot(cv_table)
    fig_cv_hist = qc.cv_histogram(cv_table)
    fig_peak_hist = qc.peak_area_histogram(peak_df, title=f"{key_prefix} — Peak Area Distribution")
    qc_log = stats_analysis.strict_log2(peak_df[qc_cols])
    state["qc_log_data"] = qc_log
    fig_corr, corr_df = qc.sample_correlation_matrix(qc_log)
    render_figure_row([
        {"title": "QC CV Distribution", "fig": fig_cv_dist,
         "download_name": f"{key_prefix}_QC_CV_Distribution", "key_prefix": f"{key_prefix}_cvdist"},
        {"title": "Feature Counts by CV Quality", "fig": fig_cv_hist,
         "download_name": f"{key_prefix}_Feature_Counts_by_Quality", "key_prefix": f"{key_prefix}_cvhist"},
        {"title": "Peak Area Distribution", "fig": fig_peak_hist,
         "caption": "Raw peak-area values across all features and samples (QC + biological), log10 scale.",
         "download_name": f"{key_prefix}_Peak_Area_Histogram", "key_prefix": f"{key_prefix}_peakhist"},
        {"title": "QC Sample Correlation", "fig": fig_corr,
         "download_name": f"{key_prefix}_QC_Sample_Correlation", "key_prefix": f"{key_prefix}_corr"},
    ])

    st.subheader("CV Filtering Table")
    st.dataframe(utils.format_df_for_display(cv_table), width='stretch', height=250)
    st.download_button(
        f"Download {key_prefix.replace(' ', '_')}_CV_Table.csv",
        utils.to_download_bytes_csv(cv_table),
        file_name=f"{key_prefix.replace(' ', '_')}_CV_Table.csv", mime="text/csv",
        key=f"{key_prefix}_dl_cv_table"
    )

    filter_cv = st.checkbox(
        "Filter out 'Variable' features (CV>20%) from downstream analysis",
        value=False, key=f"{key_prefix}_filter_cv"
    )

    if st.button("Confirm QC & Proceed", type="primary", icon="✅", key=f"{key_prefix}_confirm_qc"):
        if state.get("imputed_data") is not None:
            working_df = state["imputed_data"].copy()
            source_desc = "cleaned & imputed data from Data Cleaning & Imputation"
        else:
            working_df = peak_df[sample_cols].copy()
            source_desc = "raw biological samples (no cleaning/imputation was applied)"
        # Keep the full (pre-CV-filter) feature table too: even when the CV filter
        # below removes "Variable" features from the table used downstream, the
        # per-sample Median-IQR statistics in Normalization are computed from this
        # complete detected-feature background, not just the reduced/exported
        # subset -- so enabling the CV filter never changes the per-sample scaling
        # statistics, only which rows are kept in the final normalized table.
        state["raw_peak_df_qc_full"] = working_df
        if filter_cv:
            keep = cv_table.index[cv_table["Quality"] == "Acceptable"]
            working_df = normalization.keep_istds(working_df.loc[working_df.index.intersection(keep)],
                                                  state["raw_peak_df_qc_full"],
                                                  _state_istds(working_df.index, state))
        state["raw_peak_df_qc"] = working_df
        state["qc_figs"] = {"QC CV Distribution": fig_cv_dist, "QC Sample Correlation Matrix": fig_corr}
        state["processing_notes"].append(
            f"QC validation: {len(cv_table)} features assessed via {len(qc_cols)} QC replicates; "
            f"{'variable features (CV>20%) removed' if filter_cv else 'no CV-based filtering applied — all features kept'}, "
            f"applied on top of {source_desc}. "
            f"{len(qc_cols)} QC samples excluded from all downstream normalization/statistics/visualization "
            f"({working_df.shape[0]} features, {len(sample_cols)} biological samples proceed)."
            + (f" Median-IQR sample statistics (Normalization) will still be computed from all "
               f"{state['raw_peak_df_qc_full'].shape[0]} pre-filter features, per sample."
               if filter_cv else "")
        )
        st.success(
            f"QC step confirmed — {len(qc_cols)} QC samples excluded going forward. "
            f"Proceed to Normalization with {working_df.shape[0]} features across "
            f"{len(sample_cols)} biological samples."
        )


# ===========================================================================
# TAB 3 — QC VALIDATION
# ===========================================================================
def page_qc():
    st.header("Quality Control (QC) Validation")
    if st.session_state.data_mode == "single":
        if st.session_state.raw_peak_df is None:
            st.warning("Load data on the Data Upload page first.")
        else:
            render_qc_ui(st.session_state.raw_peak_df, st.session_state.qc_cols,
                         st.session_state.sample_cols, "single", _SingleDatasetState())
    else:
        if not st.session_state.datasets:
            st.warning("Load datasets on the Data Upload page first.")
        else:
            selected_id = render_dataset_dropdown(st.session_state.datasets, "qc", "2. QC")
            ds = st.session_state.datasets[selected_id]
            st.divider()
            st.subheader(f"QC for {ds['label']} ({ds['data_type']})")
            render_qc_ui(ds["raw_df"], ds["qc_cols"], ds["sample_cols"], selected_id, ds)
            st.caption("See **Combined QC (Multi-Dataset)** for QC diagnostics pooled across all datasets.")

# ===========================================================================
# TAB 4 — COMBINED QC
# ===========================================================================
def page_combined_qc():
    st.header("Combined QC (Multi-Dataset) — Across All Datasets")
    if st.session_state.data_mode == "single":
        st.info("Combined QC only applies in multi-dataset mode (Data Upload). In single-dataset "
                "mode, see QC Validation for QC diagnostics.")
    elif not st.session_state.datasets:
        st.warning("Load datasets on the Data Upload page first.")
    else:
        st.subheader("QC Comparison Across All Datasets")
        comparison_table = dataset_manager.qc_comparison_table(st.session_state.datasets)
        if comparison_table.empty:
            st.info("Run QC on at least one dataset in QC Validation to populate this comparison.")
        else:
            st.dataframe(utils.format_df_for_display(comparison_table), width='stretch')
            fig_compare = dataset_manager.qc_comparison_plot(st.session_state.datasets)
            if fig_compare is not None:
                render_single_figure(fig_compare, download_name="QC_Comparison_Across_Datasets",
                                      key_prefix="qc_compare")

        st.divider()
        st.subheader("Combined QC Analysis")
        st.caption("Per-dataset CV tables pooled into one, with feature names prefixed by dataset.")
        combined_cv = dataset_manager.combine_cv_tables(st.session_state.datasets)
        if combined_cv is None:
            st.info("Run QC on at least one dataset in QC Validation to populate the combined analysis.")
        else:
            cc1, cc2, cc3 = st.columns(3)
            cc1.metric("Total Features (all datasets)", len(combined_cv))
            cc2.metric("Acceptable (CV≤20%)", int((combined_cv["Quality"] == "Acceptable").sum()))
            cc3.metric("Variable (CV>20%)", int((combined_cv["Quality"] == "Variable").sum()))

            st.divider()
            st.subheader("Apply CV Filter to All Datasets")
            st.caption(
                "Sets the CV filter for every QC-confirmed dataset at once. Re-run Normalization "
                "(Normalization) and regenerate the combined data (Combined Normalization Data (Multi-Dataset)) afterward if this changes a "
                "dataset's feature count."
            )
            combined_filter_cv = st.checkbox(
                "Filter out 'Variable' features (CV>20%) from downstream analysis",
                value=False, key="combined_filter_cv"
            )
            eligible_datasets = {
                label: ds for label, ds in st.session_state.datasets.items()
                if ds.get("raw_peak_df_qc_full") is not None and ds.get("cv_table") is not None
            }
            if not eligible_datasets:
                st.info("Confirm QC (QC Validation → 'Confirm QC & Proceed') for at least one dataset first.")
            else:
                st.caption(f"{len(eligible_datasets)}/{len(st.session_state.datasets)} dataset(s) "
                           "have completed QC confirmation and are eligible below.")
                if st.button("Apply to All Eligible Datasets", type="primary", icon="⚙️",
                              key="combined_apply_cv_filter"):
                    summary_rows = []
                    for label, ds in eligible_datasets.items():
                        full_df = ds["raw_peak_df_qc_full"]
                        cv_table_ds = ds["cv_table"]
                        if combined_filter_cv:
                            keep = cv_table_ds.index[cv_table_ds["Quality"] == "Acceptable"]
                            new_working = normalization.keep_istds(
                                full_df.loc[full_df.index.intersection(keep)], full_df,
                                _state_istds(full_df.index, ds))
                        else:
                            new_working = full_df.copy()
                        before_n = (
                            ds["raw_peak_df_qc"].shape[0]
                            if ds.get("raw_peak_df_qc") is not None else full_df.shape[0]
                        )
                        ds["raw_peak_df_qc"] = new_working
                        # Invalidate everything computed from the old feature set -- it must be
                        # re-derived (Normalization, or Quick Start) before Statistics.
                        ds["log2_only"] = None
                        ds["log2_full"] = None
                        ds["log2_data"] = None
                        ds["fig_dist_log2"] = None
                        ds["fig_dist_iqr"] = None
                        ds["istd_normalized"] = None
                        ds["istd_mapping"] = None
                        ds["iqr_normalized"] = None
                        ds["processing_notes"].append(
                            f"{ds['label']}: CV filter re-applied from Combined QC "
                            + ("(variable features removed); " if combined_filter_cv
                               else "(no CV-based filtering — all features kept); ")
                            + f"{before_n} → {new_working.shape[0]} features. Re-run Normalization "
                            "(Normalization) or Quick Start before Statistics."
                        )
                        summary_rows.append({
                            "Dataset": ds["label"], "Features Before": before_n,
                            "Features After": new_working.shape[0]
                        })
                    st.session_state.combined_done = False
                    st.session_state.log2_data = None
                    st.session_state.log2_only = None
                    st.session_state.fig_dist_log2 = None
                    st.session_state.fig_dist_iqr = None
                    reset_downstream_analysis_state()
                    st.success(
                        f"CV filter applied to {len(eligible_datasets)} dataset(s). Re-run "
                        "Normalization for each affected dataset, or Quick Start again, "
                        "then regenerate the combined normalized data (Combined Normalization Data (Multi-Dataset)) before Statistics."
                    )
                    st.dataframe(utils.format_df_for_display(pd.DataFrame(summary_rows)), width='stretch', hide_index=True)

            fig_cv_dist_c = qc.cv_distribution_plot(combined_cv)
            fig_cv_hist_c = qc.cv_histogram(combined_cv)
            combined_qc_log = dataset_manager.combine_qc_log_data(
                st.session_state.datasets, st.session_state.qc_cols
            )
            fig_corr_c, corr_df_c = (
                qc.sample_correlation_matrix(combined_qc_log)
                if combined_qc_log is not None else (None, None)
            )
            combined_imputed, _ = dataset_manager.combine_imputed_datasets(
                st.session_state.datasets, st.session_state.sample_cols
            )
            fig_peak_hist_c = (
                qc.peak_area_histogram(combined_imputed, title="Combined — Peak Area Distribution (all datasets)")
                if combined_imputed is not None else None
            )

            qc_row1 = st.columns(2)
            with qc_row1[0]:
                st.markdown("**QC CV Distribution**")
                st.pyplot(fig_cv_dist_c, width="content")
                st.download_button(
                    "Download Combined_QC_CV_Table.csv", utils.to_download_bytes_csv(combined_cv),
                    file_name="Combined_QC_CV_Table.csv", mime="text/csv", key="dl_combined_cv_dist_data"
                )
                dataset_manager.render_figure_download(st, fig_cv_dist_c, "Combined_QC_CV_Distribution",
                                                        key_prefix="combined_cvdist")
            with qc_row1[1]:
                st.markdown("**Feature Counts by CV Quality**")
                st.pyplot(fig_cv_hist_c, width="content")
                quality_counts = combined_cv["Quality"].value_counts().rename_axis("Quality").reset_index(name="Count")
                st.download_button(
                    "Download Combined_Feature_Counts_by_Quality.csv",
                    utils.to_download_bytes_csv(quality_counts),
                    file_name="Combined_Feature_Counts_by_Quality.csv", mime="text/csv",
                    key="dl_combined_cv_hist_data"
                )
                dataset_manager.render_figure_download(st, fig_cv_hist_c, "Combined_Feature_Counts_by_Quality",
                                                        key_prefix="combined_cvhist")

            qc_row2 = st.columns(2)
            with qc_row2[0]:
                st.markdown("**QC Sample Correlation Matrix**")
                if fig_corr_c is None:
                    st.info("Run QC on at least one dataset in QC Validation to populate the correlation matrix.")
                else:
                    st.pyplot(fig_corr_c, width="content")
                    st.download_button(
                        "Download Combined_QC_Sample_Correlation.csv",
                        utils.to_download_bytes_csv(corr_df_c),
                        file_name="Combined_QC_Sample_Correlation.csv", mime="text/csv",
                        key="dl_combined_corr_data"
                    )
                    dataset_manager.render_figure_download(st, fig_corr_c, "Combined_QC_Sample_Correlation",
                                                            key_prefix="combined_corr")
            with qc_row2[1]:
                st.markdown("**Combined Peak Area Distribution**")
                st.caption("Raw peak areas across all features, samples, and datasets pooled together.")
                if fig_peak_hist_c is None:
                    st.info("Run Cleaning & Imputation on at least one dataset first.")
                else:
                    st.pyplot(fig_peak_hist_c, width="content")
                    dataset_manager.render_figure_download(st, fig_peak_hist_c, "Combined_Peak_Area_Histogram",
                                                            key_prefix="combined_peakhist")

# ===========================================================================
# ---------------------------------------------------------------------------
# Reusable Data Cleaning & Imputation UI, shared between single-dataset mode
# and each dataset in multi-dataset mode.
# ---------------------------------------------------------------------------
def render_cleaning_ui(source_df, key_prefix, state):
    treat_zero_as_missing = st.checkbox(
        "Treat exact-zero values as missing (common LC-MS convention for 'not detected')",
        value=True, key=f"{key_prefix}_treat_zero"
    )

    st.subheader("Step 1: Assess & Filter by Missingness")
    missingness_table = imputation_module.compute_missingness(source_df, treat_zero_as_missing)
    state["missingness_table"] = missingness_table
    counts = missingness_table["Quality"].value_counts()

    # Compact side-by-side layout: parameters/metrics on the left, a small
    # thumbnail-sized distribution figure (with its own download) on the right
    # -- keeps this step to roughly one screen instead of a tall, figure-heavy
    # stack of full-width elements.
    param_col, fig_col = st.columns([3, 2])

    with param_col:
        st.caption(
            "Guidance: **<20%** missing → keep as-is · **20-50%** → keep, impute carefully · "
            "**50-70%** → usually remove · **>70-80%** → remove unless biologically essential."
        )
        mm1, mm2 = st.columns(2)
        mm1.metric("Keep (<20%)", int(counts.get("Keep", 0)))
        mm2.metric("Careful imputation (20-50%)", int(counts.get("Keep (careful imputation)", 0)))
        mm3, mm4 = st.columns(2)
        mm3.metric("Usually remove (50-70%)", int(counts.get("Usually remove", 0)))
        mm4.metric("Remove unless essential (>70%)", int(counts.get("Remove unless essential", 0)))

        max_pct_missing = st.slider(
            "Remove features with missingness above this threshold (%)",
            min_value=0, max_value=100, value=30, step=5, key=f"{key_prefix}_max_pct"
        )
        cleaned_preview = imputation_module.filter_by_missingness(source_df, missingness_table, max_pct_missing)
        cleaned_preview = normalization.keep_istds(cleaned_preview, source_df,
                                                   _state_istds(source_df.index, state))
        st.caption(f"At this threshold: **{cleaned_preview.shape[0]} of {source_df.shape[0]}** "
                   f"features retained ({source_df.shape[0] - cleaned_preview.shape[0]} removed).")

    with fig_col:
        fig_miss, ax_miss = _plt.subplots(figsize=(3.2, 2.0))
        ax_miss.hist(missingness_table["Pct_Missing"], bins=30, color="#4C72B0", edgecolor="white")
        for x, lbl in [(20, "20%"), (50, "50%"), (70, "70%")]:
            ax_miss.axvline(x, color="red", linestyle="--", lw=0.8)
        ax_miss.set_xlabel("% Missing", fontsize=7)
        ax_miss.set_ylabel("Number of features", fontsize=7)
        ax_miss.tick_params(labelsize=6)
        ax_miss.set_title("Missingness Distribution", fontsize=8, fontweight="bold", loc="center")
        fig_miss.tight_layout()
        st.pyplot(fig_miss, width="content")
        dataset_manager.render_figure_download(
            st, fig_miss, f"{key_prefix.replace(' ', '_')}_Missingness_Distribution",
            key_prefix=f"{key_prefix}_missdist"
        )

    with st.expander("View missingness table"):
        st.dataframe(utils.format_df_for_display(missingness_table), width='stretch', height=250)

    if st.button("Apply Missingness Filter", type="primary", icon="🧹", key=f"{key_prefix}_apply_filter"):
        state["cleaned_data"] = cleaned_preview
        state["processing_notes"].append(
            f"Missingness filtering: removed {source_df.shape[0] - cleaned_preview.shape[0]} of "
            f"{source_df.shape[0]} features exceeding {max_pct_missing}% missingness; "
            f"{cleaned_preview.shape[0]} features retained."
        )
        st.success(f"Filter applied — {cleaned_preview.shape[0]} features retained.")

    if state.get("cleaned_data") is not None:
        remaining_missing = imputation_module.count_missing(state["cleaned_data"], treat_zero_as_missing)

        st.subheader("Step 2: Choose an Imputation Method")
        if remaining_missing == 0:
            st.success("No missing values remain in the filtered data — imputation isn't needed. "
                       "You can skip ahead.")
            state["imputed_data"] = state["cleaned_data"]
        else:
            st.caption(f"{remaining_missing} missing values remain across the filtered features.")
            _imputation_methods = list(imputation_module.METHOD_INFO.keys())
            _default_method_idx = (
                _imputation_methods.index("K-Nearest Neighbors (KNN)")
                if "K-Nearest Neighbors (KNN)" in _imputation_methods else 0
            )
            method = st.selectbox("Imputation method", _imputation_methods,
                                   index=_default_method_idx, key=f"{key_prefix}_method")
            st.caption(imputation_module.METHOD_INFO[method])

            kwargs = {}
            if method == "K-Nearest Neighbors (KNN)":
                kwargs["n_neighbors"] = st.slider("Number of neighbors (k)", 2, 15, 5, key=f"{key_prefix}_knn_k")
            elif method == "Random Forest (MissForest-style)":
                c1, c2 = st.columns(2)
                kwargs["n_estimators"] = c1.slider("Number of trees", 5, 50, 10, key=f"{key_prefix}_rf_n")
                kwargs["max_iter"] = c2.slider("Iterations", 1, 10, 3, key=f"{key_prefix}_rf_iter")
            elif method == "BPCA (approximated)":
                kwargs["max_iter"] = st.slider("Iterations", 1, 15, 5, key=f"{key_prefix}_bpca_iter")

            if st.button("Apply Imputation", type="primary", icon="🧩", key=f"{key_prefix}_apply_impute"):
                with st.spinner(f"Running {method}... this may take a moment for Random Forest or BPCA."):
                    imputed = imputation_module.impute(
                        state["cleaned_data"], method, treat_zero_as_missing=treat_zero_as_missing, **kwargs
                    )
                state["imputed_data"] = imputed
                state["imputation_method_used"] = method
                state["processing_notes"].append(
                    f"Missing value imputation: {method} applied to "
                    f"{state['cleaned_data'].shape[0]} features ({remaining_missing} missing values filled)."
                )
                st.success(f"Imputation complete using {method}.")

        if state.get("imputed_data") is not None:
            st.subheader("Result: Complete Post-Imputation Data")
            st.caption(f"{state['imputed_data'].shape[0]} features × {state['imputed_data'].shape[1]} samples.")
            st.dataframe(utils.format_df_for_display(state["imputed_data"]), width='stretch', height=400)
            imputed_filename = (f"{key_prefix.replace(' ', '_')}_Imputed_Data.csv"
                                 if key_prefix != "single" else "Imputed_Data.csv")
            st.download_button(
                f"Download {imputed_filename}",
                utils.to_download_bytes_csv(state["imputed_data"]),
                file_name=imputed_filename, mime="text/csv", key=f"{key_prefix}_dl_imputed"
            )


# ===========================================================================
# TAB 2 — DATA CLEANING & MISSING VALUE IMPUTATION
# ===========================================================================
def page_cleaning():
    st.header("Data Cleaning & Missing Value Imputation")
    st.caption("Optional. Runs on biological samples only (QC replicates excluded).")
    if st.session_state.data_mode == "single":
        if st.session_state.raw_peak_df is None:
            st.warning("Load data on the Data Upload page first.")
        else:
            source_df = st.session_state.raw_peak_df[st.session_state.sample_cols]
            render_cleaning_ui(source_df, "single", _SingleDatasetState())
    else:
        if not st.session_state.datasets:
            st.warning("Load datasets on the Data Upload page first.")
        else:
            selected_id = render_dataset_dropdown(st.session_state.datasets, "clean", "1. Cleaning")
            ds = st.session_state.datasets[selected_id]
            st.divider()
            st.subheader(f"Configuring: {ds['label']} ({ds['data_type']})")
            source_df = ds["raw_df"][ds["sample_cols"]]
            render_cleaning_ui(source_df, selected_id, ds)

# ===========================================================================
# ---------------------------------------------------------------------------
# Reusable Normalization UI, shared between single-dataset mode and each
# dataset in multi-dataset mode.
# ---------------------------------------------------------------------------
def _render_istd_step(base_df, data_type, key_prefix, state, detected, required):
    """
    Choose the dataset's ONE internal standard, preview it and apply ISTD normalization:
    every feature is divided by the ISTD peak area in the same sample. Returns [istd] or [].
    """
    row_annotations = st.session_state.get("row_annotations")
    label = state.get("label")
    options = base_df.index.tolist()
    key = f"{key_prefix}_istd_choice"
    if key in st.session_state and st.session_state[key] not in options:
        del st.session_state[key]                         # a different table was loaded
    if key not in st.session_state:
        default = normalization.default_istd(options, row_annotations, label, data_type)
        if default is not None:
            st.session_state[key] = default

    if detected:
        st.caption(f"Internal standard detected: **{st.session_state.get(key) or detected[0]}** "
                   "(named 'ISTD…' or listed in the row annotations' **ISTD** column), pre-selected. "
                   "Every feature is divided by this ISTD's peak area in the same sample.")
    else:
        st.warning("No internal standard was detected in this dataset. Select the ISTD row below. Name it "
                   "'ISTD_…', or list it in an **ISTD** column of the row annotation file, to have it "
                   "selected automatically (see Help → Data preparation).")
    istd = st.selectbox("Internal standard (ISTD)", options, index=None, key=key,
                        placeholder="Select the internal standard row",
                        help="The ISTD row is removed from the normalized table: "
                             "normalized = feature peak area / ISTD peak area, sample by sample.")
    if istd is None:
        state["istd_normalized"] = None
        state["istd_mapping"] = None
        if required:
            st.warning("ISTD normalization is required for targeted assays — select the internal "
                       "standard above before you can proceed to Log2 transformation.")
        return []

    others = [d for d in detected if d != istd]
    mapping = normalization.build_istd_mapping([f for f in options if f not in others], [istd],
                                               row_annotations, label, data_type)
    summary = normalization.istd_summary(mapping, [istd], state.get("cv_table"), base_df)
    m1, m2, m3 = st.columns(3)
    m1.metric("Features normalized", len(mapping))
    cv = summary["QC CV (%)"].iloc[0] if "QC CV (%)" in summary.columns else np.nan
    m2.metric("ISTD QC CV", f"{cv:.1f}%" if pd.notna(cv) else "—",
              help=None if pd.notna(cv) else "Confirm QC (QC Validation) to see the ISTD's QC CV.")
    m3.metric("Zero/missing ISTD values", int(summary["Zero/missing samples"].iloc[0]))
    with st.expander(f"Metabolite → ISTD ({len(mapping)} features)"):
        st.dataframe(utils.format_df_for_display(mapping[["ISTD", "Annotated ISTD"]]), width='stretch', height=300)
        st.download_button(
            "Download ISTD_Assignments.csv", utils.to_download_bytes_csv(mapping[["ISTD", "Annotated ISTD"]]),
            file_name=f"{'' if key_prefix == 'single' else key_prefix.replace(' ', '_') + '_'}ISTD_Assignments.csv",
            mime="text/csv", key=f"{key_prefix}_dl_istd_map")
    other_named = mapping[(mapping["Annotated ISTD"] != "")
                          & (mapping["Annotated ISTD"].map(normalization._norm_key) != normalization._norm_key(istd))]
    if len(other_named):
        st.warning(f"{len(other_named)} feature(s) name a different ISTD in the row annotations "
                   f"(e.g. {other_named['Annotated ISTD'].iloc[0]}); they are normalized to {istd} here.")
    if others:
        st.caption(f"Other internal-standard rows are left out of the normalized table: {', '.join(others)}.")
    if pd.notna(cv) and cv > 20:
        st.warning(f"The QC CV of {istd} is {cv:.1f}% (above 20%): check the spiking and injection of the standard.")
    if summary["Zero/missing samples"].iloc[0] > 0:
        st.warning(f"{istd} is zero or missing in {int(summary['Zero/missing samples'].iloc[0])} sample(s); "
                   "all features become missing (NaN) in those samples.")

    if st.button("Apply ISTD Normalization", type="primary", icon="⚖️", key=f"{key_prefix}_apply_istd"):
        try:
            working, mapping = normalization.istd_normalize_mapped(
                base_df, [istd], row_annotations, label=label, method=data_type, drop_features=others)
            state["istd_normalized"] = working
            state["istd_mapping"] = mapping
            if required:
                state["log2_data"] = None                 # Step 2 must be re-run on the new ratios
            state["processing_notes"].append(
                f"ISTD normalization applied using '{istd}': each of {working.shape[0]} features divided by "
                "the ISTD peak area in the same sample; ISTD row removed.")
            st.success(f"ISTD normalization complete using {istd} ({working.shape[0]} features).")
        except Exception as e:
            st.error(str(e))
    return [istd]


def render_normalization_ui(base_df, data_type, key_prefix, state, meta):
    is_targeted = "Targeted" in data_type
    # Both untargeted analysis types (Untargeted Metabolomics AND Untargeted Lipidomics)
    # share the same strict Log2 -> Median-IQR pipeline (see the dedicated block below)
    # and never require ISTD normalization; internal-standard rows are left out of them.
    # Targeted Metabolomics and Targeted Lipidomics REQUIRE ISTD normalization before Log2
    # -- it cannot be skipped -- using the dataset's one internal standard.
    detected = _state_istds(base_df.index, state)

    if is_targeted:
        st.subheader(f"{data_type}: ISTD Normalization → Log2 Pipeline")
    else:
        st.subheader(f"{data_type}: Log2 → Median-IQR Normalization Pipeline")

    override = False
    if not is_targeted:
        override = st.checkbox(
            "Advanced: add an optional ISTD normalization step before Log2", value=False,
            key=f"{key_prefix}_override",
            help="Untargeted assays use the strict Log2 → Median-IQR pipeline. Check this to first divide "
                 "each feature by an internal standard (e.g. when standards were spiked into every "
                 "sample); the Log2 → Median-IQR pipeline below is then applied to the ISTD-normalized data."
        )

    if is_targeted:
        st.markdown("**Step 1: ISTD Normalization** *(required for targeted assays)*")
        _render_istd_step(base_df, data_type, key_prefix, state, detected, required=True)

        working = state.get("istd_normalized")
        st.markdown("**Step 2: Strict Log2 Transformation**")
        if working is None:
            st.info("Apply ISTD normalization above first (required for targeted assays) before "
                    "Log2 transformation.")
        elif st.button("Apply Strict Log2 Transformation", type="primary", icon="📐", key=f"{key_prefix}_log2_targeted_btn"):
            log2_df = stats_analysis.strict_log2(working)
            state["log2_data"] = log2_df
            state["log2_constant"] = 0.0
            n_invalid = int((pd.to_numeric(working.stack(), errors="coerce") <= 0).sum())
            state["processing_notes"].append(
                "Strict log2 transformation applied to ISTD-normalized data: log2(x), "
                "with no pseudocount, constant, or shifting."
            )
            if n_invalid:
                st.warning(
                    f"Strict log2 transformation found {n_invalid:,} non-positive value(s). "
                    "These values were set to missing (NaN); no pseudocount or shift was used."
                )
            else:
                st.success("Strict log2 transformation complete: log2(x), with no pseudocount or shift.")
            fig_dist = normalization.distribution_plots(working, log2_df)
            render_single_figure(fig_dist)

    if not is_targeted:
        use_istd = False
        if override:
            st.markdown("**Optional: ISTD Normalization**")
            use_istd = st.checkbox("Apply internal standard (ISTD) normalization", value=False,
                                   key=f"{key_prefix}_use_istd")
            if use_istd:
                chosen = _render_istd_step(base_df, data_type, key_prefix, state, detected, required=False)
            else:
                state["istd_normalized"] = None
        if use_istd and state.get("istd_normalized") is not None:
            untgt_input = state["istd_normalized"]
            st.caption("Step 1 uses the ISTD-normalized data.")
        else:
            untgt_input = normalization.drop_istds(base_df, detected)
            if detected:
                st.caption(f"{len(detected)} internal standard row(s) are left out (not biological features): "
                           f"{', '.join(map(str, detected[:8]))}{'…' if len(detected) > 8 else ''}.")
        st.markdown("**Step 1: Strict Log2 Transformation**")
        if st.button("Apply Strict Log2 Transformation", type="primary", icon="📐", key=f"{key_prefix}_log2_untgt_btn"):
            log2_only = stats_analysis.strict_log2(untgt_input)
            state["log2_only"] = log2_only
            # Also log2-transform the full (pre-CV-filter) feature table, if it
            # differs from base_df (i.e. the optional CV filter removed some
            # features) -- Step 2 uses this as the basis for per-sample
            # median/IQR statistics, regardless of which features were filtered
            # out of the table being normalized here.
            full_df = state.get("raw_peak_df_qc_full")
            if full_df is None:
                full_df = base_df
            if use_istd and state.get("istd_normalized") is not None:
                full_df, _ = normalization.istd_normalize_mapped(
                    full_df, [c for c in chosen if c in full_df.index], st.session_state.get("row_annotations"),
                    label=state.get("label"), method=data_type,
                    drop_features=[d for d in detected if d not in chosen])
            else:
                full_df = normalization.drop_istds(full_df, detected)
            state["log2_full"] = stats_analysis.strict_log2(full_df)
            numeric_base = untgt_input.apply(pd.to_numeric, errors="coerce")
            n_invalid = int((numeric_base <= 0).sum().sum())
            state["processing_notes"].append(
                ("Strict log2 transformation applied to ISTD-normalized data: "
                 if use_istd and state.get("istd_normalized") is not None else
                 "Strict log2 transformation applied to raw (QC-validated) peak areas: ")
                + "log2(x), with no pseudocount, constant, or shifting."
                + (f" {len(detected)} internal standard row(s) excluded." if detected and not use_istd else "")
            )
            if n_invalid:
                st.warning(
                    f"Strict log2 transformation found {n_invalid:,} non-positive raw value(s). "
                    "These values were set to missing (NaN); no pseudocount or shift was used."
                )
            else:
                st.success("Strict log2 transformation complete: log2(x), with no pseudocount or shift.")
            state["fig_dist_log2"] = normalization.distribution_plots(untgt_input, log2_only)

        st.markdown("**Step 2: Median-IQR Normalization (Robust Scaling)**")
        iqr_axis = st.radio("Normalization axis", ["sample", "feature", "batch"], horizontal=True,
                             index=0, key=f"{key_prefix}_iqr_axis")
        if state.get("log2_only") is None:
            st.info("Run Step 1 (Strict Log2 Transformation) first.")
        elif st.button("Apply Median-IQR Normalization", type="primary", icon="⚖️", key=f"{key_prefix}_iqr_untgt_btn"):
            batch_map = meta["Batch"] if iqr_axis == "batch" else None
            stats_source = state.get("log2_full") if iqr_axis == "sample" else None
            log2_full = state.get("log2_full")
            n_stats_features = (
                log2_full.shape[0] if stats_source is not None and log2_full is not None
                else state["log2_only"].shape[0]
            )
            log2_iqr = normalization.iqr_normalize(
                state["log2_only"], axis=iqr_axis, batch_map=batch_map, avoid_nan=True,
                stats_source=stats_source
            )
            state["log2_data"] = log2_iqr
            state["log2_constant"] = 0.0
            state["fig_dist_iqr"] = normalization.distribution_plots(state["log2_only"], log2_iqr)
            state["processing_notes"].append(
                f"Median-IQR robust scaling applied to Log2-transformed data ({iqr_axis}-based), "
                "with zero-spread features left centered at 0 rather than converted to missing."
                + (f" median_j/IQR_j computed from all {n_stats_features} detected features per "
                   "sample (pre-CV-filter), applied to the "
                   f"{state['log2_only'].shape[0]} features kept for this normalized matrix."
                   if iqr_axis == "sample" else "")
            )
            st.success(f"Median-IQR normalization complete ({iqr_axis}-based). "
                       "This is the final normalized matrix for this dataset.")

        if state.get("fig_dist_log2") is not None or state.get("fig_dist_iqr") is not None:
            st.markdown("**Distribution Diagnostics**")
            if state.get("fig_dist_log2") is not None:
                st.caption("Row 1 — After Log2 Transformation (raw vs. log2, distribution + box plot)")
                render_single_figure(
                    state["fig_dist_log2"],
                    download_name=f"{key_prefix.replace(' ', '_')}_Distribution_AfterLog2",
                    key_prefix=f"{key_prefix}_distlog2", width_ratio=None
                )
            if state.get("fig_dist_iqr") is not None:
                st.caption("Row 2 — After Median-IQR Normalization (log2 vs. log2+IQR, distribution + box plot)")
                render_single_figure(
                    state["fig_dist_iqr"],
                    download_name=f"{key_prefix.replace(' ', '_')}_Distribution_AfterIQR",
                    key_prefix=f"{key_prefix}_distiqr", width_ratio=None
                )

    if is_targeted:
        prelog2_data = state.get("istd_normalized")
        log2_only_data = None
    else:
        prelog2_data = None
        log2_only_data = state.get("log2_only")

    if prelog2_data is not None or log2_only_data is not None or state.get("log2_data") is not None:
        st.markdown("**Current Working Matrix**")
        file_prefix = "" if key_prefix == "single" else f"{key_prefix.replace(' ', '_')}_"

        if prelog2_data is not None:
            st.markdown("*Normalized — Without Log2 Transformation*")
            st.dataframe(utils.format_df_for_display(prelog2_data), width='stretch', height=300)
            st.download_button(
                f"Download {file_prefix}Normalized_NoLog2.csv",
                utils.to_download_bytes_csv(prelog2_data),
                file_name=f"{file_prefix}Normalized_NoLog2.csv", mime="text/csv",
                key=f"{key_prefix}_dl_norm_nolog2"
            )
        elif log2_only_data is not None:
            st.markdown("*Log2-Transformed — Without Median-IQR Normalization*")
            st.dataframe(utils.format_df_for_display(log2_only_data), width='stretch', height=300)
            st.download_button(
                f"Download {file_prefix}Log2_Only.csv",
                utils.to_download_bytes_csv(log2_only_data),
                file_name=f"{file_prefix}Log2_Only.csv", mime="text/csv",
                key=f"{key_prefix}_dl_log2_only"
            )
        else:
            st.info("Run Step 1 above to enable the intermediate download.")

        if state.get("log2_data") is not None:
            st.markdown("*Final: Log2 + Median-IQR Normalized*" if not is_targeted
                         else "*Normalized — With Log2 Transformation*")
            st.dataframe(utils.format_df_for_display(state["log2_data"]), width='stretch', height=300)
            st.download_button(
                f"Download {file_prefix}Normalized_Log2.csv",
                utils.to_download_bytes_csv(state["log2_data"]),
                file_name=f"{file_prefix}Normalized_Log2.csv", mime="text/csv",
                key=f"{key_prefix}_dl_norm_log2"
            )
        else:
            st.info("Run Step 2 above to enable the final download.")


# ===========================================================================
# TAB 5 — NORMALIZATION
# ===========================================================================
def page_normalization():
    if st.session_state.data_mode == "single":
        if st.session_state.raw_peak_df is None:
            st.warning("Load data on the Data Upload page first.")
        else:
            if st.session_state.get("raw_peak_df_qc") is not None:
                base_df = st.session_state.raw_peak_df_qc
            else:
                base_df = st.session_state.raw_peak_df[st.session_state.sample_cols]
                st.info("QC wasn't confirmed yet — proceeding with biological samples only "
                         "(QC columns excluded automatically). Visit QC Validation to review QC first.")
            render_normalization_ui(base_df, st.session_state.data_type, "single",
                                     _SingleDatasetState(), st.session_state.meta)
    else:
        if not st.session_state.datasets:
            st.warning("Load datasets on the Data Upload page first.")
        else:
            selected_id = render_dataset_dropdown(st.session_state.datasets, "norm", "3. Normalization")
            ds = st.session_state.datasets[selected_id]
            st.divider()
            st.subheader(f"Configuring: {ds['label']} ({ds['data_type']})")
            if ds.get("raw_peak_df_qc") is not None:
                base_df = ds["raw_peak_df_qc"]
            else:
                base_df = ds["raw_df"][ds["sample_cols"]]
                st.info(f"QC wasn't confirmed yet for {ds['label']} — proceeding with "
                         "biological samples only.")
            render_normalization_ui(base_df, ds["data_type"], selected_id, ds, st.session_state.meta)
            st.caption("Once every dataset is normalized, combine them in **Combined Normalization Data (Multi-Dataset)**.")

# ===========================================================================
# TAB 6 — COMBINED NORMALIZED DATA
# ===========================================================================
def page_combined_normalized():
    st.header("Combined Normalization Data (Multi-Dataset)")
    if st.session_state.data_mode == "single":
        st.info("Combining only applies in multi-dataset mode (Data Upload). In single-dataset mode, "
                "Normalization's output already feeds every downstream page directly.")
    elif not st.session_state.datasets:
        st.warning("Load datasets on the Data Upload page first.")
    else:
        st.caption(
            "Combines every dataset's normalized matrix into one table (metabolite names "
            "prefixed by dataset, e.g. `Targeted Metabolomics::Glucose`) — the single input for "
            "every downstream page."
        )
        summary = dataset_manager.normalization_summary_table(st.session_state.datasets)
        st.dataframe(utils.format_df_for_display(summary), width='stretch')

        n_ready = (summary["Status"] == "Complete").sum()
        n_total = len(summary)
        if n_ready < n_total:
            st.warning(
                f"⏳ {n_ready} of {n_total} datasets are normalized. **Generate Combined "
                f"Normalized Data** stays disabled until all {n_total} are complete — finish "
                "each remaining dataset's Log2 Transformation step in Normalization (select it from "
                "the dropdown), or use the 🚀 Quick Start button in Data Upload to normalize all of "
                "them at once."
            )
        generate_clicked = st.button(
            "🔗 Generate Combined Normalized Data", type="primary", disabled=(n_ready < n_total)
        )
        if generate_clicked:
            # Intermediate reference table: each dataset's matrix from immediately before
            # its own final pipeline step (pre-Log2 for pipelines that end in Log2, or
            # Log2-only-pre-IQR for the Untargeted Metabolomics pipeline, which runs Log2
            # before IQR-normalization). Kept only as a downloadable reference alongside
            # the final combined matrix below.
            combined_prelog2, prelog2_counts = dataset_manager.combine_normalized_datasets(
                st.session_state.datasets, st.session_state.sample_cols, use_log2=False
            )
            # Build the final combined matrix directly from each dataset's own already-
            # finalized 'log2_data' (normalized AND Log2-transformed, in whichever order
            # that dataset's pipeline applies them). Log2 is computed elementwise per
            # dataset in Normalization, independently of every other row/dataset, so combining
            # each dataset's finished matrix here gives the same result as combining
            # pre-Log2 matrices and Log2-transforming them together -- while also
            # correctly supporting pipelines (Untargeted Metabolomics) where Log2 runs
            # before the final normalization step, for which re-applying Log2 to an
            # already Log2-transformed, IQR-centered (possibly negative) combined table
            # would be incorrect.
            combined, counts = dataset_manager.combine_normalized_datasets(
                st.session_state.datasets, st.session_state.sample_cols, use_log2=True
            )
            st.session_state.log2_data = combined
            st.session_state.combined_normalized_prelog2 = combined_prelog2
            st.session_state.combined_done = True
            reset_downstream_analysis_state()

            # Aggregate each dataset's QC figures into the flat structure downstream
            # exports read from, prefixed by dataset label.
            combined_qc_figs = {}
            for ds2 in st.session_state.datasets.values():
                for fig_name, fig in ds2.get("qc_figs", {}).items():
                    combined_qc_figs[f"{ds2['label']}: {fig_name}"] = fig
            st.session_state.qc_figs = combined_qc_figs

            st.session_state.processing_notes.append(
                f"Generated combined normalized data from {len(counts)} independently-"
                f"normalized datasets into one final Log2 matrix ({combined.shape[0]} total "
                "features): " + ", ".join(f"{lbl} ({n})" for lbl, n in counts.items()) + ". "
                f"An intermediate pre-final-step combined matrix "
                f"({combined_prelog2.shape[0] if combined_prelog2 is not None else 0} "
                "total features) was also generated for reference. This combined table is "
                "now the single input for every downstream analysis page."
            )
            st.success(f"Combined into {combined.shape[0]} total features across "
                       f"{len(counts)} datasets. Every downstream page (PCA, Statistics, "
                       "Volcano, Biomarker Discovery, Heatmap, Boxplot) now uses this table.")

        if st.session_state.get("combined_done") and st.session_state.log2_data is not None:
            st.markdown("### Combined Normalized Data Table")
            st.markdown("**Normalized Without Log2 Transformation**")
            prelog2 = st.session_state.get("combined_normalized_prelog2")
            if prelog2 is not None:
                st.dataframe(utils.format_df_for_display(prelog2), width='stretch', height=400)
                st.download_button(
                    "Download Combined_Normalized_NoLog2.csv",
                    utils.to_download_bytes_csv(prelog2),
                    file_name="Combined_Normalized_NoLog2.csv", mime="text/csv",
                    key="dl_combined_norm_nolog2"
                )
            else:
                st.info("No pre-Log2 normalized data available to combine yet.")

            st.markdown("**Normalized With Log2 Transformation**")
            st.dataframe(utils.format_df_for_display(st.session_state.log2_data), width='stretch', height=400)
            st.download_button(
                "Download Combined_Normalized_Strict_Log2.csv",
                utils.to_download_bytes_csv(st.session_state.log2_data),
                file_name="Combined_Normalized_Strict_Log2.csv", mime="text/csv",
                key="dl_combined_norm_log2"
            )
            st.caption(
                f"{st.session_state.log2_data.shape[0]} total features (dataset-prefixed) × "
                f"{st.session_state.log2_data.shape[1]} samples. This is the table every "
                "downstream page now uses — head to **PCA** to continue."
            )

# ===========================================================================
# TAB 8 — STATISTICS
# ===========================================================================
def page_statistics():
    st.header("Statistical Comparison")
    if st.session_state.log2_data is None:
        st.warning(
            "Complete Normalization (strict log2(x) transformation) first."
            if st.session_state.data_mode == "single" else
            "Complete Normalization for every dataset, then click **🔗 Generate Combined "
            "Normalized Data** first — that combined table is what every "
            "downstream page (PCA, Statistics, Volcano, Biomarker Discovery, Heatmap, "
            "Boxplot) analyzes."
        )
    else:
        if st.session_state.data_mode == "multi":
            st.caption("Analyzing the **combined normalized data** across all datasets.")
        fdr_label = current_fdr_label()
        log2_df = utils.with_display_feature_names(st.session_state.log2_data)
        meta = st.session_state.meta
        sample_cols = log2_df.columns.tolist()
        stats_cat_cols = utils.get_categorical_metadata_columns(meta, sample_cols)
        grouping_var = st.selectbox(
            "Grouping variable:", stats_cat_cols,
            index=stats_cat_cols.index("Group") if "Group" in stats_cat_cols else 0,
            key="stats_grouping_var"
        )
        groups_available = meta.loc[meta.index.intersection(sample_cols), grouping_var].unique().tolist()

        mode = st.radio("Comparison type", ["Two-group comparison", "ANOVA (≥3 groups)"], horizontal=True)

        if mode == "Two-group comparison":
            c1, c2, c3 = st.columns(3)
            default_a = groups_available.index("Untreated") if "Untreated" in groups_available else 0
            default_b = groups_available.index("IR Day 2") if "IR Day 2" in groups_available else min(1, len(groups_available) - 1)
            group_a = c1.selectbox("Group A", groups_available, index=default_a)
            group_b = c2.selectbox("Group B", groups_available, index=default_b)
            method = c3.selectbox("Method", ["Welch's t-test", "Wilcoxon rank-sum"])
            method_key = "ttest" if method.startswith("Welch") else "wilcoxon"

            if st.button("Run Two-Group Test", type="primary", icon="📊"):
                a_samples = meta.index[meta[grouping_var] == group_a].tolist()
                b_samples = meta.index[meta[grouping_var] == group_b].tolist()
                a_samples = [s for s in a_samples if s in log2_df.columns]
                b_samples = [s for s in b_samples if s in log2_df.columns]

                result = stats_analysis.two_group_test(
                    log2_df, a_samples, b_samples, method=method_key, fdr_label=fdr_label
                )
                st.session_state.stats_result = result
                st.session_state.stats_result_groups = (group_a, group_b)
                st.session_state.stats_result_group_col = grouping_var
                n_significant = int(
                    ((result["p-value"] < 0.05) & (result[fdr_label] < 0.25)).sum()
                )
                st.session_state.processing_notes.append(
                    f"Statistical comparison ({grouping_var}): {group_a} vs {group_b} using {method} "
                    f"(Benjamini-Hochberg p.adjust(method=\"BH\") correction, shown as '{fdr_label}'); "
                    f"{n_significant} significant metabolites (p<0.05 & {fdr_label}<0.25)."
                )
                st.success(f"Test complete: {n_significant} significant metabolites found.")

            if st.session_state.stats_result is not None and st.session_state.stats_result_groups is not None:
                result = st.session_state.stats_result
                g_a, g_b = st.session_state.stats_result_groups
                comparison_name = f"{g_a}_vs_{g_b}"
                st.dataframe(utils.format_df_for_display(result), width='stretch', height=400)
                st.download_button(
                    f"Download {comparison_name}_Statistics.csv",
                    utils.to_download_bytes_csv(result),
                    file_name=f"{comparison_name}_Statistics.csv", mime="text/csv",
                    key="dl_stats_twogroup"
                )

        else:
            if len(groups_available) < 3:
                st.warning(f"Need at least 3 values of '{grouping_var}' (excluding QC) for ANOVA. "
                           "Pick a different grouping variable, or add more groups in metadata.")
            else:
                anova_groups = st.multiselect(
                    f"{grouping_var} values to include in ANOVA (select 3 or more)",
                    groups_available, default=groups_available, key="anova_group_select"
                )
                posthoc_method = st.selectbox("Post-hoc test", ["tukey", "dunnett", "pairwise"])

                if len(anova_groups) < 3:
                    st.warning(f"Select at least 3 groups to run ANOVA (currently {len(anova_groups)} selected).")
                elif st.button("Run ANOVA", type="primary", icon="📊"):
                    group_map = meta.loc[sample_cols, grouping_var]
                    group_map = group_map[group_map.isin(anova_groups)]
                    group_map = group_map[group_map.index.isin(log2_df.columns)]
                    anova_table, posthoc_results = stats_analysis.anova_test(
                        log2_df[group_map.index], group_map, posthoc=posthoc_method,
                        fdr_label=fdr_label
                    )
                    st.session_state.anova_result = anova_table
                    st.session_state.posthoc_result = posthoc_results
                    st.session_state.anova_groups_used = list(anova_groups)
                    st.session_state.stats_result_group_col = grouping_var
                    n_sig_anova = int((anova_table[fdr_label] < 0.25).sum())
                    st.session_state.processing_notes.append(
                        f"One-way ANOVA across {len(anova_groups)} selected {grouping_var} values "
                        f"({', '.join(anova_groups)}) with {posthoc_method} post-hoc; "
                        f"{n_sig_anova} significant metabolites ({fdr_label}<0.25)."
                    )
                    st.success(f"ANOVA complete: {n_sig_anova} significant metabolites ({fdr_label}<0.25).")

                if st.session_state.anova_result is not None:
                    anova_table = st.session_state.anova_result
                    anova_groups_used = st.session_state.get("anova_groups_used") or anova_groups
                    anova_comparison_name = "ANOVA_" + "_vs_".join(anova_groups_used)
                    st.dataframe(utils.format_df_for_display(anova_table), width='stretch', height=350)
                    st.download_button(
                        f"Download {anova_comparison_name}.csv",
                        utils.to_download_bytes_csv(anova_table),
                        file_name=f"{anova_comparison_name}.csv", mime="text/csv",
                        key="dl_stats_anova"
                    )

                    posthoc_results = st.session_state.get("posthoc_result")
                    if posthoc_results:
                        feat_choice = st.selectbox("View post-hoc results for feature:", list(posthoc_results.keys()))
                        st.dataframe(utils.format_df_for_display(posthoc_results[feat_choice]), width='stretch')
                        st.download_button(
                            f"Download {anova_comparison_name}_Posthoc_{feat_choice}.csv",
                            utils.to_download_bytes_csv(posthoc_results[feat_choice]),
                            file_name=f"{anova_comparison_name}_Posthoc_{feat_choice}.csv", mime="text/csv",
                            key="dl_stats_posthoc_one"
                        )
                        all_posthoc = pd.concat(
                            [df.assign(Feature=feat) for feat, df in posthoc_results.items()],
                            ignore_index=True
                        )
                        st.download_button(
                            f"Download {anova_comparison_name}_Posthoc_All_Features.csv",
                            utils.to_download_bytes_csv(all_posthoc),
                            file_name=f"{anova_comparison_name}_Posthoc_All_Features.csv", mime="text/csv",
                            key="dl_stats_posthoc_all"
                        )

# ===========================================================================
# TAB 7 — PCA
# ===========================================================================
def page_pca():
    st.header("PCA Visualization (Biological Samples Only)")
    if st.session_state.log2_data is None:
        st.warning(
            "Complete Normalization (strict log2(x) transformation) first."
            if st.session_state.data_mode == "single" else
            "Complete Normalization for every dataset, then click **🔗 Generate Combined "
            "Normalized Data** first — that combined table is what every "
            "downstream page analyzes."
        )
    else:
        if st.session_state.data_mode == "multi":
            st.caption("Analyzing the **combined normalized data** across all datasets.")
        log2_df_full = utils.with_display_feature_names(st.session_state.log2_data)
        meta = st.session_state.meta
        pca_cat_cols = utils.get_categorical_metadata_columns(meta, log2_df_full.columns)
        pca_group_col = st.selectbox(
            "Color/group samples by:", pca_cat_cols,
            index=pca_cat_cols.index("Group") if "Group" in pca_cat_cols else 0,
            key="pca_group_col"
        )
        groups_all_pca = meta.loc[meta.index.intersection(log2_df_full.columns), pca_group_col].unique().tolist()

        st.subheader("Select Groups")
        selected_pca_groups = st.multiselect(
            f"{pca_group_col} values to include in PCA (choose any subset — 2, 3, 4, or more)",
            groups_all_pca, default=groups_all_pca, key="pca_group_select"
        )
        pca_cols = [c for c in log2_df_full.columns if meta.loc[c, pca_group_col] in selected_pca_groups]
        log2_df = log2_df_full[pca_cols]

        if len(selected_pca_groups) < 1 or len(pca_cols) < 3:
            st.warning("Select at least one group with enough samples (≥3 total) to run PCA.")
        else:
            n_comp = st.slider("Number of components", 2, min(10, log2_df.shape[1] - 1 if log2_df.shape[1] > 2 else 2), 5)

            st.subheader("Customize Appearance")
            c1, c2 = st.columns(2)
            palette_choice = c1.selectbox("Color palette", list(pca_module.PALETTES.keys()))
            show_ellipse = c2.checkbox("Show 95% confidence ellipses", value=True)

            st.caption("Optional: assign a marker style per group (defaults to circles for all).")
            marker_map = {}
            marker_cols = st.columns(min(4, len(selected_pca_groups)) or 1)
            for i, g in enumerate(selected_pca_groups):
                with marker_cols[i % len(marker_cols)]:
                    marker_map[g] = st.selectbox(f"Marker: {g}", pca_module.MARKER_STYLES, key=f"marker_{g}")

            pca, scores_df, cols = pca_module.run_pca(log2_df, n_components=n_comp)
            fig_score = pca_module.pca_score_plot(pca, scores_df, meta, palette=palette_choice,
                                                    marker_map=marker_map, show_ellipse=show_ellipse,
                                                    group_col=pca_group_col)
            fig_loading, top_loadings = pca_module.pca_loading_plot(pca, log2_df.index, top_n=20)
            fig_var = pca_module.pca_variance_plot(pca)

            st.session_state.viz_figs["PCA Score Plot"] = fig_score
            st.session_state.viz_figs["PCA Loading Plot"] = fig_loading
            st.session_state.viz_figs["PCA Variance Plot"] = fig_var

            render_figure_row([
                {"title": "PCA Score Plot", "fig": fig_score,
                 "download_name": "PCA_Score_Plot", "key_prefix": "pca_score"},
                {"title": "PCA Loading Plot", "fig": fig_loading,
                 "download_name": "PCA_Loading_Plot", "key_prefix": "pca_loading"},
                {"title": "PCA Variance Plot", "fig": fig_var,
                 "download_name": "PCA_Variance_Plot", "key_prefix": "pca_var"},
            ])

            st.subheader("Top Contributing Metabolites (Loadings)")
            st.dataframe(utils.format_df_for_display(top_loadings), width='stretch')

            st.divider()
            st.subheader("Download PCA Data")
            st.caption(
                "The underlying data behind the plots above: sample scores per principal "
                "component (with group assignment), % variance explained per component, "
                "the metabolite loadings, and the exact normalized data matrix PCA was run on."
            )
            explained_var_df = pd.DataFrame({
                "Component": [f"PC{i+1}" for i in range(len(pca.explained_variance_ratio_))],
                "Variance Explained (%)": (pca.explained_variance_ratio_ * 100).round(3),
                "Cumulative Variance Explained (%)": (np.cumsum(pca.explained_variance_ratio_) * 100).round(3),
            })
            scores_with_group = scores_df.copy()
            scores_with_group.insert(0, pca_group_col, meta.loc[scores_with_group.index, pca_group_col].values)

            pca_tag = "_".join(_safe_tag(g) for g in selected_pca_groups) or "AllGroups"
            if st.session_state.data_mode == "multi":
                pca_tag = f"Combined_{pca_tag}"
            dl1, dl2, dl3, dl4 = st.columns(4)
            dl1.download_button(
                f"Download {pca_tag}_PCA_Scores.csv", utils.to_download_bytes_csv(scores_with_group),
                file_name=f"{pca_tag}_PCA_Scores.csv", mime="text/csv", key="dl_pca_scores"
            )
            dl2.download_button(
                f"Download {pca_tag}_PCA_Explained_Variance.csv", utils.to_download_bytes_csv(explained_var_df),
                file_name=f"{pca_tag}_PCA_Explained_Variance.csv", mime="text/csv", key="dl_pca_variance"
            )
            dl3.download_button(
                f"Download {pca_tag}_PCA_Loadings.csv", utils.to_download_bytes_csv(top_loadings),
                file_name=f"{pca_tag}_PCA_Loadings.csv", mime="text/csv", key="dl_pca_loadings"
            )
            dl4.download_button(
                f"Download {pca_tag}_PCA_Input_Data.csv", utils.to_download_bytes_csv(log2_df),
                file_name=f"{pca_tag}_PCA_Input_Data.csv", mime="text/csv", key="dl_pca_input"
            )

# ===========================================================================
# TAB 9 — VOLCANO PLOT
# ===========================================================================
def page_volcano():
    st.header("Volcano Plot")
    if st.session_state.stats_result is None:
        st.warning("Run a two-group statistical comparison in Statistics first.")
    else:
        result = st.session_state.stats_result
        fdr_label = stats_fdr_label(result)
        groups_used = st.session_state.get("stats_result_groups")
        group_col_used = st.session_state.get("stats_result_group_col") or "Group"
        if groups_used:
            st.caption(
                f"Reflects the two {group_col_used} values compared in Statistics: "
                f"**{groups_used[0]}** vs **{groups_used[1]}**"
                + (" (combined normalized data, all datasets)."
                   if st.session_state.data_mode == "multi" else ".")
                + " To compare a different pair or variable, go back to Statistics and "
                "re-run the two-group test."
            )

        # ---- Quick Style: the four settings that do the most to make the figure look
        # publication-ready, surfaced up front instead of buried in an expander below. ----
        st.markdown("##### 🎨 Quick Style")
        qs1, qs2, qs3, qs4 = st.columns(4)
        theme_choice = qs1.selectbox("Journal theme", list(volcano.THEMES.keys()),
                                      help="Stylistic approximation (font, grid, spine, palette), not an "
                                           "official journal figure spec.")
        palette = qs2.selectbox("Color palette", list(volcano.COLORBLIND_PALETTES.keys()))
        legend_position = qs3.selectbox("Legend position", volcano.LEGEND_POSITION_CHOICES, index=1,
                                        help="Right / Left: pick top, middle or bottom alignment beside the plot. "
                                             "Top / Bottom: centered above or below the plot.")
        size_preset = qs4.selectbox("Figure size", list(volcano.FIGURE_SIZE_PRESETS.keys()), index=2)
        if size_preset == "Custom":
            cw, ch = st.columns(2)
            fig_width_in = cw.number_input("Width (in)", value=7.2, min_value=2.0, max_value=20.0)
            fig_height_in = ch.number_input("Height (in)", value=6.4, min_value=2.0, max_value=20.0)
        else:
            fig_width_in, fig_height_in = volcano.FIGURE_SIZE_PRESETS[size_preset]

        # ---- Statistical thresholds ----
        with st.expander("🎯 Statistical Thresholds", expanded=True):
            c1, c2, c3, c4 = st.columns(4)
            y_metric_label = c1.radio("Y-axis metric", ["p-value", fdr_label], horizontal=False)
            y_metric = "pvalue" if y_metric_label == "p-value" else "fdr"
            sig_cutoff = c2.number_input("Significance cutoff", value=0.05 if y_metric == "pvalue" else 0.25,
                                          min_value=0.0001, max_value=1.0, step=0.01)
            fc_threshold = c3.number_input("Fold change threshold (log2 units)", value=0.0,
                                            min_value=0.0, max_value=5.0, step=0.1)
            use_fdr_secondary = c4.checkbox(f"Also require {fdr_label} <", value=False)
            fdr_cutoff = c4.number_input(f"Secondary {fdr_label} cutoff", value=0.25, min_value=0.0001, max_value=1.0,
                                          step=0.01, disabled=not use_fdr_secondary) if use_fdr_secondary else None
            c5, c6 = st.columns(2)
            use_fc_axis = c5.checkbox("Show linear FC on x-axis instead of log2FC", value=False)
            threshold_line_style = c6.selectbox("Threshold line style", list(volcano.LINE_STYLES.keys()))

        # ---- Point style & color overrides ----
        with st.expander("⚫ Point Style & Color Overrides"):
            override_colors = st.checkbox("Override individual colors (default: follow the palette above)",
                                           value=False)
            up_color = down_color = ns_color = None
            if override_colors:
                cc1, cc2, cc3 = st.columns(3)
                up_color = cc1.color_picker("Upregulated color", volcano.COLORBLIND_PALETTES[palette]["up"])
                down_color = cc2.color_picker("Downregulated color", volcano.COLORBLIND_PALETTES[palette]["down"])
                ns_color = cc3.color_picker("Not significant color", volcano.COLORBLIND_PALETTES[palette]["ns"])
            c3, c4, c5 = st.columns(3)
            alpha = c3.slider("Point transparency (alpha)", 0.2, 1.0, 0.85, 0.05)
            point_size = c4.slider("Point size", 1, 100, 18, 1)
            point_shape = c5.selectbox("Point shape", list(volcano.MARKER_SHAPES.keys()))
            c6, c7 = st.columns(2)
            edge_width = c6.slider("Point border width", 0.0, 2.0, 0.3, 0.1)
            edge_color = c7.color_picker("Point border color", "#333333") if edge_width > 0 else "none"

        # ---- Labels ----
        with st.expander("🏷️ Metabolite Labels"):
            label_mode_label = st.radio(
                "Label mode",
                ["Top N significant", "All significant", "Manually selected", "None"],
                horizontal=True
            )
            label_mode = {"Top N significant": "top_n", "All significant": "significant_only",
                          "Manually selected": "manual", "None": "none"}[label_mode_label]
            top_label_n = 10
            manual_labels = []
            if label_mode == "top_n":
                top_label_n = st.number_input("Label top N metabolites", value=10, min_value=0, max_value=100, step=1)
            elif label_mode == "manual":
                manual_labels = st.multiselect("Search and select metabolites to label by name",
                                                result.index.tolist())
                st.caption("Note: labeling by HMDB ID isn't available — this dataset identifies "
                           "metabolites by name only (no compound-database annotation).")
            c1, c2, c3 = st.columns(3)
            label_font_size = c1.slider("Label font size", 5.0, 16.0, 7.5, 0.5)
            label_bold = c2.checkbox("Bold labels", value=False)
            label_italic = c3.checkbox("Italic labels", value=False)
            c4, c5 = st.columns(2)
            override_label_color = c4.checkbox("Override label color (default: match point color)", value=False)
            label_color = c5.color_picker("Label color", "#333333") if override_label_color else None
            repel_labels = st.checkbox("Repel overlapping labels automatically", value=True)

        # ---- Legend label format (position is set above, under Quick Style) ----
        with st.expander("📋 Legend Label Format"):
            legend_format = st.radio("Legend label format",
                                      ["Full (\"Significantly Upregulated (n=5)\")", "Short (\"Up (5)\")"], horizontal=False)
            legend_format_key = "full" if legend_format.startswith("Full") else "short"

        # ---- Axes and title ----
        with st.expander("📐 Axes & Title"):
            st.markdown("**Axes**")
            c1, c2, c3 = st.columns(3)
            axis_font_size = c1.slider("Axis font size", 8.0, 20.0, 12.0, 0.5)
            axis_bold = c2.checkbox("Bold axis titles", value=False)
            y_decimals = c3.selectbox("Y-axis decimals", [None, 0, 1, 2, 3], index=0,
                                      format_func=lambda x: "Auto" if x is None else str(x))
            c4, c5 = st.columns(2)
            custom_x_limits = c4.checkbox("Set custom X-axis limits", value=False)
            x_limits = None
            if custom_x_limits:
                xlo, xhi = st.columns(2)
                x_limits = (xlo.number_input("X min", value=-6.0), xhi.number_input("X max", value=6.0))
            custom_y_limits = c5.checkbox("Set custom Y-axis limits", value=False)
            y_limits = None
            if custom_y_limits:
                ylo, yhi = st.columns(2)
                y_limits = (ylo.number_input("Y min", value=0.0), yhi.number_input("Y max", value=10.0))
            x_tick_spacing = st.number_input("X-axis tick spacing (0 = auto)", value=0.0, min_value=0.0, step=0.5)
            x_tick_spacing = x_tick_spacing or None

            st.markdown("**Title**")
            c6, c7 = st.columns(2)
            default_title = f"{groups_used[0]} vs {groups_used[1]}" if groups_used else "Volcano Plot"
            custom_title = c6.text_input("Title", value=default_title)
            subtitle = c7.text_input("Subtitle (optional)", value="")
            c8, c9, c10, c11 = st.columns(4)
            title_bold = c8.checkbox("Bold title", value=True)
            title_italic = c9.checkbox("Italic title", value=False)
            title_align = c10.selectbox("Title alignment", ["center", "left", "right"])
            hide_title = c11.checkbox("Hide title", value=False)

        # ---- Gridlines and threshold line style ----
        with st.expander("▦ Gridlines & Threshold Line Style"):
            c1, c2, c3 = st.columns(3)
            grid_mode = c1.selectbox("Gridlines", ["none", "major", "both"],
                                     format_func=lambda x: {"none": "No grid", "major": "Major only",
                                                              "both": "Major + minor"}[x])
            grid_color = c2.color_picker("Grid color", "#D9D9D9")
            grid_style = c3.selectbox("Grid line style", list(volcano.LINE_STYLES.keys()), index=0)
            c4, c5 = st.columns(2)
            threshold_line_color = c4.color_picker("Threshold line color", "#808080")
            threshold_line_width = c5.slider("Threshold line width", 0.2, 3.0, 0.7, 0.1)

        # ---- Highlight specific metabolites ----
        with st.expander("⭐ Highlight Specific Metabolites"):
            highlight_names = st.multiselect(
                "Search and select metabolites to highlight (e.g. known markers)",
                result.index.tolist()
            )
            c1, c2, c3 = st.columns(3)
            highlight_color = c1.color_picker("Highlight color", "#FFD700")
            highlight_size_mult = c2.slider("Highlight size multiplier", 1.0, 5.0, 2.0, 0.25)
            highlight_shape = c3.selectbox("Highlight shape", ["Star", "Circle", "Triangle", "Square", "Diamond"])

        # ---- Background (figure size is set above, under Quick Style) ----
        with st.expander("🖼️ Background"):
            background_choice = st.selectbox("Background", ["White", "Transparent", "Gray", "Custom"])
            background = {"White": "white", "Transparent": "transparent", "Gray": "gray"}.get(background_choice)
            if background_choice == "Custom":
                background = st.color_picker("Custom background color", "#FFFFFF")

        # ---- Statistics overlay (theme is set above, under Quick Style) ----
        with st.expander("📊 Statistics Overlay"):
            show_stats_box = st.checkbox("Show statistics box on figure (total/up/down/cutoffs)", value=False)

        fig_volc, annotated = volcano.volcano_plot(
            result,
            use_fc_not_log2=use_fc_axis, fc_threshold=fc_threshold, y_metric=y_metric,
            sig_cutoff=sig_cutoff, fdr_cutoff=fdr_cutoff, threshold_line_style=threshold_line_style,
            palette=palette, up_color=up_color, down_color=down_color, ns_color=ns_color,
            alpha=alpha, edge_color=edge_color, edge_width=edge_width,
            point_size=point_size, point_shape=point_shape,
            label_mode=label_mode, top_label_n=top_label_n, manual_labels=manual_labels,
            label_font_size=label_font_size, label_color=label_color, label_bold=label_bold,
            label_italic=label_italic, repel_labels=repel_labels,
            legend_position=legend_position, legend_format=legend_format_key,
            axis_font_size=axis_font_size, axis_bold=axis_bold,
            x_limits=x_limits, y_limits=y_limits, x_tick_spacing=x_tick_spacing, y_decimals=y_decimals,
            title=custom_title, subtitle=subtitle or None, title_bold=title_bold, title_italic=title_italic,
            title_align=title_align, hide_title=hide_title,
            grid_mode=grid_mode, grid_style=grid_style, grid_color=grid_color,
            threshold_line_color=threshold_line_color, threshold_line_width=threshold_line_width,
            highlight_names=highlight_names, highlight_color=highlight_color,
            highlight_size_mult=highlight_size_mult, highlight_shape=highlight_shape,
            background=background, fig_width_in=fig_width_in, fig_height_in=fig_height_in,
            show_stats_box=show_stats_box, theme=theme_choice, fdr_label=fdr_label,
        )
        st.session_state.viz_figs[f"Volcano Plot ({y_metric_label})"] = fig_volc
        st.session_state.volcano_fig = fig_volc
        st.session_state.volcano_annotated = annotated
        st.session_state.volcano_settings = {
            "y_metric": y_metric, "sig_cutoff": sig_cutoff, "fc_threshold": fc_threshold,
            "fdr_cutoff": fdr_cutoff, "palette": palette, "point_size": point_size,
            "point_shape": point_shape, "label_mode": label_mode, "legend_position": legend_position,
            "theme": theme_choice, "groups": groups_used,
        }
        render_single_figure(fig_volc, width_ratio=(1, 6, 1))

        volcano_tag = (f"{_safe_tag(groups_used[0])}_vs_{_safe_tag(groups_used[1])}"
                       if groups_used else "Comparison")
        if groups_used and st.session_state.data_mode == "multi":
            volcano_tag = f"Combined_{volcano_tag}"

        ec1, ec2, ec3 = st.columns(3)
        export_fmt = ec1.selectbox("Format", ["PNG", "PDF", "SVG", "EPS", "JPEG", "TIFF"],
                                    key="volcano_export_fmt")
        export_dpi = ec2.selectbox("Resolution (DPI)", [300, 600, 1200], index=0,
                                    disabled=export_fmt in ("PDF", "SVG", "EPS"), key="volcano_export_dpi")
        mime_map = {"PNG": "image/png", "PDF": "application/pdf", "SVG": "image/svg+xml",
                    "EPS": "application/postscript", "JPEG": "image/jpeg", "TIFF": "image/tiff"}
        ext_map = {"PNG": "png", "PDF": "pdf", "SVG": "svg", "EPS": "eps", "JPEG": "jpg", "TIFF": "tiff"}
        fmt_key = "jpeg" if export_fmt == "JPEG" else export_fmt.lower()
        file_bytes = volcano.export_figure(fig_volc, fmt=fmt_key, dpi=export_dpi)
        volcano_plot_filename = f"{volcano_tag}_VolcanoPlot.{ext_map[export_fmt]}"
        with ec3:
            st.write("")
            st.write("")
            st.download_button(
                f"Download {volcano_plot_filename}", file_bytes,
                file_name=volcano_plot_filename, mime=mime_map[export_fmt]
            )

        sig_col = "p-value" if y_metric == "pvalue" else fdr_label

        st.subheader("Top 20 Biomarkers")
        st.dataframe(utils.format_df_for_display(volcano.top_biomarker_labels(annotated, sig_col=sig_col, n=20)), width='stretch')

        st.subheader("Export Data")
        dc1, dc2, dc3 = st.columns(3)
        up_tbl = volcano.get_direction_table(annotated, "Up")
        down_tbl = volcano.get_direction_table(annotated, "Down")
        sig_tbl = volcano.get_direction_table(annotated, None)
        with dc1:
            up_filename = f"{volcano_tag}_Upregulated_Metabolites.csv"
            st.download_button(f"Download {up_filename}", utils.to_download_bytes_csv(up_tbl),
                                up_filename, "text/csv")
        with dc2:
            down_filename = f"{volcano_tag}_Downregulated_Metabolites.csv"
            st.download_button(f"Download {down_filename}", utils.to_download_bytes_csv(down_tbl),
                                down_filename, "text/csv")
        with dc3:
            complete_filename = f"{volcano_tag}_Complete_Volcano_Data.csv"
            st.download_button(f"Download {complete_filename}", utils.to_download_bytes_csv(annotated),
                                complete_filename, "text/csv")

# ===========================================================================
# TAB 10 — BIOMARKER DISCOVERY
# ===========================================================================
def page_biomarker():
    st.header("Biomarker Discovery")
    if st.session_state.stats_result is None:
        st.warning("Run a two-group statistical comparison in Statistics first.")
    else:
        result = st.session_state.stats_result
        fdr_label = stats_fdr_label(result)
        groups_used = st.session_state.get("stats_result_groups")
        group_col_used = st.session_state.get("stats_result_group_col") or "Group"
        if groups_used:
            st.caption(
                f"Reflects the two {group_col_used} values compared in Statistics: "
                f"**{groups_used[0]}** vs **{groups_used[1]}**"
                + (" (combined normalized data, all datasets)."
                   if st.session_state.data_mode == "multi" else ".")
                + " To compare a different pair or variable, go back to Statistics and "
                "re-run the two-group test."
            )
        c1, c2, c3 = st.columns(3)
        criterion = c1.selectbox("Filtering criterion", ["combined", "pvalue", "fdr"],
                                  format_func=lambda x: {"combined": f"p<0.05 AND {fdr_label}<0.25",
                                                          "pvalue": "p<0.05 only", "fdr": f"{fdr_label}<0.25 only"}[x])
        p_cutoff = c2.number_input("p-value cutoff", value=0.05, min_value=0.0001, max_value=1.0, step=0.01)
        fdr_cutoff = c3.number_input(f"{fdr_label} cutoff", value=0.25, min_value=0.0001, max_value=1.0, step=0.01)

        if st.button("Discover Biomarkers", type="primary", icon="🎯"):
            biomarkers = biomarker.discover_biomarkers(result, criterion=criterion,
                                                        p_cutoff=p_cutoff, fdr_cutoff=fdr_cutoff,
                                                        fdr_label=fdr_label)
            st.session_state.biomarkers = biomarkers
            st.success(f"{len(biomarkers)} biomarker candidates identified.")

        if st.session_state.get("biomarkers") is not None:
            biomarkers = st.session_state.biomarkers
            st.dataframe(utils.format_df_for_display(biomarkers), width='stretch', height=450)
            comparison_name = (f"{groups_used[0]}_vs_{groups_used[1]}" if groups_used else "Comparison")
            biomarker_filename = f"{comparison_name}_Biomarkers.csv"
            st.download_button(
                f"Download {biomarker_filename}", utils.to_download_bytes_csv(biomarkers),
                file_name=biomarker_filename, mime="text/csv", key="dl_biomarkers"
            )

# ===========================================================================
# TAB 11 — HEATMAP
# ===========================================================================
def page_heatmap():
    st.header("Clustered Heatmap")
    if st.session_state.stats_result is None:
        st.warning(
            "Run a two-group comparison in Statistics first — the heatmap needs "
            "that result to know which metabolites are significant."
        )
    else:
        log2_df_full = utils.with_display_feature_names(st.session_state.log2_data)
        result = st.session_state.stats_result
        fdr_label = stats_fdr_label(result)
        meta = st.session_state.meta

        heat_cat_cols = utils.get_categorical_metadata_columns(meta, log2_df_full.columns)
        default_group_col = st.session_state.get("stats_result_group_col") or "Group"
        heat_group_col = st.selectbox(
            "Filter samples by:", heat_cat_cols,
            index=heat_cat_cols.index(default_group_col) if default_group_col in heat_cat_cols else 0,
            key="heatmap_group_col"
        )
        groups_all_heat = meta.loc[meta.index.intersection(log2_df_full.columns), heat_group_col].unique().tolist()
        st.subheader("Select Groups")
        selected_heat_groups = st.multiselect(
            f"{heat_group_col} values to include in the heatmap (choose any subset — 2, 3, 4, or more)",
            groups_all_heat, default=groups_all_heat, key="heatmap_group_select"
        )
        heat_cols = [c for c in log2_df_full.columns if meta.loc[c, heat_group_col] in selected_heat_groups]
        log2_df = log2_df_full[heat_cols]

        st.subheader("Column (Sample) Annotation")
        st.caption(
            "Pick any number of metadata columns to show as stacked annotation tracks above the "
            "heatmap — categorical columns (Diagnosis, Gender, Treatment, ...) render as discrete "
            "color blocks; numeric columns with many distinct values (Age, Body Weight, ...) render "
            "as a color gradient with its own scale."
        )
        col_annot_all_options = [c for c in meta.columns if c != "IsQC"]
        default_col_annot = [heat_group_col] if heat_group_col in col_annot_all_options else []
        col_annot_cols = st.multiselect(
            "Column annotation tracks:", col_annot_all_options, default=default_col_annot,
            key="heatmap_col_annot_cols"
        )

        row_annotations = st.session_state.get("row_annotations")
        row_annot_cols = []
        group_rows_by = None
        if row_annotations is not None:
            st.subheader("Row (Metabolite) Annotation")
            row_annot_cols = st.multiselect(
                "Row annotation tracks:", list(row_annotations.columns),
                default=list(row_annotations.columns)[:1], key="heatmap_row_annot_cols"
            )
            if row_annot_cols:
                n_matched = heatmap_module.count_row_annotation_matches(row_annotations, log2_df_full.index)
                st.caption(f"{n_matched}/{len(log2_df_full.index)} features matched to row annotations "
                           f"by name; unmatched features show as '(unannotated)' in the row bar(s).")

            group_by_options = ["(none — cluster normally)"] + [
                c for c in row_annotations.columns
                if heatmap_module.classify_annotation_series(row_annotations[c]) == "categorical"
            ]
            group_rows_by_choice = st.selectbox(
                "🧩 Group metabolites into module bands by:", group_by_options,
                help="Instead of one global dendrogram, bands metabolites into contiguous, "
                     "labeled blocks by this category (e.g. pathway/class — all Lipids together, "
                     "all Amino acids together, ...) — a 'metabolic fingerprint × clinical feature "
                     "matrix' style view. Rows still cluster within each block if row clustering "
                     "is on below."
            )
            group_rows_by = None if group_rows_by_choice.startswith("(none") else group_rows_by_choice

        # Stage-then-apply color model: editing pickers below never itself
        # triggers a full heatmap re-render (which can be slow with many tracks) --
        # only "Apply Colors" copies the staged picks into what the heatmap below
        # actually uses. "Reset Colors" clears both the applied colors AND each
        # widget's own state, so pickers visually revert to their auto defaults --
        # since widgets already exist earlier in this same run by the time the
        # buttons are clicked, that reversion needs a rerun to take visual effect,
        # which is why reset sets a pending flag + calls st.rerun() rather than
        # popping keys immediately.
        all_annotation_tracks = (
            [(c, meta, "(missing)") for c in col_annot_cols]
            + [(c, row_annotations, "(unannotated)") for c in row_annot_cols]
        )
        if st.session_state.get("heatmap_colors_reset_pending"):
            for track_col, lookup_df, _ in all_annotation_tracks:
                kind = heatmap_module.classify_annotation_series(lookup_df[track_col])
                if kind == "continuous":
                    st.session_state.pop(f"heatmap_annot_cmap_{track_col}", None)
                else:
                    for val in lookup_df[track_col].dropna().unique():
                        st.session_state.pop(f"heatmap_annot_color_{track_col}_{str(val)}", None)
            st.session_state.heatmap_colors_reset_pending = False
            st.session_state.heatmap_annotation_colors_applied = {}

        if all_annotation_tracks:
            with st.expander("🎨 Annotation Colors", expanded=False):
                n_tracks = len(all_annotation_tracks)
                chunk_size = min(6, max(2, n_tracks)) if n_tracks <= 6 else 3
                annotation_colors_live = {}
                for chunk_start in range(0, n_tracks, chunk_size):
                    chunk = all_annotation_tracks[chunk_start:chunk_start + chunk_size]
                    grid_cols = st.columns(chunk_size)
                    for i, (track_col, lookup_df, missing_label) in enumerate(chunk):
                        kind = heatmap_module.classify_annotation_series(lookup_df[track_col])
                        with grid_cols[i]:
                            if kind == "continuous":
                                row = st.columns([1, 1])
                                row[0].caption(track_col)
                                cmap_choice = row[1].selectbox(
                                    track_col, heatmap_module.CONTINUOUS_ANNOT_CMAPS,
                                    key=f"heatmap_annot_cmap_{track_col}", label_visibility="collapsed"
                                )
                                annotation_colors_live[track_col] = cmap_choice
                            else:
                                values = sorted(str(v) for v in lookup_df[track_col].dropna().unique())
                                row = st.columns([1.2] + [0.5] * len(values))
                                row[0].caption(track_col)
                                track_colors = {}
                                for j, val in enumerate(values):
                                    default_hex = heatmap_module.GROUP_PALETTE[j % len(heatmap_module.GROUP_PALETTE)]
                                    track_colors[val] = row[j + 1].color_picker(
                                        val, value=default_hex,
                                        key=f"heatmap_annot_color_{track_col}_{val}",
                                        label_visibility="collapsed"
                                    )
                                annotation_colors_live[track_col] = track_colors

                fc1, fc2, _ = st.columns([1, 1, 4])
                apply_clicked = fc1.button("Apply Colors", type="primary", key="heatmap_apply_colors")
                reset_clicked = fc2.button("Reset Colors", key="heatmap_reset_colors")
                if apply_clicked:
                    st.session_state.heatmap_annotation_colors_applied = annotation_colors_live
                if reset_clicked:
                    st.session_state.heatmap_colors_reset_pending = True
                    st.session_state.heatmap_annotation_colors_applied = {}
                    st.rerun()

        annotation_colors = st.session_state.get("heatmap_annotation_colors_applied") or {}

        CUTOFF_OPTIONS = {
            f"{fdr_label} ≤ 1 (no filter)": (fdr_label, 1.0),
            f"{fdr_label} ≤ 0.25": (fdr_label, 0.25),
            "P-value ≤ 1 (no filter)": ("p-value", 1.0),
            "P-value ≤ 0.05": ("p-value", 0.05),
        }

        st.subheader("Feature Selection")
        c1, c2 = st.columns(2)
        cutoff_label = c1.selectbox("Significant feature cutoff", list(CUTOFF_OPTIONS.keys()), index=1)
        col_preview, cutoff_preview = CUTOFF_OPTIONS[cutoff_label]
        n_preview = (result[col_preview] <= cutoff_preview).sum()
        with c2:
            st.metric("Metabolites at this cutoff", n_preview)
        if n_preview > heatmap_module.HIDE_ROW_LABELS_ABOVE:
            st.warning(
                f"{n_preview} metabolites selected — beyond {heatmap_module.HIDE_ROW_LABELS_ABOVE}, "
                f"row labels become unreadable and will be hidden automatically. Tighten the cutoff "
                f"above for a labeled figure, or proceed for an unlabeled overview heatmap."
            )

        st.subheader("Clustering Options")
        c3, c4, c5 = st.columns(3)
        cluster_mode = c3.selectbox("Clustering", ["Rows only", "Both rows and columns", "Columns only", "No clustering"])
        distance = c4.selectbox("Distance metric", ["euclidean", "correlation"])
        linkage_m = c5.selectbox("Linkage method", ["ward", "average", "complete"])
        cluster_rows = cluster_mode in ("Both rows and columns", "Rows only")
        cluster_cols = cluster_mode in ("Both rows and columns", "Columns only")

        st.subheader("Color Scale Customization")
        c6, c7 = st.columns(2)
        use_custom_gradient = c6.checkbox("Use a custom color chart instead of a preset palette", value=False)
        reverse_cmap = c7.checkbox("Reverse color chart", value=False)

        if use_custom_gradient:
            gc1, gc2, gc3 = st.columns(3)
            color_low = gc1.color_picker("Low color", "#2166AC")
            color_mid = gc2.color_picker("Mid color", "#FFFFFF")
            color_high = gc3.color_picker("High color", "#B2182B")
            custom_colors = [color_low, color_mid, color_high]
            cmap_name = None
        else:
            cmap_name = st.selectbox("Predefined color palette", heatmap_module.PREDEFINED_PALETTES)
            custom_colors = None

        c8, c9 = st.columns(2)
        vmin = c8.number_input("Z-score minimum", value=-2.5, step=0.1)
        vmax = c9.number_input("Z-score maximum", value=2.5, step=0.1)

        use_breakpoints = st.checkbox("Use custom discrete color breakpoints (instead of a continuous scale)", value=False)
        breakpoints = None
        if use_breakpoints:
            bp_text = st.text_input("Comma-separated breakpoints (e.g. -3,-1,0,1,3)",
                                     value=f"{vmin},{vmin/2:.2g},0,{vmax/2:.2g},{vmax}")
            try:
                breakpoints = [float(x.strip()) for x in bp_text.split(",") if x.strip()]
            except ValueError:
                st.warning("Couldn't parse breakpoints — using a continuous scale instead.")
                breakpoints = None

        st.subheader("Font")
        cf1, cf2 = st.columns(2)
        heatmap_font_family = cf1.selectbox("Font family", ["sans-serif", "serif", "monospace"],
                                             key="heatmap_font_family")
        heatmap_font_size = cf2.slider("Base font size", 6, 20, 10, 1, key="heatmap_font_size")
        st.caption("Title, tick labels, legends, and annotation track labels all scale from this "
                   "one size.")

        st.subheader("Figure Size")
        c_size1, c_size2, c_size3 = st.columns(3)
        use_custom_size = c_size1.checkbox("Customize heatmap panel size", value=False)
        heatmap_width_in = heatmap_height_in = None
        if use_custom_size:
            heatmap_width_in = c_size2.number_input("Width (inches)", value=8.0, min_value=2.0, max_value=30.0, step=0.5)
            heatmap_height_in = c_size3.number_input("Height (inches)", value=6.0, min_value=2.0, max_value=30.0, step=0.5)
            st.caption("All other elements (dendrograms, annotation bar, legend, colorbar, labels) "
                       "scale automatically with the size you set here.")

        if len(heat_cols) < 2:
            st.warning("Select at least one group with 2+ samples to build a heatmap.")
        else:
            col, cutoff = CUTOFF_OPTIONS[cutoff_label]
            sig_feats = result[result[col] <= cutoff].index

            if len(sig_feats) < 2:
                st.warning(
                    f"Only {len(sig_feats)} significant feature(s) at this cutoff — "
                    "clustering needs at least 2. Try relaxing the threshold in the Statistics or "
                    "Biomarker Discovery."
                )
            else:
                try:
                    fig_heat, z_ordered, notes = heatmap_module.clustered_heatmap(
                        log2_df, sig_feats, meta=meta, col_annot_cols=col_annot_cols,
                        cluster_rows=cluster_rows, cluster_cols=cluster_cols,
                        distance=distance, linkage_method=linkage_m,
                        cmap_name=cmap_name, custom_colors=custom_colors, reverse_cmap=reverse_cmap,
                        vmin=vmin, vmax=vmax, breakpoints=breakpoints,
                        heatmap_width_in=heatmap_width_in, heatmap_height_in=heatmap_height_in,
                        row_meta=row_annotations, row_annot_cols=row_annot_cols,
                        annotation_colors=annotation_colors, font_family=heatmap_font_family,
                        font_size=heatmap_font_size,
                        title_suffix=f"{len(sig_feats)} features, {cutoff_label}",
                        group_rows_by=group_rows_by,
                    )
                    st.session_state.viz_figs["Clustered Heatmap"] = fig_heat
                    st.session_state.heatmap_fig = fig_heat
                    st.session_state.heatmap_zscore = z_ordered
                    for note in notes:
                        st.info(f"ℹ️ {note}")
                    render_single_figure(fig_heat, width_ratio=(1, 8, 1))
                    st.caption(f"{len(sig_feats)} metabolites shown ({cutoff_label}), "
                               f"{len(heat_cols)} samples across {len(selected_heat_groups)} selected group(s), "
                               f"row-scaled (z-score).")

                    heatmap_tag = "_".join(_safe_tag(g) for g in selected_heat_groups) or "AllGroups"
                    if st.session_state.data_mode == "multi":
                        heatmap_tag = f"Combined_{heatmap_tag}"

                    st.subheader("Export Heatmap")
                    ec1, ec2, ec3 = st.columns(3)
                    export_fmt = ec1.selectbox("Format", ["PNG", "PDF", "SVG", "JPEG", "TIFF"])
                    export_dpi = ec2.selectbox("Resolution (DPI)", [150, 300, 600], index=1)
                    st.caption(
                        "Note: unlike typical vector plots, this heatmap's cells are drawn as an image "
                        "internally, so DPI affects sharpness even for PDF/SVG — use 300+ to avoid a "
                        "blurry/smeared appearance when zoomed in or printed."
                    )
                    mime_map = {"PNG": "image/png", "PDF": "application/pdf", "SVG": "image/svg+xml",
                                "JPEG": "image/jpeg", "TIFF": "image/tiff"}
                    ext_map = {"PNG": "png", "PDF": "pdf", "SVG": "svg", "JPEG": "jpg", "TIFF": "tiff"}
                    fmt_key = "jpeg" if export_fmt == "JPEG" else export_fmt.lower()
                    file_bytes = heatmap_module.export_figure(fig_heat, fmt=fmt_key, dpi=export_dpi)
                    heatmap_filename = f"{heatmap_tag}_Heatmap.{ext_map[export_fmt]}"
                    with ec3:
                        st.write("")
                        st.write("")
                        st.download_button(
                            f"Download {heatmap_filename}", file_bytes,
                            file_name=heatmap_filename, mime=mime_map[export_fmt]
                        )

                    st.subheader("Z-score Data Table")
                    st.caption(
                        "The exact row-scaled Z-score values used to render the heatmap above "
                        "(rows/columns in the same clustered order shown)."
                    )
                    st.dataframe(utils.format_df_for_display(z_ordered), width='stretch', height=300)
                    zscore_filename = f"{heatmap_tag}_Zscore_Table.csv"
                    st.download_button(
                        f"Download {zscore_filename}", utils.to_download_bytes_csv(z_ordered),
                        file_name=zscore_filename, mime="text/csv"
                    )
                except ValueError as e:
                    st.warning(str(e))

# ===========================================================================
# TAB 12 — BOXPLOT OF METABOLITES
# ===========================================================================
def page_boxplot():
    st.header("Boxplot of Metabolites")
    if st.session_state.log2_data is None:
        st.warning(
            "Complete Normalization (strict log2(x) transformation) first."
            if st.session_state.data_mode == "single" else
            "Complete Normalization for every dataset, then click **🔗 Generate Combined "
            "Normalized Data** in Normalization first — that combined table is what "
            "every downstream page analyzes."
        )
    else:
        if st.session_state.data_mode == "multi":
            st.caption("Analyzing the **combined normalized data** across all datasets.")
        fdr_label = current_fdr_label()
        log2_df_full = utils.with_display_feature_names(st.session_state.log2_data)
        meta = st.session_state.meta

        st.subheader("Select Metabolites & Groups")
        selected_metabolites = st.multiselect(
            "Search and select one or more metabolites",
            log2_df_full.index.tolist()
        )
        box_cat_cols = utils.get_categorical_metadata_columns(meta, log2_df_full.columns)
        box_group_col = st.selectbox(
            "Grouping variable:", box_cat_cols,
            index=box_cat_cols.index("Group") if "Group" in box_cat_cols else 0,
            key="boxplot_group_col"
        )
        groups_all_box = meta.loc[meta.index.intersection(log2_df_full.columns), box_group_col].unique().tolist()
        selected_box_groups = st.multiselect(
            f"{box_group_col} values to include in the comparison (choose any subset — 2, 3, 4, or more)",
            groups_all_box, default=groups_all_box, key="boxplot_group_select"
        )

        if selected_metabolites and len(selected_box_groups) >= 2:
            test_name = "Welch's t-test" if len(selected_box_groups) == 2 else "one-way ANOVA"
            st.caption(
                f"Statistics computed on log2-transformed data: {test_name} "
                f"across the {len(selected_box_groups)} selected group(s). {fdr_label} here is corrected "
                f"across only the {len(selected_metabolites)} metabolite(s) shown, not the full feature "
                f"panel — for a panel-wide {fdr_label}, use Statistics."
            )

        st.subheader("Customize Appearance")
        c1, c2, c3 = st.columns(3)
        font_family = c1.selectbox("Font family", ["sans-serif", "serif", "monospace"])
        font_size = c2.slider("Font size", 6, 20, 10, 1)
        ncols = c3.number_input("Panels per row", min_value=1, max_value=6, value=3, step=1)

        c4, c5 = st.columns(2)
        use_custom_box_size = c4.checkbox("Customize figure size", value=False)
        fig_width_in = fig_height_in = None
        if use_custom_box_size:
            fig_width_in = c5.number_input("Width (inches)", value=10.0, min_value=3.0, max_value=30.0, step=0.5)
            fig_height_in = st.number_input("Height (inches)", value=6.0, min_value=3.0, max_value=30.0, step=0.5)

        c6, c7, c8 = st.columns(3)
        show_points = c6.checkbox("Show individual data points", value=True)
        show_mean = c7.checkbox("Show mean (dashed line)", value=False)
        show_median = c8.checkbox("Show median (solid line)", value=True)

        st.caption("Optional: assign a color per group (defaults to a standard palette).")
        group_colors = {}
        if selected_box_groups:
            color_cols = st.columns(min(4, len(selected_box_groups)) or 1)
            for i, g in enumerate(selected_box_groups):
                with color_cols[i % len(color_cols)]:
                    default_c = boxplot_module.GROUP_PALETTE[i % len(boxplot_module.GROUP_PALETTE)]
                    group_colors[g] = st.color_picker(f"Color: {g}", default_c, key=f"boxcolor_{g}")

        if not selected_metabolites:
            st.info("Select at least one metabolite above to generate a boxplot.")
        elif len(selected_box_groups) < 2:
            st.warning("Select at least 2 groups to compare.")
        else:
            box_cols = [c for c in log2_df_full.columns if meta.loc[c, box_group_col] in selected_box_groups]
            log2_df_box = log2_df_full[box_cols]

            stats_table = boxplot_module.compute_stats_for_metabolites(
                log2_df_box, meta, selected_metabolites, selected_box_groups, group_col=box_group_col,
                fdr_label=fdr_label
            )
            fig_box = boxplot_module.boxplot_metabolites(
                log2_df_box, meta, selected_metabolites, selected_box_groups,
                stats_table=stats_table, fig_width_in=fig_width_in, fig_height_in=fig_height_in,
                font_size=font_size, font_family=font_family, group_colors=group_colors,
                show_points=show_points, show_mean=show_mean, show_median=show_median, ncols=ncols,
                group_col=box_group_col, fdr_label=fdr_label,
            )
            st.session_state.viz_figs["Boxplot of Metabolites"] = fig_box
            st.session_state.boxplot_fig = fig_box
            st.session_state.boxplot_stats = stats_table
            render_single_figure(fig_box, width_ratio=(1, 8, 1))

            boxplot_tag = "_".join(_safe_tag(g) for g in selected_box_groups) or "AllGroups"
            if len(selected_metabolites) == 1:
                boxplot_tag = f"{_safe_tag(selected_metabolites[0])}_{boxplot_tag}"
            if st.session_state.data_mode == "multi":
                boxplot_tag = f"Combined_{boxplot_tag}"

            st.subheader("Statistical Results")
            st.dataframe(utils.format_df_for_display(stats_table), width='stretch')
            boxplot_stats_filename = f"{boxplot_tag}_Boxplot_Statistics.csv"
            st.download_button(
                f"Download {boxplot_stats_filename}", utils.to_download_bytes_csv(stats_table),
                file_name=boxplot_stats_filename, mime="text/csv", key="dl_boxplot_stats"
            )

            st.subheader("Export")
            ec1, ec2, ec3 = st.columns(3)
            export_fmt = ec1.selectbox("Format", ["PNG", "PDF", "SVG", "JPEG", "TIFF"], key="boxplot_export_fmt")
            export_dpi = ec2.selectbox("Resolution (DPI)", [150, 300, 600], index=1,
                                        disabled=export_fmt in ("PDF", "SVG"), key="boxplot_export_dpi")
            mime_map = {"PNG": "image/png", "PDF": "application/pdf", "SVG": "image/svg+xml",
                        "JPEG": "image/jpeg", "TIFF": "image/tiff"}
            ext_map = {"PNG": "png", "PDF": "pdf", "SVG": "svg", "JPEG": "jpg", "TIFF": "tiff"}
            fmt_key = "jpeg" if export_fmt == "JPEG" else export_fmt.lower()
            file_bytes = boxplot_module.export_figure(fig_box, fmt=fmt_key, dpi=export_dpi)
            boxplot_filename = f"{boxplot_tag}_Boxplot.{ext_map[export_fmt]}"
            with ec3:
                st.write("")
                st.write("")
                st.download_button(
                    f"Download {boxplot_filename}", file_bytes,
                    file_name=boxplot_filename, mime=mime_map[export_fmt]
                )

# ===========================================================================
# TAB 13 — METABOLITE-SET ENRICHMENT ANALYSIS (MSEA) and
# TAB 14 — METABOLITE–GENE PATHWAY ANALYSIS (MGPA)
# (self-contained module; both read existing Statistics results only — see metabolomics_modules/pathway_analysis.py)
# ===========================================================================
def page_msea():
    pathway_analysis.render_msea_tab(st, st.session_state, render_single_figure)

def page_mgpa():
    pathway_analysis.render_mgpa_tab(st, st.session_state, render_single_figure)


# ===========================================================================
# NETWORK ANALYSIS — CORRELATION REWIRING MAP
# (self-contained: reads the normalized log2 data, writes only rewire_* keys —
#  see metabolomics_modules/rewiring.py for the statistics and metabolomics_modules/rewiring_ui.py for the UI)
# ===========================================================================
def _rewiring_annotations(orig_index, display_index, row_annotations):
    """Row annotations re-indexed onto the display feature names used on this page.
    Each feature is looked up by its full name first, then by the part after the
    multi-dataset 'Dataset::' prefix (the same rule the Heatmap uses), so plain and
    prefixed annotation files both match."""
    if row_annotations is None or row_annotations.empty:
        return None
    cols = {c: rewiring._match_pathways(list(orig_index), row_annotations[c]) for c in row_annotations.columns}
    out = pd.DataFrame(cols, index=pd.Index(display_index, name="Metabolite"))
    out = out.dropna(axis=1, how="all")
    return out if not out.empty else None


def page_rewiring():
    st.header("Correlation Rewiring Map")
    st.markdown(
        "Finds metabolite pairs whose **relationship** changes between two groups, even when their levels "
        "don't. Couplings can be *lost*, *gained* or *flipped*. Pairs are linked to the enzymatic reactions that "
        "connect them (Human-GEM), pathways are tested as sets with the same libraries as MSEA, and the hypothesis "
        "cards check each module for robustness and suggest the experiment that would test it."
    )
    log2_data = st.session_state.log2_data
    if log2_data is None:
        st.info(
            "Complete **Normalization** first; this page analyses the normalized log2 data"
            + (" (in multi-dataset mode, the combined normalized data)." if st.session_state.data_mode == "multi"
               else ".")
            + " You can try the method now on a simulated dataset with known answers."
        )
        if st.toggle("Explore with the simulated rewiring demo (Control vs Tumor, 30 + 30 samples)",
                     key="rewire_use_demo"):
            peak, meta, ann = rewiring_ui.load_simulated_demo()
            st.caption("43 targeted metabolites with planted rewiring (glutamine–TCA decoupling, a tryptophan–"
                       "kynurenine sign flip, a lactate–kynurenine coupling gained in Tumor), two batches, and "
                       "three outlier samples that create a decoy correlation. Your loaded data is not affected.")
            rewiring_ui.render_workspace(st, peak, meta, ann, key_prefix="rewire_demo",
                                         source_label="Simulated demo", data_is_log2=False)
        return

    if st.session_state.data_mode == "multi":
        st.caption("Analyzing the **combined normalized data** across all datasets (already log2).")
    else:
        st.caption("Analyzing the **normalized log2 data** from the Normalization step.")
    display = utils.with_display_feature_names(log2_data)
    ann = _rewiring_annotations(log2_data.index, display.index, st.session_state.get("row_annotations"))
    rewiring_ui.render_workspace(
        st, display, st.session_state.meta, ann, key_prefix="rewire_main",
        source_label=f"Omnix Metabolomics · {st.session_state.data_type}" if st.session_state.data_mode == "single"
        else "Omnix Metabolomics · combined datasets",
        data_is_log2=True,
    )


# ===========================================================================
# Grouped navigation: the SIDEBAR lists group names only; the MAIN PAGE shows
# that group's individual pages as a row of links, and the selected page's
# content renders below the row.
# ===========================================================================
NAV_GROUPS = [
    ("Setup", [
        ("Data Upload", page_data_upload),
    ]),
    ("Preprocessing", [
        ("Data Cleaning & Imputation", page_cleaning),
        ("QC Validation", page_qc),
        ("Combined QC (Multi-Dataset)", page_combined_qc),
    ]),
    ("Normalization", [
        ("Normalization", page_normalization),
        ("Combined Normalization Data (Multi-Dataset)", page_combined_normalized),
    ]),
    ("Statistical Analysis", [
        ("Statistics", page_statistics),
    ]),
    ("Exploratory Analysis", [
        ("PCA", page_pca),
        ("Volcano Plot", page_volcano),
        ("Heatmap", page_heatmap),
        ("Boxplot", page_boxplot),
    ]),
    ("Biomarker Discovery", [
        ("Biomarker Discovery", page_biomarker),
    ]),
    ("Pathway Analysis", [
        ("Metabolite-Set Enrichment Analysis (MSEA)", page_msea),
        ("Metabolite–Gene Pathway Analysis (MGPA)", page_mgpa),
    ]),
    ("Network Analysis", [
        ("Correlation Rewiring Map", page_rewiring),
    ]),
]
_GROUP_BY_NAME = {name: pages for name, pages in NAV_GROUPS}


def _slug(text):
    """CSS-class-safe, unique-enough key fragment for a group/page label."""
    return re.sub(r"[^a-z0-9]+", "-", text.lower()).strip("-")


if "nav_group" not in st.session_state:
    st.session_state.nav_group = NAV_GROUPS[0][0]
if "nav_page" not in st.session_state:
    st.session_state.nav_page = NAV_GROUPS[0][1][0][0]

# --- Sidebar: group names only ---
st.sidebar.title(APP_NAME)
st.sidebar.caption(f"LC-MS based {APP_FAMILY} Analysis (Untargeted & Targeted)" if APP_FAMILY
                   else "LC-MS based Metabolomics & Lipidomics Analysis")
with st.sidebar.container(key="nav_group_box"):
    for group_name, _pages in NAV_GROUPS:
        if st.button(group_name, key=f"navgrp__{_slug(group_name)}", use_container_width=True):
            st.session_state.nav_group = group_name
            st.session_state.nav_page = _pages[0][0]
            st.rerun()

# Highlight the active group: solid filled pill, teal (the workflow-step color) —
# matching MetaboAnalyst's numbered left-hand step list.
st.markdown(
    f"""<style>
    [class*="st-key-nav_group_box"] [class*="st-key-navgrp__{_slug(st.session_state.nav_group)}"] [data-testid="stButton"] button,
    [class*="st-key-nav_group_box"] [class*="st-key-navgrp__{_slug(st.session_state.nav_group)}"] [data-testid="stButton"] button:focus:not(:active) {{
        background: #00695c !important;
        border-left-color: #e65100 !important;
    }}
    [class*="st-key-nav_group_box"] [class*="st-key-navgrp__{_slug(st.session_state.nav_group)}"] [data-testid="stButton"] button p,
    [class*="st-key-nav_group_box"] [class*="st-key-navgrp__{_slug(st.session_state.nav_group)}"] [data-testid="stButton"] button:focus:not(:active) p {{
        color: #ffffff !important;
    }}
    [class*="st-key-nav_group_box"] [class*="st-key-navgrp__{_slug(st.session_state.nav_group)}"] [data-testid="stButton"] button:hover {{
        background: #004d40 !important;
    }}
    [class*="st-key-nav_group_box"] [class*="st-key-navgrp__{_slug(st.session_state.nav_group)}"] [data-testid="stButton"] button:hover p {{
        color: #ffffff !important;
    }}
    </style>""",
    unsafe_allow_html=True,
)

# --- Main page: the active group's pages, as a row of links ---
active_pages = _GROUP_BY_NAME[st.session_state.nav_group]
if st.session_state.nav_page not in dict(active_pages):
    st.session_state.nav_page = active_pages[0][0]

with st.container(key="nav_page_row"):
    cols = st.columns(len(active_pages))
    for col, (page_title, _page_fn) in zip(cols, active_pages):
        with col:
            if st.button(page_title, key=f"navpg__{_slug(st.session_state.nav_group)}__{_slug(page_title)}", use_container_width=False):
                st.session_state.nav_page = page_title
                st.rerun()

# Sub-step (page) tabs get the second MetaboAnalyst accent — orange — so the two nav
# levels (workflow step vs. page within it) read as visually distinct tiers.
st.markdown(
    f"""<style>
    [class*="st-key-nav_page_row"] [class*="st-key-navpg__{_slug(st.session_state.nav_group)}__{_slug(st.session_state.nav_page)}"] [data-testid="stButton"] button,
    [class*="st-key-nav_page_row"] [class*="st-key-navpg__{_slug(st.session_state.nav_group)}__{_slug(st.session_state.nav_page)}"] [data-testid="stButton"] button:focus:not(:active) {{
        background: #e65100 !important;
        border-color: #e65100 !important;
    }}
    [class*="st-key-nav_page_row"] [class*="st-key-navpg__{_slug(st.session_state.nav_group)}__{_slug(st.session_state.nav_page)}"] [data-testid="stButton"] button p,
    [class*="st-key-nav_page_row"] [class*="st-key-navpg__{_slug(st.session_state.nav_group)}__{_slug(st.session_state.nav_page)}"] [data-testid="stButton"] button:focus:not(:active) p {{
        color: #ffffff !important;
        font-weight: 700 !important;
    }}
    [class*="st-key-nav_page_row"] [class*="st-key-navpg__{_slug(st.session_state.nav_group)}__{_slug(st.session_state.nav_page)}"] [data-testid="stButton"] button:hover {{
        background: #bf360c !important;
        border-color: #bf360c !important;
    }}
    [class*="st-key-nav_page_row"] [class*="st-key-navpg__{_slug(st.session_state.nav_group)}__{_slug(st.session_state.nav_page)}"] [data-testid="stButton"] button:hover p {{
        color: #ffffff !important;
    }}
    </style>""",
    unsafe_allow_html=True,
)

st.divider()

dict(active_pages)[st.session_state.nav_page]()
