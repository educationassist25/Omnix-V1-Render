"""
Omnix Transcriptomics — Bulk RNA-seq Differential Expression & Reporting Platform
Streamlit application entry point.

A single-dataset workflow (one expression matrix + one clinical metadata
table) taken from raw counts (or pre-normalized logCPM) straight through
cleaning, QC, DESeq2 normalization/statistics, and biological interpretation.
"""

import os
import re
import sys
import numpy as np
import pandas as pd
import streamlit as st

sys.path.insert(0, os.path.dirname(__file__))
from transcriptomics_modules import (
    utils, qc, rnaseq_io, filtering_rnaseq, normalization_deseq2,
    deseq2_stats, stats_analysis, pca_module, volcano, biomarker, heatmap_module,
    boxplot_module, enrichment,
)

# Inside the Omnix portal (omnix_app.py) the portal sets the page title and header; run on its own,
# the app sets its own.
OMNIX = globals().get("OMNIX", False)
if not OMNIX:
    st.set_page_config(page_title="Omnix · Transcriptomics", layout="wide", page_icon="🧬")

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

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
# Bundled demo gene-set library. No longer offered as a built-in GSEA option, but kept
# in the package so it can still be selected through "Upload custom .gmt file" when
# working offline or with the Rich Clinical Demo dataset.
BUILTIN_GMT_PATH = os.path.join(BASE_DIR, "gene_sets_demo.gmt")

# ---------------------------------------------------------------------------
# Session state initialization
# ---------------------------------------------------------------------------
for key, default in [
    ("input_mode", "Raw Counts"), ("raw_expr_df", None), ("meta", None), ("row_annotations", None),
    ("qc_cols", []), ("sample_cols", []),
    ("filtered_counts", None), ("filter_report", None),
    ("cleaned_logcpm", None),
    ("qc_outliers", None), ("pca_before_after", None),
    ("stats_transformed_result", None), ("stats_multi_deseq_contrasts", None),
    ("log2_data", None), ("dds", None), ("dds_design_factors", None), ("size_factors", None),
    ("norm_counts", None),
    ("stats_result", None), ("stats_group_col", None), ("stats_contrast_label", None),
    ("stats_all_contrasts", None), ("stats_result_groups", None), ("stats_result_group_col", None),
    ("anova_result", None), ("posthoc_result", None), ("anova_groups_used", None),
    ("volcano_settings", None), ("gsea_method_used", None), ("gsea_library_meta", None),
    ("volcano_fig", None), ("volcano_annotated", None),
    ("boxplot_fig", None), ("boxplot_stats", None),
    ("biomarkers", None), ("biomarker_criterion_label", None),
    ("heatmap_fig", None),
    ("gsea_result", None), ("gsea_gene_sets", None), ("gsea_ranked_scores", None),
    ("gsea_run_meta", None), ("processing_notes", []),
]:
    if key not in st.session_state:
        st.session_state[key] = default


def reset_downstream_state():
    for k in ("filtered_counts", "filter_report", "cleaned_logcpm",
              "log2_data", "dds", "dds_design_factors",
              "size_factors", "norm_counts", "qc_outliers", "pca_before_after",
              "stats_result", "stats_group_col", "stats_transformed_result",
              "stats_multi_deseq_contrasts",
              "stats_contrast_label", "stats_all_contrasts", "stats_result_groups",
              "stats_result_group_col", "anova_result", "posthoc_result", "anova_groups_used",
              "volcano_fig", "volcano_annotated", "volcano_settings",
              "boxplot_fig", "boxplot_stats", "biomarkers", "biomarker_criterion_label", "heatmap_fig",
              "gsea_result", "gsea_gene_sets", "gsea_ranked_scores", "gsea_run_meta"):
        st.session_state[k] = None if k != "processing_notes" else []


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
# Navigation (same layout as MetaboAI Pro): workflow steps in the sidebar, the selected
# step's pages as a row of buttons on the main page. Each page's block below runs only
# when that page is open (the number is the page's position in the original tab order).
# ---------------------------------------------------------------------------
NAV_GROUPS = [
    ('Setup', [
        ('Data Upload', 0),
    ]),
    ('Preprocessing', [
        ('Data Cleaning', 1),
        ('QC Validation', 2),
    ]),
    ('Normalization', [
        ('Normalization', 3),
    ]),
    ('Statistical Analysis', [
        ('Statistics', 5),
    ]),
    ('Exploratory Analysis', [
        ('PCA', 4),
        ('Volcano Plot', 6),
        ('Heatmap', 8),
        ('Boxplot', 9),
    ]),
    ('Biomarker Discovery', [
        ('Biomarker Discovery', 7),
    ]),
    ('Pathway Analysis', [
        ('GSEA', 10),
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
st.sidebar.title("Transcriptomics")
st.sidebar.caption("RNA-seq Differential Expression Analysis")
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
    for col, (page_title, _page_idx) in zip(cols, active_pages):
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

ACTIVE_PAGE = dict(active_pages)[st.session_state.nav_page]

# ===========================================================================
# TAB 1 — DATA UPLOAD
# ===========================================================================
if ACTIVE_PAGE == 0:
    st.header("Data Upload & Study Configuration")
    st.caption(
        "A bulk RNA-seq study here is a single dataset — one gene expression matrix "
        "plus one clinical metadata table — cleaned, QC'd, normalized, and analyzed "
        "straight through."
    )

    st.session_state.input_mode = st.selectbox(
        "Input data type",
        ["Raw Counts", "Pre-normalized logCPM"],
        help="**Raw Counts**: integer read counts per gene per sample (STAR+featureCounts, "
             "Salmon+tximport, RSEM, etc.). This is the ONLY mode that supports pyDESeq2's "
             "negative-binomial Wald test in the Statistics page — DESeq2 fits its own size "
             "factors and dispersions from raw counts; feeding it already-normalized values "
             "violates the model. **Pre-normalized logCPM**: already log2(CPM+1)-style values "
             "(e.g. exported from edgeR/limma-voom). No raw counts to model, so Statistics "
             "falls back to classical Welch's t-test / one-way ANOVA on the supplied values."
    )
    is_raw = st.session_state.input_mode == "Raw Counts"

    use_demo = st.checkbox(
        "Use built-in Rich Clinical RNA-seq Demo dataset",
        help="2,178 genes x 46 samples: 4 diagnosis groups (Healthy Control, Prediabetic, "
             "Type 2 Diabetes, Metabolic Syndrome, 10 each) + 6 pooled QC replicates. "
             "Metadata includes Diagnosis, Age, Gender, Treatment, Ethnicity, Body Weight, "
             "and Batch. 300 genes carry real gene symbols across 15 curated pathways "
             "(insulin signaling, inflammation, OXPHOS, adipogenesis, etc.) with genuine "
             "diagnosis-driven differential expression, so Statistics/Volcano/Heatmap/GSEA "
             "all return meaningful, reproducible results out of the box."
    )

    expr_file = meta_file = annot_file = None
    if not use_demo:
        c1, c2 = st.columns(2)
        expr_file = c1.file_uploader(
            "Expression matrix (CSV/XLSX) — genes as rows, samples as columns",
            type=["csv", "xlsx", "xls"], key="expr_upload"
        )
        meta_file = c2.file_uploader(
            "Clinical metadata (CSV/XLSX) — Sample, Group, IsQC, Batch, + any clinical columns",
            type=["csv", "xlsx", "xls"], key="meta_upload"
        )
        annot_file = st.file_uploader(
            "Optional: gene row-annotation file (Gene, Pathway, ...) — unlocks Heatmap row tracks",
            type=["csv", "xlsx", "xls"], key="annot_upload"
        )

    load_clicked = st.button("Load Dataset", type="primary", key="load_dataset_btn")

    if load_clicked:
        try:
            if use_demo:
                expr_raw = pd.read_csv(os.path.join(BASE_DIR, "sample_counts_matrix_richdemo.csv"))
                meta_raw = pd.read_csv(os.path.join(BASE_DIR, "sample_metadata_richdemo.csv"))
                annot_raw = pd.read_csv(os.path.join(BASE_DIR, "gene_row_annotations_richdemo.csv"))
                st.session_state.input_mode = "Raw Counts"
                is_raw = True
            else:
                if expr_file is None or meta_file is None:
                    st.error("Please upload both an expression matrix and a metadata file (or check the demo box).")
                    st.stop()
                expr_raw = utils.load_table(expr_file)
                meta_raw = utils.load_table(meta_file)
                annot_raw = utils.load_table(annot_file) if annot_file is not None else None

            expr_df = rnaseq_io.validate_expression_matrix(expr_raw, mode=st.session_state.input_mode)
            meta = rnaseq_io.validate_clinical_metadata(meta_raw, expr_df.columns)
            qc_cols, sample_cols = rnaseq_io.split_qc_and_samples(expr_df, meta)

            reset_downstream_state()
            st.session_state.raw_expr_df = expr_df
            st.session_state.meta = meta
            st.session_state.qc_cols = qc_cols
            st.session_state.sample_cols = sample_cols
            st.session_state.row_annotations = (
                annot_raw.set_index(annot_raw.columns[0]) if annot_raw is not None else None
            )
            st.success(
                f"Loaded {expr_df.shape[0]} genes x {expr_df.shape[1]} samples "
                f"({len(sample_cols)} biological, {len(qc_cols)} QC)."
            )
        except rnaseq_io.DataValidationError as e:
            st.error(str(e))

    if st.session_state.raw_expr_df is not None:
        st.divider()
        st.subheader("Loaded Data Preview")
        c1, c2 = st.columns(2)
        c1.markdown("**Expression matrix** (first 10 genes)")
        c1.dataframe(st.session_state.raw_expr_df.head(10), width='stretch')
        c2.markdown("**Clinical metadata**")
        c2.dataframe(st.session_state.meta, width='stretch', height=350)

        cat_cols = utils.get_categorical_metadata_columns(
            st.session_state.meta, st.session_state.sample_cols
        )
        st.caption(
            f"Detected categorical/grouping variables: {', '.join(cat_cols) if cat_cols else '(none)'}. "
            "Any of these can drive grouping in PCA, Statistics, Heatmap, and Boxplot."
        )
        st.download_button(
            "Download Raw_Expression_Matrix.csv",
            utils.to_download_bytes_csv(st.session_state.raw_expr_df),
            file_name="Raw_Expression_Matrix.csv", mime="text/csv", key="dl_raw_matrix_tab1", help="As loaded, before any cleaning."
        )

# ===========================================================================
# TAB 2 — DATA CLEANING
# ===========================================================================
if ACTIVE_PAGE == 1:
    st.header("Data Cleaning")
    if st.session_state.raw_expr_df is None:
        st.info("Load a dataset in the Data Upload page first.")
    else:
        bio_cols = st.session_state.sample_cols
        bio_df = st.session_state.raw_expr_df[bio_cols]

        if st.session_state.input_mode == "Raw Counts":
            st.caption(
                "Standard pyDESeq2/DESeq2 preprocessing — no imputation step: bulk "
                "RNA-seq raw counts are essentially never 'missing' (an unmapped gene "
                "is a true zero, not missing data), so cleaning here is the standard "
                "low-expression gene filter used before any DESeq2/edgeR/limma "
                "analysis. Any stray non-numeric cell in the uploaded matrix is "
                "treated as zero (undetected), matching that convention."
            )
            n_missing_cells = int(bio_df.isna().values.sum())
            filter_input = bio_df.fillna(0)
            if n_missing_cells:
                st.caption(f"({n_missing_cells} non-numeric cell(s) treated as zero.)")

            c1, c2 = st.columns(2)
            min_cpm = c1.slider("Minimum CPM", 0.1, 10.0, 1.0, 0.1, key="min_cpm_slider")
            min_frac = c2.slider("Minimum fraction of samples meeting CPM threshold", 0.05, 1.0, 0.2, 0.05,
                                  key="min_frac_slider")

            if st.button("Run Gene Filtering", type="primary", key="run_filter_btn"):
                filtered, report = filtering_rnaseq.filter_low_expression_genes(
                    filter_input, min_cpm=min_cpm, min_fraction_samples=min_frac
                )
                st.session_state.filtered_counts = filtered
                st.session_state.filter_report = report
                st.session_state.processing_notes.append(
                    f"Filtered {report.shape[0] - filtered.shape[0]} of {report.shape[0]} genes "
                    f"below {min_cpm} CPM in >={min_frac*100:.0f}% of samples; {filtered.shape[0]} genes retained."
                )

            if st.session_state.filter_report is not None:
                st.success(
                    f"{st.session_state.filtered_counts.shape[0]} of "
                    f"{st.session_state.filter_report.shape[0]} genes retained."
                )
                fig = filtering_rnaseq.filtering_summary_plot(st.session_state.filter_report, min_cpm)
                utils.render_figure(st, fig)
                utils.render_figure_download(st, fig, "Gene_Filtering_Summary", key_prefix="filter_summary")
                st.dataframe(utils.format_for_display(st.session_state.filter_report), width='stretch', height=300)
                dlc1, dlc2 = st.columns(2)
                dlc1.download_button(
                    "Download Filter_Report.csv", utils.to_download_bytes_csv(st.session_state.filter_report),
                    file_name="Filter_Report.csv", mime="text/csv", key="dl_filter_report"
                )
                dlc2.download_button(
                    "Download Filtered_Counts_Matrix.csv",
                    utils.to_download_bytes_csv(st.session_state.filtered_counts),
                    file_name="Filtered_Counts_Matrix.csv", mime="text/csv", key="dl_filtered_counts"
                )
        else:
            st.caption(
                "Pre-normalized logCPM input: used as supplied — no imputation is "
                "applied. Any missing values are excluded gene-by-gene, automatically, "
                "wherever a downstream statistic (t-test/ANOVA) can't be computed with them."
            )
            n_missing = int(bio_df.isna().values.sum())
            if n_missing:
                st.warning(f"**{n_missing}** missing value(s) present in the supplied logCPM matrix "
                           f"across {bio_df.shape[0]} genes x {bio_df.shape[1]} samples.")
            else:
                st.success("No missing values in the supplied logCPM matrix.")
            st.session_state.cleaned_logcpm = bio_df
            st.dataframe(bio_df.head(10), width='stretch')
            st.download_button(
                "Download Cleaned_LogCPM.csv", utils.to_download_bytes_csv(bio_df),
                file_name="Cleaned_LogCPM.csv", mime="text/csv", key="dl_cleaned_logcpm"
            )

# ===========================================================================
# TAB 3 — QC VALIDATION
# ===========================================================================
if ACTIVE_PAGE == 2:
    st.header("Quality Control Validation")
    if st.session_state.raw_expr_df is None:
        st.info("Load a dataset in the Data Upload page first.")
    elif st.session_state.input_mode == "Raw Counts" and st.session_state.filtered_counts is None:
        st.info("Run gene filtering in the Data Cleaning page first.")
    elif st.session_state.input_mode == "Pre-normalized logCPM" and st.session_state.cleaned_logcpm is None:
        st.info("Visit the Data Cleaning page first.")
    else:
        bio_cols = st.session_state.sample_cols
        qc_cols = st.session_state.qc_cols

        if st.session_state.input_mode == "Raw Counts":
            counts = st.session_state.filtered_counts
            # Library size uses the FULL (unfiltered) gene set — an accurate read on
            # sequencing depth, not just depth among genes that survived filtering.
            full_cols = bio_cols + qc_cols if qc_cols else bio_cols
            all_counts_full = st.session_state.raw_expr_df[full_cols]
            lib_summary = rnaseq_io.library_size_summary(all_counts_full)
            fig_libsize = filtering_rnaseq.library_size_only_plot(lib_summary)
            fig_genedet = filtering_rnaseq.genes_detected_plot(lib_summary)

            log_for_qc = normalization_deseq2.counts_to_log2cpm(all_counts_full)
        else:
            full_cols = bio_cols + qc_cols if qc_cols else bio_cols
            all_counts_full = st.session_state.raw_expr_df.loc[
                st.session_state.cleaned_logcpm.index, full_cols
            ]
            log_for_qc = all_counts_full
            lib_summary = None
            fig_libsize = None
            fig_genedet = None

        bio_log_data = log_for_qc[bio_cols]

        # -----------------------------------------------------------------
        # Build every QC figure up front, then lay them all out together in
        # one balanced 4-then-4 grid (see requirements: 8 figures -> row 1 =
        # 4, row 2 = 4), with each figure's supporting table/download kept
        # immediately below the grid in the same reading order. Library Size
        # and Genes Detected are rendered as two independent single-axis
        # figures (rather than one combined 2-subplot figure) so every panel
        # in this grid shares the same figsize and gets its own controls.
        # -----------------------------------------------------------------
        qc_panels = []
        if fig_libsize is not None:
            qc_panels.append({"title": "**Library Size per Sample**", "fig": fig_libsize,
                               "download_name": "Library_Size_QC", "download_key": "qc_libsize"})
        if fig_genedet is not None:
            qc_panels.append({"title": "**Genes Detected per Sample**", "fig": fig_genedet,
                               "download_name": "Genes_Detected_QC", "download_key": "qc_genedet"})

        cv_table = None
        fig_corr = corr = None
        if qc_cols:
            qc_log_data = log_for_qc[qc_cols]
            cv_table = qc.calculate_cv(2 ** qc_log_data if st.session_state.input_mode == "Raw Counts" else qc_log_data)
            fig_cv_dist = qc.cv_distribution_plot(cv_table)
            fig_cv_hist = qc.cv_histogram(cv_table)
            fig_corr, corr = qc.sample_correlation_matrix(qc_log_data)
            qc_panels.append({"title": "**QC CV Distribution**", "fig": fig_cv_dist,
                               "download_name": "QC_CV_Distribution", "download_key": "qc_cvdist"})
            qc_panels.append({"title": "**Feature Count by CV Quality**", "fig": fig_cv_hist,
                               "download_name": "QC_CV_Quality_Counts", "download_key": "qc_cvhist"})
            qc_panels.append({"title": "**QC Sample Correlation Matrix**", "fig": fig_corr,
                               "download_name": "QC_Sample_Correlation", "download_key": "qc_corr"})

        fig_dist_qc = qc.expression_distribution_plot(
            bio_log_data,
            title="Per-Sample Expression Distribution (log2 CPM, pre-normalization)"
            if st.session_state.input_mode == "Raw Counts"
            else "Per-Sample Expression Distribution (logCPM, as supplied)"
        )
        qc_panels.append({"title": "**Count / Expression Distribution Across Samples**", "fig": fig_dist_qc,
                           "download_name": "Expression_Distribution", "download_key": "qc_dist"})

        fig_all_corr, all_corr = qc.all_sample_correlation_matrix(bio_log_data)
        qc_panels.append({"title": "**Sample Correlation Analysis (All Biological Samples)**", "fig": fig_all_corr,
                           "download_name": "Sample_Correlation_Matrix", "download_key": "qc_allcorr"})

        fig_gene_hist = qc.gene_expression_histogram(
            all_counts_full, title="Gene Expression Distribution (raw values, log10 scale)"
        )
        qc_panels.append({"title": "**Gene Expression Distribution**", "fig": fig_gene_hist,
                           "download_name": "Gene_Expression_Distribution", "download_key": "qc_genehist"})

        st.subheader(f"QC Figures ({len(qc_panels)} total)")
        # 8 figures is the expected count (Library Size, Genes Detected, QC CV
        # Distribution, Feature Count by CV Quality, QC Sample Correlation,
        # Expression Distribution, Sample Correlation (all samples), Gene
        # Expression Distribution) -> row 1 = 4, row 2 = 4, all sharing the
        # same figsize. If QC replicate samples aren't present, 3 of those
        # figures are skipped and the automatic near-square layout is used
        # instead.
        qc_row_layout = [4, 4] if len(qc_panels) == 8 else None
        utils.render_figure_grid(st, qc_panels, row_layout=qc_row_layout)

        if not qc_cols:
            st.info("No QC replicate samples (IsQC=True) in this dataset — "
                     "QC CV Distribution, Feature Count by CV Quality, and QC "
                     "Sample Correlation Matrix are unavailable.")

        if lib_summary is not None:
            st.markdown("**Library Size & Gene Detection — summary table**")
            st.dataframe(utils.format_for_display(lib_summary), width='stretch', height=250)
            st.download_button(
                "Download Library_Size_Summary.csv", utils.to_download_bytes_csv(lib_summary),
                file_name="Library_Size_Summary.csv", mime="text/csv", key="dl_lib_summary"
            )

        if cv_table is not None:
            st.markdown(f"**QC Replicate Reproducibility ({len(qc_cols)} replicates) — CV table**")
            st.dataframe(utils.format_for_display(cv_table), width='stretch', height=250)
            st.download_button(
                "Download QC_CV_Table.csv", utils.to_download_bytes_csv(cv_table),
                file_name="QC_CV_Table.csv", mime="text/csv", key="dl_qc_cv"
            )
            st.download_button(
                "Download QC_Sample_Correlation_Matrix.csv", utils.to_download_bytes_csv(corr),
                file_name="QC_Sample_Correlation_Matrix.csv", mime="text/csv", key="dl_qc_corr_table"
            )

        st.download_button(
            "Download Sample_Correlation_Matrix.csv", utils.to_download_bytes_csv(all_corr),
            file_name="Sample_Correlation_Matrix.csv", mime="text/csv", key="dl_all_corr_table"
        )

        # -----------------------------------------------------------------
        # Outlier sample detection
        # -----------------------------------------------------------------
        st.subheader("Potential Outlier Samples")
        st.caption(
            "Flags samples with (a) unusually low mean correlation to the rest of the "
            "cohort, or (b) unusually large distance from the PC1/PC2 centroid — "
            "either can indicate a technical failure, mislabeling, or genuine biological "
            "outlier. Inspect flagged samples before proceeding to Normalization; nothing "
            "is removed automatically."
        )
        outlier_table = qc.detect_outlier_samples(bio_log_data)
        st.session_state.qc_outliers = outlier_table
        n_flagged = int(outlier_table["Is_Outlier"].sum())
        if n_flagged:
            st.warning(f"**{n_flagged}** sample(s) flagged as potential outliers.")
        else:
            st.success("No samples flagged as potential outliers.")
        st.dataframe(utils.format_for_display(outlier_table), width='stretch', height=300)
        st.download_button(
            "Download QC_Outlier_Report.csv", utils.to_download_bytes_csv(outlier_table),
            file_name="QC_Outlier_Report.csv", mime="text/csv", key="dl_qc_outliers"
        )

# ===========================================================================
# TAB 4 — NORMALIZATION
# ===========================================================================
if ACTIVE_PAGE == 3:
    st.header("Normalization")
    if st.session_state.raw_expr_df is None:
        st.info("Load a dataset in the Data Upload page first.")
    elif st.session_state.input_mode == "Raw Counts" and st.session_state.filtered_counts is None:
        st.info("Run gene filtering in the Data Cleaning page first.")
    elif st.session_state.input_mode == "Pre-normalized logCPM" and st.session_state.cleaned_logcpm is None:
        st.info("Visit the Data Cleaning page first.")
    else:
        meta = st.session_state.meta
        bio_cols = st.session_state.sample_cols

        if st.session_state.input_mode == "Raw Counts":
            st.caption(
                "pyDESeq2 computes median-of-ratios **size factors** (correcting for "
                "sequencing depth and RNA composition) and fits per-gene dispersions, "
                "then produces a **variance-stabilizing transform (VST)** — the "
                "recommended matrix for PCA, clustering, and heatmaps."
            )
            cat_cols = utils.get_categorical_metadata_columns(meta, bio_cols)
            c1, c2 = st.columns(2)
            primary_col = c1.selectbox("Primary grouping variable (design factor)", cat_cols, key="norm_primary_col")
            covariate_options = [c for c in cat_cols if c != primary_col]
            covariates = c2.multiselect(
                "Additional covariates to control for (optional, e.g. Batch)",
                covariate_options, key="norm_covariates"
            )
            design_factors = covariates + [primary_col]

            if st.button("Run DESeq2 Normalization", type="primary", key="run_deseq2_norm_btn"):
                try:
                    with st.spinner("Fitting DESeq2 (size factors, dispersions, GLM)..."):
                        dds = normalization_deseq2.fit_deseq2(
                            st.session_state.filtered_counts, meta, design_factors
                        )
                    st.session_state.dds = dds
                    st.session_state.dds_design_factors = design_factors
                    st.session_state.size_factors = normalization_deseq2.size_factors_table(dds)
                    st.session_state.norm_counts = normalization_deseq2.normalized_counts(dds)
                    with st.spinner("Computing variance-stabilizing transform..."):
                        st.session_state.log2_data = normalization_deseq2.vst_transform(dds)
                    st.session_state.processing_notes.append(
                        f"DESeq2 fit with design ~ {' + '.join(design_factors)}."
                    )
                    st.success("DESeq2 normalization complete.")
                except normalization_deseq2.NormalizationError as e:
                    st.error(str(e))
        else:
            st.caption(
                "Data was already supplied as pre-normalized logCPM — no further "
                "normalization is applied; this matrix is used directly for PCA, "
                "Statistics (classical t-test/ANOVA), Heatmap, and Boxplot."
            )
            if st.button("Use logCPM Matrix As-Is", type="primary", key="use_logcpm_asis_btn"):
                st.session_state.log2_data = st.session_state.cleaned_logcpm
                st.success("logCPM matrix set as the active expression matrix.")

        if st.session_state.log2_data is not None:
            st.divider()

            # ---------------------------------------------------------
            # Build every Normalization figure up front, then lay them all
            # out together in one balanced grid (5 figures -> row 1 = 3,
            # row 2 = 2): DESeq2 Size Factors, Before/After Distribution
            # (raw + VST, as two separate panels), PCA Before/After.
            # ---------------------------------------------------------
            norm_panels = []
            if st.session_state.size_factors is not None:
                fig_sf = normalization_deseq2.size_factor_plot(st.session_state.size_factors)
                norm_panels.append({"title": "**DESeq2 Size Factors**", "fig": fig_sf,
                                     "download_name": "DESeq2_Size_Factors", "download_key": "norm_sf"})

            if st.session_state.input_mode == "Raw Counts":
                fig_before_dist = normalization_deseq2.raw_counts_distribution_plot(
                    st.session_state.filtered_counts
                )
                fig_after_dist = normalization_deseq2.vst_distribution_plot(st.session_state.log2_data)
                norm_panels.append({"title": "**Before: Raw Counts Distribution**", "fig": fig_before_dist,
                                     "download_name": "Normalization_Distribution_Before",
                                     "download_key": "norm_dist_before"})
                norm_panels.append({"title": "**After: VST-Normalized Distribution**", "fig": fig_after_dist,
                                     "download_name": "Normalization_Distribution_After",
                                     "download_key": "norm_dist_after"})

            before_after_scores = None
            if st.session_state.input_mode == "Raw Counts" and st.session_state.filtered_counts is not None:
                before_log2 = normalization_deseq2.counts_to_log2cpm(st.session_state.filtered_counts)
                after_log2 = st.session_state.log2_data
                common_cols = [c for c in before_log2.columns if c in after_log2.columns]
                if len(common_cols) >= 4:
                    pca_before, scores_before, _ = pca_module.run_pca(before_log2[common_cols], n_components=2)
                    pca_after, scores_after, _ = pca_module.run_pca(after_log2[common_cols], n_components=2)
                    pca_cat_cols_n = utils.get_categorical_metadata_columns(meta, common_cols)
                    group_col_n = "Diagnosis" if "Diagnosis" in pca_cat_cols_n else (
                        "Group" if "Group" in pca_cat_cols_n else (pca_cat_cols_n[0] if pca_cat_cols_n else None))
                    if group_col_n:
                        fig_before = pca_module.pca_score_plot(pca_before, scores_before, meta,
                                                                 show_ellipse=False, group_col=group_col_n)
                        fig_before.axes[0].set_title(f"Before Normalization (colored by {group_col_n})")
                        fig_after = pca_module.pca_score_plot(pca_after, scores_after, meta,
                                                                show_ellipse=False, group_col=group_col_n)
                        fig_after.axes[0].set_title(f"After Normalization (colored by {group_col_n})")
                        norm_panels.append({"title": "**PCA — Before Normalization**", "fig": fig_before,
                                             "download_name": "PCA_Before_Normalization",
                                             "download_key": "pca_before"})
                        norm_panels.append({"title": "**PCA — After Normalization**", "fig": fig_after,
                                             "download_name": "PCA_After_Normalization",
                                             "download_key": "pca_after"})
                        before_after_scores = scores_before.add_prefix("Before_").join(
                            scores_after.add_prefix("After_")
                        )
                        st.session_state.pca_before_after = before_after_scores

            st.subheader(f"Normalization Figures ({len(norm_panels)} total)")
            if st.session_state.input_mode == "Raw Counts":
                st.caption(
                    "PCA (Before/After) panels: Before = PCA on simple log2(CPM+1) of the filtered "
                    "raw counts (no size-factor correction); After = PCA on the DESeq2 "
                    "variance-stabilized (VST) matrix. Samples that cluster by technical "
                    "batch/library-size before, but by biological group after, confirm "
                    "normalization worked as intended."
                )
            # 5 figures is the expected count -> row 1 = 3, row 2 = 2. If fewer
            # figures are available (e.g. logCPM input mode, or <4 samples for
            # the PCA comparison), fall back to the automatic near-square layout.
            norm_row_layout = [3, 2] if len(norm_panels) == 5 else None
            utils.render_figure_grid(st, norm_panels, row_layout=norm_row_layout)

            if st.session_state.size_factors is not None:
                st.markdown("**DESeq2 Size Factors — table**")
                st.dataframe(utils.format_for_display(st.session_state.size_factors), width='stretch', height=250)

            dlA, dlB = st.columns(2)
            with dlA:
                st.download_button(
                    "Download Normalized_Counts.csv",
                    utils.to_download_bytes_csv(st.session_state.norm_counts)
                    if st.session_state.norm_counts is not None
                    else utils.to_download_bytes_csv(st.session_state.log2_data),
                    file_name="Normalized_Counts.csv", mime="text/csv", key="dl_norm_counts_linear", help="Linear scale, DESeq2 size-factor-corrected.",
                    disabled=st.session_state.norm_counts is None,
                )
            with dlB:
                st.download_button(
                    "Download Normalized_Expression_Matrix.csv",
                    utils.to_download_bytes_csv(st.session_state.log2_data),
                    file_name="Normalized_Expression_Matrix.csv", mime="text/csv", key="dl_norm_matrix", help="VST / log2-transformed, for plotting."
                )

            if before_after_scores is not None:
                st.download_button(
                    "Download PCA_Before_After_Scores.csv",
                    utils.to_download_bytes_csv(before_after_scores),
                    file_name="PCA_Before_After_Scores.csv", mime="text/csv",
                    key="dl_pca_before_after"
                )
            elif st.session_state.input_mode == "Raw Counts" and st.session_state.filtered_counts is not None:
                st.info("Need at least 4 samples to compare PCA before/after normalization.")

# ===========================================================================
# TAB 5 — PCA
# ===========================================================================
if ACTIVE_PAGE == 4:
    st.header("Principal Component Analysis")
    if st.session_state.log2_data is None:
        st.info("Run Normalization first.")
    else:
        log2_df_full = st.session_state.log2_data
        meta = st.session_state.meta
        pca_cat_cols = utils.get_categorical_metadata_columns(meta, log2_df_full.columns)
        pca_group_col = st.selectbox(
            "Color/group samples by:", pca_cat_cols,
            index=pca_cat_cols.index("Diagnosis") if "Diagnosis" in pca_cat_cols
            else (pca_cat_cols.index("Group") if "Group" in pca_cat_cols else 0),
            key="pca_group_col"
        )
        groups_all_pca = meta.loc[meta.index.intersection(log2_df_full.columns), pca_group_col].dropna().unique().tolist()

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
            n_comp = st.slider("Number of components", 2,
                                min(10, log2_df.shape[1] - 1 if log2_df.shape[1] > 2 else 2), 5)

            st.subheader("Customize Appearance")
            c1, c2 = st.columns(2)
            color_mode = c1.radio("Color mode", ["Preset palette", "Custom colors (pick each group)"],
                                   key="pca_color_mode")
            show_ellipse = c2.checkbox("Show 95% confidence ellipses", value=True, key="pca_ellipse")

            if color_mode == "Preset palette":
                palette_choice = st.selectbox("Color palette", list(pca_module.PALETTES.keys()), key="pca_palette")
            else:
                st.caption("Pick an exact color for each group with the color picker below.")
                default_swatches = pca_module.PALETTES["Default (tab10)"]
                palette_choice = {}
                pick_cols = st.columns(min(4, len(selected_pca_groups)) or 1)
                for i, g in enumerate(selected_pca_groups):
                    with pick_cols[i % len(pick_cols)]:
                        palette_choice[g] = st.color_picker(
                            f"Color: {g}", default_swatches[i % len(default_swatches)], key=f"pca_color_{g}"
                        )

            st.caption("Optional: assign a marker style per group (defaults to circles for all).")
            marker_map = {}
            marker_cols = st.columns(min(4, len(selected_pca_groups)) or 1)
            for i, g in enumerate(selected_pca_groups):
                with marker_cols[i % len(marker_cols)]:
                    marker_map[g] = st.selectbox(f"Marker: {g}", pca_module.MARKER_STYLES, key=f"pca_marker_{g}")

            pca_obj, scores_df, cols = pca_module.run_pca(log2_df, n_components=n_comp)
            fig_score = pca_module.pca_score_plot(pca_obj, scores_df, meta, palette=palette_choice,
                                                   marker_map=marker_map, show_ellipse=show_ellipse,
                                                   group_col=pca_group_col)
            fig_loading, top_loadings = pca_module.pca_loading_plot(pca_obj, log2_df.index, top_n=20)
            fig_var = pca_module.pca_variance_plot(pca_obj)

            utils.render_figure_grid(st, [
                {"fig": fig_score, "download_name": "PCA_Score_Plot", "download_key": "pca_score"},
                {"fig": fig_loading, "download_name": "PCA_Loading_Plot", "download_key": "pca_loading"},
                {"fig": fig_var, "download_name": "PCA_Variance_Plot", "download_key": "pca_var"},
            ])

            st.subheader("Top Contributing Genes (Loadings)")
            st.dataframe(utils.format_for_display(top_loadings), width='stretch')

            st.divider()
            st.subheader("Download PCA Data")
            st.caption(
                "The underlying data behind the plots above: sample scores per principal "
                "component (with group assignment), % variance explained per component, "
                "the gene loadings, and the exact normalized expression matrix PCA was run on."
            )
            explained_var_df = pd.DataFrame({
                "Component": [f"PC{i+1}" for i in range(len(pca_obj.explained_variance_ratio_))],
                "Variance Explained (%)": (pca_obj.explained_variance_ratio_ * 100).round(3),
                "Cumulative Variance Explained (%)": (np.cumsum(pca_obj.explained_variance_ratio_) * 100).round(3),
            })
            scores_with_group = scores_df.copy()
            scores_with_group.insert(0, pca_group_col, meta.loc[scores_with_group.index, pca_group_col].values)

            pca_tag = "_".join(_safe_tag(g) for g in selected_pca_groups) or "AllGroups"

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
# TAB 6 — STATISTICS
# ===========================================================================
if ACTIVE_PAGE == 5:
    st.header("Statistics")
    if st.session_state.log2_data is None:
        st.info("Run Normalization first.")
    else:
        log2_df = st.session_state.log2_data
        meta = st.session_state.meta
        sample_cols = log2_df.columns.tolist()
        stats_cat_cols = utils.get_categorical_metadata_columns(meta, sample_cols)
        grouping_var = st.selectbox(
            "Grouping variable:", stats_cat_cols,
            index=stats_cat_cols.index("Diagnosis") if "Diagnosis" in stats_cat_cols
            else (stats_cat_cols.index("Group") if "Group" in stats_cat_cols else 0),
            key="stats_grouping_var"
        )
        groups_available = meta.loc[meta.index.intersection(sample_cols), grouping_var].dropna().unique().tolist()

        st.caption(
            "Two groups → Wald test (pyDESeq2's negative-binomial GLM for raw counts, or a "
            "linear-model Wald test for logCPM — same statistic and output schema either way); "
            "three or more groups → one-way ANOVA. Choose the comparison type below."
        )

        mode = st.radio("Comparison type", ["Two-group comparison", "ANOVA (≥3 groups)"], horizontal=True,
                         key="stats_mode_radio")

        # =================================================================
        # TWO-GROUP COMPARISON
        # =================================================================
        if mode == "Two-group comparison":
            c1, c2, c3 = st.columns(3)
            group_a = c1.selectbox("Group A", groups_available, index=0, key="stats_group_a")
            # Default Group B to the strongest/most clinically distinct contrast against
            # Group A when it's present (e.g. "Type 2 Diabetes" vs a "Healthy Control"
            # Group A on the Rich Clinical Demo), rather than always the 2nd group in
            # discovery order -- so a first-time user's default comparison is a
            # biologically meaningful one with real signal, not the most subtle pair.
            _default_b_candidates = ["Type 2 Diabetes", "Metabolic Syndrome"]
            _default_b_idx = next(
                (groups_available.index(g) for g in _default_b_candidates
                 if g in groups_available and g != group_a),
                min(1, len(groups_available) - 1)
            )
            group_b = c2.selectbox("Group B", groups_available, index=_default_b_idx,
                                    key="stats_group_b")

            st.markdown("**Significance Threshold**")
            sc1, sc2, sc3 = st.columns(3)
            sig_metric = sc1.selectbox(
                "Metric", ["FDR", "P-value"], key="stats_sig_metric",
                help="Which statistic the 'Significant' column and the filter below are based on."
            )
            sig_fdr_threshold = sc2.number_input(
                "Threshold", min_value=0.0001, max_value=1.0, value=0.05, step=0.01,
                key="stats_sig_threshold",
                help="Genes with the selected metric ≤ this value are marked Significant."
            )
            sig_apply_filter = sc3.checkbox(
                "Apply significance filter", value=False, key="stats_sig_apply_filter",
                help="Shows only genes meeting the selected significance criterion. Applies to "
                     "every two-group result table on this page (DESeq2 Wald test, logCPM Wald "
                     "test, and the Welch's t-test) and their CSV downloads. Leave unchecked to see all "
                     "genes."
            )

            if st.session_state.input_mode == "Raw Counts":
                method = c3.selectbox("Method", ["pyDESeq2 Wald test (recommended for raw counts)"],
                                       key="stats_method_raw")
                covariate_options = [c for c in stats_cat_cols if c != grouping_var]
                covariates = st.multiselect("Covariates to control for (optional, e.g. Batch)",
                                             covariate_options, key="stats_covariates")
                design_factors = covariates + [grouping_var]
                st.caption(
                    "pyDESeq2's Wald test on raw, gene-filtered counts is the count-appropriate "
                    "equivalent of a t-test here — DESeq2 models counts with its own "
                    "negative-binomial GLM rather than assuming normally-distributed values."
                )

                if st.button("Run Two-Group Test", key="run_two_group_raw") and group_a != group_b:
                    try:
                        with st.spinner(f"Fitting DESeq2 for {grouping_var} and testing {group_a} vs {group_b}..."):
                            # reference_level=group_b sets group_b as the model's baseline
                            # category, which is required for apeGLM LFC shrinkage (the
                            # coefficient for group_a must exist as "<factor>[T.<group_a>]").
                            dds = normalization_deseq2.fit_deseq2(
                                st.session_state.filtered_counts, meta, design_factors,
                                reference_level=group_b
                            )
                            result = deseq2_stats.deseq2_contrast(dds, grouping_var, group_a, group_b,
                                                                   shrink_lfc=True,
                                                                   fdr_threshold=sig_fdr_threshold)
                            result = utils.apply_significance_filter(
                                result, sig_metric, sig_fdr_threshold, apply_filter=False)
                        st.session_state.stats_result = result
                        st.session_state.stats_result_groups = (group_a, group_b)
                        st.session_state.stats_result_group_col = grouping_var
                        shrunk_note = " (LFC shrunk via apeGLM prior)" if result.attrs.get("lfc_shrunk") else ""
                        st.session_state.processing_notes.append(
                            f"DESeq2 Wald test ({grouping_var}): {group_a} vs {group_b} "
                            f"(design ~ {' + '.join(design_factors)}){shrunk_note}; "
                            f"{int(result['Significant'].sum())} significant genes ({sig_metric} ≤ {sig_fdr_threshold:g})."
                        )
                        st.success(f"Test complete: {int(result['Significant'].sum())} significant genes found.")
                        if result.attrs.get("lfc_shrunk"):
                            st.caption(
                                "Log2FC/Linear_FC above are pyDESeq2's **apeGLM-shrunk** estimates "
                                "(its own recommended value for ranking/plotting) — the original, "
                                "unshrunk MLE estimate is kept in Log2FC_MLE_Unshrunk / "
                                "Linear_FC_MLE_Unshrunk for reference. Shrinkage does not change the "
                                "p-value or FDR."
                            )
                    except deseq2_stats.StatsError as e:
                        st.error(str(e))

                run_welch_too = st.checkbox(
                    "Also compute a Welch's two-sample t-test on the transformed (VST/log2) expression "
                    "data for this same comparison — a descriptive statistic kept separate from the "
                    "DESeq2 model-based result above.",
                    value=False, key="stats_run_welch_on_raw"
                )
                if run_welch_too and st.button("Run Welch t-test (transformed data)", key="run_welch_raw_btn") \
                        and group_a != group_b and st.session_state.log2_data is not None:
                    a_samples = [s for s in log2_df.columns if meta.loc[s, grouping_var] == group_a]
                    b_samples = [s for s in log2_df.columns if meta.loc[s, grouping_var] == group_b]
                    welch_result = deseq2_stats.classical_two_group(
                        log2_df, a_samples, b_samples, method="ttest",
                        group_a_label=group_a, group_b_label=group_b,
                        fdr_threshold=sig_fdr_threshold
                    )
                    welch_result = utils.apply_significance_filter(
                        welch_result, sig_metric, sig_fdr_threshold, apply_filter=False)
                    st.session_state.stats_transformed_result = welch_result
                    st.session_state.processing_notes.append(
                        f"Welch's t-test on VST ({grouping_var}): {group_a} vs {group_b}; "
                        f"{int(welch_result['Significant'].sum())} significant genes ({sig_metric} ≤ {sig_fdr_threshold:g})."
                    )
            else:
                method = c3.selectbox("Method", ["Wald test (linear model on logCPM)"],
                                       key="stats_method_logcpm")
                covariate_options = [c for c in stats_cat_cols if c != grouping_var]
                covariates = st.multiselect("Covariates to control for (optional, e.g. Batch)",
                                             covariate_options, key="stats_covariates_logcpm")
                st.caption(
                    "Both input modes use a Wald test (coefficient ÷ standard error, referenced to a "
                    "normal distribution) and return the identical result schema. They differ only in "
                    "the model the input allows: raw counts get pyDESeq2's negative-binomial GLM with "
                    "dispersion shrinkage, while logCPM — no longer counts — gets a Gaussian linear "
                    "model, as limma/voom-style workflows do with already-normalized values."
                )

                if st.button("Run Two-Group Test", key="run_two_group_logcpm") and group_a != group_b:
                    try:
                        result = deseq2_stats.wald_test_logcpm(
                            log2_df, meta, grouping_var, group_a, group_b, covariates=covariates,
                            fdr_threshold=sig_fdr_threshold
                        )
                        result = utils.apply_significance_filter(
                            result, sig_metric, sig_fdr_threshold, apply_filter=False)
                        st.session_state.stats_result = result
                        st.session_state.stats_result_groups = (group_a, group_b)
                        st.session_state.stats_result_group_col = grouping_var
                        design_desc = " + ".join(covariates + [grouping_var])
                        st.session_state.processing_notes.append(
                            f"Wald test on logCPM ({grouping_var}): {group_a} vs {group_b} "
                            f"(design ~ {design_desc}); "
                            f"{int(result['Significant'].sum())} significant genes ({sig_metric} ≤ {sig_fdr_threshold:g})."
                        )
                        st.success(f"Test complete: {int(result['Significant'].sum())} significant genes found.")
                    except deseq2_stats.StatsError as e:
                        st.error(str(e))

            if st.session_state.stats_result is not None and st.session_state.get("stats_result_groups"):
                result = st.session_state.stats_result
                g_a, g_b = st.session_state.stats_result_groups
                st.session_state.stats_contrast_label = f"{g_a} vs {g_b}"
                comparison_name = f"{g_a}_vs_{g_b}"
                # Display/export only the 6 requested columns for the primary two-group
                # test output (Gene, p-value, FDR, Log2FC, Linear_FC, Significant). This
                # is a display-only view built on a copy -- the full `result` (with its
                # fuller set of intermediate columns, e.g. baseMean, lfcSE, stat,
                # Group_A/Group_B, Higher_In, Log2FC_MLE_Unshrunk) is left untouched in
                # session state for Volcano/Biomarker/Heatmap and other downstream tabs.
                # Dynamic Significance Threshold filter (Metric/Threshold/Apply widgets above):
                # recomputes 'Significant' against the currently-selected metric+threshold and,
                # when the checkbox is on, drops every row that doesn't meet it -- applied to
                # both the on-screen table and its CSV download below.
                result_filtered = utils.apply_significance_filter(
                    result, sig_metric, sig_fdr_threshold, sig_apply_filter)
                result_display = utils.two_group_display_view(result_filtered)
                st.dataframe(utils.format_for_display(result_display), width='stretch', height=400)
                st.download_button(
                    f"Download {comparison_name}_Statistics.csv",
                    utils.to_download_bytes_csv(result_display),
                    file_name=f"{comparison_name}_Statistics.csv", mime="text/csv",
                    key="dl_stats_twogroup"
                )
                st.caption(f"Filename includes the comparison ({comparison_name}) so results from "
                           "different comparisons stay distinguishable. Showing "
                           f"{len(result_display)} of {len(result)} genes"
                           f"{' (filtered by ' + sig_metric + ' ≤ ' + format(sig_fdr_threshold, 'g') + ')' if sig_apply_filter else ''}.")

                if st.session_state.get("stats_transformed_result") is not None:
                    st.divider()
                    st.markdown(
                        "**Transformed-data statistics (Welch's t-test on VST/log2 values)** — "
                        "descriptive only; kept separate from the DESeq2 model-based result above "
                        "(different model, different Log2FC/p-value)."
                    )
                    tresult = st.session_state.stats_transformed_result
                    # Same 6-column display/export restriction and Significance Threshold
                    # filter as the primary result above.
                    tresult_filtered = utils.apply_significance_filter(
                        tresult, sig_metric, sig_fdr_threshold, sig_apply_filter)
                    tresult_display = utils.two_group_display_view(tresult_filtered)
                    st.dataframe(utils.format_for_display(tresult_display), width='stretch', height=350)
                    st.download_button(
                        f"Download {comparison_name}_Welch_Transformed_Statistics.csv",
                        utils.to_download_bytes_csv(tresult_display),
                        file_name=f"{comparison_name}_Welch_Transformed_Statistics.csv", mime="text/csv",
                        key="dl_stats_welch_transformed"
                    )

        # =================================================================
        # ANOVA (≥3 groups) — classical one-way ANOVA on the normalized/VST
        # expression matrix (DESeq2 has no direct omnibus-ANOVA equivalent for
        # an arbitrary number of groups without a likelihood-ratio-test design,
        # so both input modes share this well-tested path here).
        # =================================================================
        else:
            multi_mode = "ANOVA on transformed data (VST/logCPM)"
            if st.session_state.input_mode == "Raw Counts":
                multi_mode = st.radio(
                    "Multi-group comparison type",
                    ["ANOVA on transformed data (VST/logCPM)",
                     "PyDESeq2 count-based contrasts (model-based, recommended for raw counts)"],
                    key="multi_group_mode_radio"
                )
                if multi_mode.startswith("PyDESeq2"):
                    st.caption(
                        "Fits pyDESeq2's negative-binomial GLM once across the selected groups (the "
                        "appropriate multi-factor/contrast framework for count data — ordinary ANOVA is "
                        "never applied directly to raw counts here), then draws pairwise Wald contrasts "
                        "from that single fit, each with its own Benjamini-Hochberg FDR."
                    )
                else:
                    st.caption(
                        "Computed on the DESeq2 variance-stabilized (VST) expression matrix from the "
                        "Normalization page — a classical one-way ANOVA on transformed values. For "
                        "count-native multi-group contrasts, choose the PyDESeq2 option above instead."
                    )
            if len(groups_available) < 3:
                st.warning(f"Need at least 3 values of '{grouping_var}' (excluding QC) for a multi-group comparison. "
                           "Pick a different grouping variable, or add more groups in metadata.")
            elif multi_mode.startswith("PyDESeq2"):
                # =============================================================
                # PyDESeq2 count-based multi-group contrast framework
                # =============================================================
                deseq_groups = st.multiselect(
                    f"{grouping_var} values to include (select 3 or more)",
                    groups_available, default=groups_available, key="deseq_multi_group_select"
                )
                contrast_style = st.radio(
                    "Contrasts", ["Every group vs a reference (Dunnett-style)", "Every pairwise combination"],
                    key="deseq_multi_contrast_style"
                )
                reference_level = None
                if contrast_style.startswith("Every group vs"):
                    reference_level = st.selectbox("Reference group", deseq_groups, key="deseq_multi_reference")
                covariate_options_m = [c for c in stats_cat_cols if c != grouping_var]
                covariates_m = st.multiselect("Covariates to control for (optional, e.g. Batch)",
                                               covariate_options_m, key="deseq_multi_covariates")
                design_factors_m = covariates_m + [grouping_var]
                st.markdown("**Significance Threshold**")
                mc1, mc2, mc3 = st.columns(3)
                sig_metric_m = mc1.selectbox(
                    "Metric", ["FDR", "P-value"], key="deseq_multi_sig_metric",
                    help="Which statistic the 'Significant' column and the filter below are based on."
                )
                sig_fdr_threshold_m = mc2.number_input(
                    "Threshold", min_value=0.0001, max_value=1.0, value=0.05, step=0.01,
                    key="deseq_multi_sig_threshold",
                    help="Genes with the selected metric ≤ this value are marked Significant."
                )
                sig_apply_filter_m = mc3.checkbox(
                    "Apply significance filter", value=False, key="deseq_multi_sig_apply_filter",
                    help="Shows only genes meeting the selected significance criterion, in each "
                         "pairwise contrast below and its CSV download."
                )

                if len(deseq_groups) < 3:
                    st.warning(f"Select at least 3 groups (currently {len(deseq_groups)} selected).")
                elif st.button("Run PyDESeq2 Contrasts", key="run_deseq_multi_btn"):
                    try:
                        with st.spinner("Fitting DESeq2 and computing pairwise contrasts..."):
                            sub_meta_mask = meta.loc[meta.index.intersection(sample_cols), grouping_var].isin(deseq_groups)
                            sub_samples = sub_meta_mask.index[sub_meta_mask].tolist()
                            sub_counts = st.session_state.filtered_counts[
                                [c for c in st.session_state.filtered_counts.columns if c in sub_samples]
                            ]
                            is_dunnett = contrast_style.startswith("Every group vs")
                            # reference_level lets apeGLM LFC shrinkage apply to every contrast in
                            # one shared fit when every contrast is vs the same reference (Dunnett-
                            # style). "All pairwise" has no single reference, so shrinkage is skipped
                            # there and the unshrunk MLE result is reported instead.
                            dds_m = normalization_deseq2.fit_deseq2(
                                sub_counts, meta, design_factors_m,
                                reference_level=reference_level if is_dunnett else None
                            )
                            contrasts = deseq2_stats.all_pairwise_contrasts(
                                dds_m, grouping_var, deseq_groups,
                                reference=reference_level if is_dunnett else None,
                                fdr_threshold=sig_fdr_threshold_m
                            )
                        st.session_state.stats_multi_deseq_contrasts = contrasts
                        st.session_state.anova_groups_used = list(deseq_groups)
                        st.session_state.stats_result_group_col = grouping_var
                        st.session_state.processing_notes.append(
                            f"PyDESeq2 multi-group contrasts ({grouping_var}, {len(deseq_groups)} groups, "
                            f"design ~ {' + '.join(design_factors_m)}): {len(contrasts)} pairwise contrast(s) computed."
                        )
                        st.success(f"Computed {len(contrasts)} pairwise DESeq2 contrast(s).")
                    except (normalization_deseq2.NormalizationError, deseq2_stats.StatsError) as e:
                        st.error(str(e))

                if st.session_state.get("stats_multi_deseq_contrasts"):
                    contrasts = st.session_state.stats_multi_deseq_contrasts
                    contrast_choice = st.selectbox("View contrast:", list(contrasts.keys()),
                                                    key="deseq_multi_contrast_view")
                    contrast_df = contrasts[contrast_choice]
                    st.session_state.stats_result = contrast_df
                    st.session_state.stats_result_groups = tuple(contrast_choice.split(" vs "))
                    st.session_state.stats_contrast_label = contrast_choice
                    n_sig = int(utils.apply_significance_filter(
                        contrast_df, sig_metric_m, sig_fdr_threshold_m, apply_filter=False)["Significant"].sum())
                    st.caption(f"{contrast_choice}: {n_sig} significant genes "
                               f"({sig_metric_m} ≤ {sig_fdr_threshold_m:g}).")
                    if contrast_df.attrs.get("lfc_shrunk"):
                        st.caption(
                            "Log2FC/Linear_FC are pyDESeq2's **apeGLM-shrunk** estimates; the "
                            "unshrunk MLE estimate is in Log2FC_MLE_Unshrunk / Linear_FC_MLE_Unshrunk. "
                            "Shrinkage does not change the p-value or FDR."
                        )
                    elif contrast_style.startswith("Every group vs"):
                        st.caption("LFC shrinkage was not applied to this contrast (reference-coefficient mismatch); "
                                   "showing the unshrunk MLE estimate.")
                    else:
                        st.caption(
                            "'All pairwise' contrasts aren't vs. a single reference, so apeGLM LFC "
                            "shrinkage (which shrinks one named design coefficient) doesn't apply here — "
                            "this is the unshrunk MLE estimate. Use 'Every group vs a reference' for "
                            "shrunk estimates."
                        )
                    # Each pairwise DESeq2 contrast here is itself a two-group Wald test
                    # result (same schema as the primary Two-Group Test above) -- NOT the
                    # F-statistic-based omnibus ANOVA table -- so the same 6-column
                    # display/export restriction applies. The full contrast_df (baseMean,
                    # lfcSE, stat, Group_A/Group_B, Higher_In, etc.) stays in session state
                    # for downstream Volcano/Biomarker/Heatmap use.
                    contrast_filtered = utils.apply_significance_filter(
                        contrast_df, sig_metric_m, sig_fdr_threshold_m, sig_apply_filter_m)
                    contrast_display = utils.two_group_display_view(contrast_filtered)
                    st.dataframe(utils.format_for_display(contrast_display), width='stretch', height=350)
                    safe_name = contrast_choice.replace(" ", "_")
                    st.download_button(
                        f"Download {safe_name}_DESeq2_Contrast.csv",
                        utils.to_download_bytes_csv(contrast_display),
                        file_name=f"{safe_name}_DESeq2_Contrast.csv", mime="text/csv",
                        key="dl_deseq_multi_one"
                    )
                    # Combined multi-contrast export: same 6 columns, plus a "Contrast"
                    # label column (necessary here since rows from different pairwise
                    # contrasts are concatenated and must remain distinguishable). Each
                    # contrast gets the same Significance Threshold filter as the single-
                    # contrast view above.
                    all_contrasts_df = pd.concat(
                        [utils.two_group_display_view(
                            utils.apply_significance_filter(df, sig_metric_m, sig_fdr_threshold_m, sig_apply_filter_m)
                         ).assign(Contrast=name)
                         for name, df in contrasts.items()]
                    )
                    st.download_button(
                        "Download All_DESeq2_Contrasts.csv",
                        utils.to_download_bytes_csv(all_contrasts_df),
                        file_name="All_DESeq2_Contrasts.csv", mime="text/csv",
                        key="dl_deseq_multi_all"
                    )
            else:
                anova_groups = st.multiselect(
                    f"{grouping_var} values to include in ANOVA (select 3 or more)",
                    groups_available, default=groups_available, key="anova_group_select"
                )
                posthoc_method = st.selectbox("Post-hoc test", ["tukey", "dunnett", "pairwise"],
                                               key="anova_posthoc_method")
                st.markdown("**Significance Threshold**")
                ac1, ac2, ac3 = st.columns(3)
                sig_metric_anova = ac1.selectbox(
                    "Metric", ["FDR", "P-value"], key="anova_sig_metric",
                    help="Which statistic the 'Significant' column and the filter below are based on."
                )
                sig_fdr_threshold_anova = ac2.number_input(
                    "Threshold", min_value=0.0001, max_value=1.0, value=0.05, step=0.01,
                    key="anova_sig_threshold",
                    help="Genes with the selected metric ≤ this value are marked Significant."
                )
                sig_apply_filter_anova = ac3.checkbox(
                    "Apply significance filter", value=False, key="anova_sig_apply_filter",
                    help="Shows only genes meeting the selected significance criterion in the "
                         "ANOVA table below and its CSV download."
                )

                if len(anova_groups) < 3:
                    st.warning(f"Select at least 3 groups to run ANOVA (currently {len(anova_groups)} selected).")
                elif st.button("Run ANOVA", key="run_anova_btn"):
                    group_map = meta.loc[sample_cols, grouping_var]
                    group_map = group_map[group_map.isin(anova_groups)]
                    group_map = group_map[group_map.index.isin(log2_df.columns)]
                    anova_table, posthoc_results = stats_analysis.anova_test(
                        log2_df[group_map.index], group_map, posthoc=posthoc_method,
                        fdr_threshold=sig_fdr_threshold_anova
                    )
                    anova_table = utils.apply_significance_filter(
                        anova_table, sig_metric_anova, sig_fdr_threshold_anova, apply_filter=False,
                        pval_col="ANOVA p-value")
                    st.session_state.anova_result = anova_table
                    st.session_state.posthoc_result = posthoc_results
                    st.session_state.anova_groups_used = list(anova_groups)
                    st.session_state.stats_result_group_col = grouping_var
                    st.session_state.processing_notes.append(
                        f"One-way ANOVA across {len(anova_groups)} selected {grouping_var} values "
                        f"({', '.join(anova_groups)}) with {posthoc_method} post-hoc; "
                        f"{int(anova_table['Significant'].sum())} significant genes "
                        f"({sig_metric_anova} ≤ {sig_fdr_threshold_anova:g})."
                    )
                    st.success(f"ANOVA complete: {int(anova_table['Significant'].sum())} significant genes "
                               f"({sig_metric_anova} ≤ {sig_fdr_threshold_anova:g}).")

                if st.session_state.get("anova_result") is not None:
                    anova_table = st.session_state.anova_result
                    anova_groups_used = st.session_state.get("anova_groups_used") or anova_groups
                    anova_comparison_name = "ANOVA_" + "_vs_".join(anova_groups_used)

                    # Also expose this as the active "stats_result" so the Heatmap can use the ANOVA
                    # significance. An omnibus test has no single fold change, so Log2FC/Linear_FC are
                    # left empty (not 0 / 1, which would look like "no change"); Volcano, Biomarker
                    # Discovery and GSEA ask for a two-group comparison instead.
                    stats_view = anova_table.rename(columns={"ANOVA p-value": "p-value"}).copy()
                    stats_view["Log2FC"] = np.nan
                    stats_view["Linear_FC"] = np.nan
                    st.session_state.stats_result = stats_view
                    st.session_state.stats_contrast_label = f"ANOVA across {', '.join(anova_groups_used)}"

                    # Dynamic Significance Threshold filter (Metric/Threshold/Apply widgets
                    # above), applied to both the on-screen ANOVA table and its CSV download.
                    anova_display = utils.apply_significance_filter(
                        anova_table, sig_metric_anova, sig_fdr_threshold_anova, sig_apply_filter_anova,
                        pval_col="ANOVA p-value")
                    st.dataframe(utils.format_for_display(anova_display), width='stretch', height=350)
                    st.download_button(
                        f"Download {anova_comparison_name}.csv",
                        utils.to_download_bytes_csv(anova_display),
                        file_name=f"{anova_comparison_name}.csv", mime="text/csv",
                        key="dl_stats_anova"
                    )
                    st.caption(f"Filename includes the groups compared ({', '.join(anova_groups_used)}) "
                               "so results from different ANOVA runs stay distinguishable. Showing "
                               f"{len(anova_display)} of {len(anova_table)} genes"
                               f"{' (filtered by ' + sig_metric_anova + ' ≤ ' + format(sig_fdr_threshold_anova, 'g') + ')' if sig_apply_filter_anova else ''}.")

                    posthoc_results = st.session_state.get("posthoc_result")
                    if posthoc_results:
                        feat_choice = st.selectbox("View post-hoc results for gene:", list(posthoc_results.keys()),
                                                    key="anova_posthoc_feat")
                        st.dataframe(utils.format_for_display(posthoc_results[feat_choice]), width='stretch')
                        st.download_button(
                            f"Download {anova_comparison_name}_Posthoc_{feat_choice}.csv",
                            utils.to_download_bytes_csv(posthoc_results[feat_choice]),
                            file_name=f"{anova_comparison_name}_Posthoc_{feat_choice}.csv", mime="text/csv",
                            key="dl_stats_posthoc_one"
                        )
                        all_posthoc = pd.concat(
                            [df.assign(Gene=feat) for feat, df in posthoc_results.items()],
                            ignore_index=True
                        )
                        st.download_button(
                            f"Download {anova_comparison_name}_Posthoc_All_Genes.csv",
                            utils.to_download_bytes_csv(all_posthoc),
                            file_name=f"{anova_comparison_name}_Posthoc_All_Genes.csv", mime="text/csv",
                            key="dl_stats_posthoc_all"
                        )

# ===========================================================================
# TAB 7 — VOLCANO PLOT
# ===========================================================================
if ACTIVE_PAGE == 6:
    st.header("Volcano Plot")
    if st.session_state.stats_result is None:
        st.info("Run a two-group statistical comparison in the Statistics page first.")
    elif str(st.session_state.get("stats_contrast_label") or "").startswith("ANOVA"):
        st.info("The current Statistics result is an ANOVA across three or more groups, which has no single "
                "fold change. The volcano plot needs a two-group comparison: on the Statistics page, run a two-group test "
                "or a DESeq2 contrast between two groups.")
    else:
        result = st.session_state.stats_result
        groups_used = st.session_state.get("stats_result_groups")
        group_col_used = st.session_state.get("stats_result_group_col") or "Group"
        if groups_used:
            st.caption(
                f"Reflects the two {group_col_used} values compared in the Statistics page: "
                f"**{groups_used[0]}** vs **{groups_used[1]}**."
                + " To compare a different pair or variable, go back to Statistics and "
                "re-run the two-group test."
            )

        # ---- 1. Thresholds ----
        with st.expander("1 — Statistical Thresholds", expanded=True):
            c1, c2, c3, c4 = st.columns(4)
            y_metric_label = c1.radio("Y-axis metric", ["p-value", "FDR"], horizontal=False, key="volc_ymetric_label")
            y_metric = "pvalue" if y_metric_label == "p-value" else "fdr"
            sig_cutoff = c2.number_input("Significance cutoff", value=0.05,
                                          min_value=0.0001, max_value=1.0, step=0.01, key="volc_cutoff")
            fc_threshold = c3.number_input("Fold change threshold (log2 units)", value=0.0,
                                            min_value=0.0, max_value=5.0, step=0.1, key="volc_fc")
            use_fdr_secondary = c4.checkbox("Also require FDR <", value=False, key="volc_use_fdr2")
            fdr_cutoff = c4.number_input("Secondary FDR cutoff", value=0.05, min_value=0.0001, max_value=1.0,
                                          step=0.01, disabled=not use_fdr_secondary, key="volc_fdr2") \
                if use_fdr_secondary else None
            c5, c6 = st.columns(2)
            use_fc_axis = c5.checkbox("Show linear FC on x-axis instead of log2FC", value=False, key="volc_fc_axis")
            threshold_line_style = c6.selectbox("Threshold line style", list(volcano.LINE_STYLES.keys()),
                                                 key="volc_threshold_style")

        # ---- 2 & 3. Colors and point style ----
        with st.expander("2-3 — Colors & Point Style"):
            c1, c2 = st.columns(2)
            palette = c1.selectbox("Color palette (colorblind-friendly options included)",
                                    list(volcano.COLORBLIND_PALETTES.keys()), key="volc_palette")
            override_colors = c2.checkbox("Override individual colors", value=False, key="volc_override_colors")
            up_color = down_color = ns_color = None
            if override_colors:
                cc1, cc2, cc3 = st.columns(3)
                with cc1:
                    up_color = st.color_picker("Upregulated color", volcano.COLORBLIND_PALETTES[palette]["up"],
                                                key="volcano_up_color")
                with cc2:
                    down_color = st.color_picker("Downregulated color", volcano.COLORBLIND_PALETTES[palette]["down"],
                                                  key="volcano_down_color")
                with cc3:
                    ns_color = st.color_picker("Not significant color", volcano.COLORBLIND_PALETTES[palette]["ns"],
                                                key="volcano_ns_color")
            c3, c4, c5 = st.columns(3)
            alpha = c3.slider("Point transparency (alpha)", 0.2, 1.0, 0.75, 0.05, key="volc_alpha")
            point_size = c4.slider("Point size", 1, 100, 14, 1, key="volc_point_size")
            point_shape = c5.selectbox("Point shape", list(volcano.MARKER_SHAPES.keys()), key="volc_point_shape")
            c6, c7 = st.columns(2)
            edge_width = c6.slider("Point border width", 0.0, 2.0, 0.0, 0.1, key="volc_edge_width")
            edge_color = "none"
            if edge_width > 0:
                with c7:
                    edge_color = st.color_picker("Point border color", "#000000", key="volcano_edge_color")

        # ---- 4. Labels ----
        with st.expander("4 — Gene Labels"):
            label_mode_label = st.radio(
                "Label mode",
                ["Top N significant", "All significant", "Manually selected", "None"],
                horizontal=True, key="volc_label_mode_label"
            )
            label_mode = {"Top N significant": "top_n", "All significant": "significant_only",
                          "Manually selected": "manual", "None": "none"}[label_mode_label]
            top_label_n = 10
            manual_labels = []
            if label_mode == "top_n":
                top_label_n = st.number_input("Label top N genes", value=10, min_value=0, max_value=100, step=1,
                                               key="volc_top_n")
            elif label_mode == "manual":
                manual_labels = st.multiselect("Search and select genes to label by name",
                                                result.index.tolist(), key="volc_manual_labels")
            c1, c2, c3 = st.columns(3)
            label_font_size = c1.slider("Label font size", 5.0, 16.0, 7.5, 0.5, key="volc_label_fontsize")
            label_bold = c2.checkbox("Bold labels", value=False, key="volc_label_bold")
            label_italic = c3.checkbox("Italic labels", value=False, key="volc_label_italic")
            c4, c5 = st.columns(2)
            override_label_color = c4.checkbox("Override label color (default: match point color)", value=False,
                                                key="volc_override_label_color")
            label_color = None
            if override_label_color:
                with c5:
                    label_color = st.color_picker("Label color", "#333333", key="volcano_label_color")
            repel_labels = st.checkbox("Repel overlapping labels automatically", value=True, key="volc_repel")

        # ---- 5. Legend ----
        with st.expander("5 — Legend"):
            c1, c2 = st.columns(2)
            legend_position = c1.selectbox("Legend position", ["Right", "Left", "Top", "Bottom", "Hidden"],
                                            key="volc_legend_pos")
            legend_format = c2.radio("Legend label format",
                                      ["Full (\"Significantly Upregulated (n=5)\")", "Short (\"Up (5)\")"], horizontal=False,
                                      key="volc_legend_fmt")
            legend_format_key = "full" if legend_format.startswith("Full") else "short"

        # ---- 6 & 7. Axes and title ----
        with st.expander("6-7 — Axes & Title"):
            st.markdown("**Axes**")
            c1, c2, c3 = st.columns(3)
            axis_font_size = c1.slider("Axis font size", 8.0, 20.0, 12.0, 0.5, key="volc_axis_fontsize")
            axis_bold = c2.checkbox("Bold axis titles", value=False, key="volc_axis_bold")
            y_decimals = c3.selectbox("Y-axis decimals", [None, 0, 1, 2, 3], index=0,
                                      format_func=lambda x: "Auto" if x is None else str(x), key="volc_y_decimals")
            c4, c5 = st.columns(2)
            custom_x_limits = c4.checkbox("Set custom X-axis limits", value=False, key="volc_custom_xlim")
            x_limits = None
            if custom_x_limits:
                xlo, xhi = st.columns(2)
                x_limits = (xlo.number_input("X min", value=-6.0, key="volc_xmin"),
                            xhi.number_input("X max", value=6.0, key="volc_xmax"))
            custom_y_limits = c5.checkbox("Set custom Y-axis limits", value=False, key="volc_custom_ylim")
            y_limits = None
            if custom_y_limits:
                ylo, yhi = st.columns(2)
                y_limits = (ylo.number_input("Y min", value=0.0, key="volc_ymin"),
                            yhi.number_input("Y max", value=10.0, key="volc_ymax"))
            x_tick_spacing = st.number_input("X-axis tick spacing (0 = auto)", value=0.0, min_value=0.0, step=0.5,
                                              key="volc_xtick")
            x_tick_spacing = x_tick_spacing or None

            st.markdown("**Title**")
            c6, c7 = st.columns(2)
            default_title = f"{groups_used[0]} vs {groups_used[1]}" if groups_used else "Volcano Plot"
            custom_title = c6.text_input("Title", value=default_title, key="volc_title")
            subtitle = c7.text_input("Subtitle (optional)", value="", key="volc_subtitle")
            c8, c9, c10, c11 = st.columns(4)
            title_bold = c8.checkbox("Bold title", value=True, key="volc_title_bold")
            title_italic = c9.checkbox("Italic title", value=False, key="volc_title_italic")
            title_align = c10.selectbox("Title alignment", ["center", "left", "right"], key="volc_title_align")
            hide_title = c11.checkbox("Hide title", value=False, key="volc_hide_title")

        # ---- 8 & 9. Gridlines and threshold line style ----
        with st.expander("8-9 — Gridlines & Threshold Line Style"):
            c1, c2, c3 = st.columns(3)
            grid_mode = c1.selectbox("Gridlines", ["none", "major", "both"],
                                     format_func=lambda x: {"none": "No grid", "major": "Major only",
                                                              "both": "Major + minor"}[x], key="volc_grid_mode")
            with c2:
                grid_color = st.color_picker("Grid color", "#D9D9D9", key="volcano_grid_color")
            grid_style = c3.selectbox("Grid line style", list(volcano.LINE_STYLES.keys()), index=0,
                                       key="volc_grid_style")
            c4, c5 = st.columns(2)
            with c4:
                threshold_line_color = st.color_picker("Threshold line color", "#808080",
                                                         key="volcano_threshold_color")
            threshold_line_width = c5.slider("Threshold line width", 0.2, 3.0, 0.7, 0.1, key="volc_threshold_width")

        # ---- 10. Highlight specific genes ----
        with st.expander("10 — Highlight Specific Genes"):
            highlight_names = st.multiselect(
                "Search and select genes to highlight (e.g. known markers)",
                result.index.tolist(), key="volc_highlight_names"
            )
            c1, c2, c3 = st.columns(3)
            with c1:
                highlight_color = st.color_picker("Highlight color", "#FFD700", key="volcano_highlight_color")
            highlight_size_mult = c2.slider("Highlight size multiplier", 1.0, 5.0, 2.0, 0.25, key="volc_hl_size")
            highlight_shape = c3.selectbox("Highlight shape", ["Star", "Circle", "Triangle", "Square", "Diamond"],
                                            key="volc_hl_shape")

        # ---- 12 & 13. Background and figure size ----
        with st.expander("12-13 — Background & Figure Size"):
            c1, c2 = st.columns(2)
            background_choice = c1.selectbox("Background", ["White", "Transparent", "Gray", "Custom"],
                                              key="volc_bg_choice")
            background = {"White": "white", "Transparent": "transparent", "Gray": "gray"}.get(background_choice)
            if background_choice == "Custom":
                with c2:
                    background = st.color_picker("Custom background color", "#FFFFFF", key="volcano_bg_color")
            c3, c4 = st.columns(2)
            size_preset = c3.selectbox("Figure size preset", list(volcano.FIGURE_SIZE_PRESETS.keys()), index=2,
                                        key="volc_size_preset")
            if size_preset == "Custom":
                cw, ch = st.columns(2)
                fig_width_in = cw.number_input("Width (in)", value=7.2, min_value=2.0, max_value=20.0,
                                                key="volc_width")
                fig_height_in = ch.number_input("Height (in)", value=6.4, min_value=2.0, max_value=20.0,
                                                 key="volc_height")
            else:
                fig_width_in, fig_height_in = volcano.FIGURE_SIZE_PRESETS[size_preset]

        # ---- 15. Statistics box, theme ----
        with st.expander("Statistics Display & Publication Theme"):
            c1, c2 = st.columns(2)
            show_stats_box = c1.checkbox("Show statistics box on figure (total/up/down/cutoffs)", value=False,
                                          key="volc_stats_box")
            theme_choice = c2.selectbox("Publication theme (stylistic approximation, not an official spec)",
                                        list(volcano.THEMES.keys()), key="volc_theme")

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
            show_stats_box=show_stats_box, theme=theme_choice,
        )
        st.session_state.volcano_fig = fig_volc
        st.session_state.volcano_annotated = annotated
        st.session_state.volcano_settings = {
            "y_metric": y_metric, "sig_cutoff": sig_cutoff, "fc_threshold": fc_threshold,
            "fdr_cutoff": fdr_cutoff, "palette": palette, "point_size": point_size,
            "point_shape": point_shape, "label_mode": label_mode, "legend_position": legend_position,
            "theme": theme_choice, "groups": groups_used,
        }
        utils.render_figure(st, fig_volc)

        ec1, ec2, ec3 = st.columns(3)
        export_fmt = ec1.selectbox("Format", ["PNG", "PDF", "SVG", "JPEG", "TIFF"], key="volcano_export_fmt")
        export_dpi = ec2.selectbox("Resolution (DPI)", [300, 600, 1200], index=0,
                                    disabled=export_fmt in ("PDF", "SVG"), key="volcano_export_dpi")
        mime_map = {"PNG": "image/png", "PDF": "application/pdf", "SVG": "image/svg+xml",
                    "JPEG": "image/jpeg", "TIFF": "image/tiff"}
        ext_map = {"PNG": "png", "PDF": "pdf", "SVG": "svg", "JPEG": "jpg", "TIFF": "tiff"}
        fmt_key = "jpeg" if export_fmt == "JPEG" else export_fmt.lower()
        file_bytes = volcano.export_figure(fig_volc, fmt=fmt_key, dpi=export_dpi)
        volcano_tag = (f"{_safe_tag(groups_used[0])}_vs_{_safe_tag(groups_used[1])}"
                       if groups_used else "Comparison")
        volcano_plot_filename = f"{volcano_tag}_VolcanoPlot.{ext_map[export_fmt]}"
        with ec3:
            st.write("")
            st.write("")
            st.download_button(
                f"Download {volcano_plot_filename}", file_bytes,
                file_name=volcano_plot_filename, mime=mime_map[export_fmt],
                key="dl_volcano_fig"
            )

        sig_col = "p-value" if y_metric == "pvalue" else "FDR"

        st.subheader("Top 20 Biomarkers")
        st.dataframe(utils.format_for_display(volcano.top_biomarker_labels(annotated, sig_col=sig_col, n=20)), width='stretch')

        st.subheader("Export Data")
        dc1, dc2, dc3 = st.columns(3)
        up_tbl = volcano.get_direction_table(annotated, "Up")
        down_tbl = volcano.get_direction_table(annotated, "Down")
        with dc1:
            st.download_button(f"Download {volcano_tag}_Upregulated_Genes.csv", utils.to_download_bytes_csv(up_tbl),
                                f"{volcano_tag}_Upregulated_Genes.csv", "text/csv", key="dl_volcano_up")
        with dc2:
            st.download_button(f"Download {volcano_tag}_Downregulated_Genes.csv", utils.to_download_bytes_csv(down_tbl),
                                f"{volcano_tag}_Downregulated_Genes.csv", "text/csv", key="dl_volcano_down")
        with dc3:
            st.download_button(f"Download {volcano_tag}_Complete_Volcano_Data.csv", utils.to_download_bytes_csv(annotated),
                                f"{volcano_tag}_Complete_Volcano_Data.csv", "text/csv", key="dl_volcano_all")
        st.caption("Figure settings (for reproducibility):")
        settings_json = volcano.export_settings_json(st.session_state.volcano_settings)
        st.download_button(f"Download {volcano_tag}_Volcano_Figure_Settings.json", settings_json.encode(),
                            f"{volcano_tag}_Volcano_Figure_Settings.json", "application/json", key="dl_volcano_settings")

# ===========================================================================
# TAB 8 — BIOMARKER DISCOVERY
# ===========================================================================
if ACTIVE_PAGE == 7:
    st.header("Biomarker Discovery")
    if st.session_state.stats_result is None:
        st.info("Run a comparison in the Statistics page first.")
    elif str(st.session_state.get("stats_contrast_label") or "").startswith("ANOVA"):
        st.info("The current Statistics result is an ANOVA across three or more groups, which has no single "
                "fold change. Biomarker Discovery (fold change and direction) needs a two-group comparison: on the Statistics page, run a two-group test "
                "or a DESeq2 contrast between two groups.")
    else:
        c1, c2, c3 = st.columns(3)
        # p-value/FDR cutoffs are defined before the criterion selectbox (though laid
        # out in column c1) so the criterion labels below can be built dynamically
        # from their current values — never hardcoded to "< 0.05", so raising either
        # threshold immediately updates the label text.
        p_cutoff = c2.number_input("p-value cutoff", 0.0001, 0.5, 0.05, 0.005, key="bio_pcut")
        fdr_cutoff = c3.number_input("FDR cutoff", 0.0001, 1.0, 0.05, 0.01, key="bio_fdrcut")
        criterion_labels = {
            "combined": f"FDR < {fdr_cutoff:g} and P value < {p_cutoff:g}",
            "pvalue": f"P value < {p_cutoff:g}",
            "fdr": f"FDR < {fdr_cutoff:g}",
        }
        criterion = c1.selectbox(
            "Selection criterion", ["combined", "pvalue", "fdr"],
            format_func=lambda k: criterion_labels[k], key="bio_criterion"
        )

        if st.button("Run Biomarker Discovery", type="primary", key="run_biomarker_btn"):
            st.session_state.biomarkers = biomarker.discover_biomarkers(
                st.session_state.stats_result, criterion=criterion, p_cutoff=p_cutoff, fdr_cutoff=fdr_cutoff
            )
            st.session_state.biomarker_criterion_label = criterion_labels[criterion]

        if st.session_state.biomarkers is not None:
            bm = st.session_state.biomarkers
            applied_label = st.session_state.get("biomarker_criterion_label", criterion_labels[criterion])
            st.success(f"{bm.shape[0]} candidate biomarker genes identified (criterion: {applied_label}).")

            row_annot = st.session_state.row_annotations
            # Optional gene-annotation left-join: purely informational extra columns
            # (e.g. Pathway, Gene_Symbol, Description) from the loaded gene
            # row-annotation file, if any -- never alters the statistics above.
            # Non-breaking: identical to current behavior when no annotation file
            # was loaded (row_annot is None).
            bm_display = bm.join(row_annot, how="left") if row_annot is not None else bm
            st.dataframe(utils.format_for_display(bm_display), width='stretch', height=350)
            if row_annot is not None and "Pathway" in row_annot.columns:
                merged = bm.join(row_annot[["Pathway"]], how="left")
                st.caption("Pathway breakdown of candidate biomarkers:")
                st.bar_chart(merged["Pathway"].value_counts())
            groups_used = st.session_state.get("stats_result_groups")
            comparison_name = (f"{groups_used[0]}_vs_{groups_used[1]}" if groups_used else "Comparison")
            biomarker_filename = f"{comparison_name}_Biomarkers.csv"
            st.download_button(
                f"Download {biomarker_filename}", utils.to_download_bytes_csv(bm_display),
                file_name=biomarker_filename, mime="text/csv", key="dl_biomarkers"
            )

# ===========================================================================
# TAB 9 — HEATMAP
# ===========================================================================
if ACTIVE_PAGE == 8:
    st.header("Clustered Heatmap")
    if st.session_state.stats_result is None:
        st.warning(
            "Run a comparison in the Statistics page first — the heatmap needs "
            "that result to know which genes are significant."
        )
    else:
        log2_df_full = st.session_state.log2_data
        result = st.session_state.stats_result
        meta = st.session_state.meta

        heat_cat_cols = utils.get_categorical_metadata_columns(meta, log2_df_full.columns)
        default_group_col = st.session_state.get("stats_result_group_col") or "Group"
        heat_group_col = st.selectbox(
            "Filter samples by:", heat_cat_cols,
            index=heat_cat_cols.index(default_group_col) if default_group_col in heat_cat_cols else 0,
            key="heatmap_group_col"
        )
        groups_all_heat = meta.loc[meta.index.intersection(log2_df_full.columns), heat_group_col].dropna().unique().tolist()
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
        if row_annotations is not None:
            st.subheader("Row (Gene) Annotation")
            row_annot_cols = st.multiselect(
                "Row annotation tracks:", list(row_annotations.columns),
                default=list(row_annotations.columns)[:1], key="heatmap_row_annot_cols"
            )
            if row_annot_cols:
                n_matched = heatmap_module.count_row_annotation_matches(row_annotations, log2_df_full.index)
                st.caption(f"{n_matched}/{len(log2_df_full.index)} genes matched to row annotations "
                           f"by name; unmatched genes show as '(unannotated)' in the row bar(s).")

        annotation_colors = {}
        if col_annot_cols or row_annot_cols:
            with st.expander("🎨 Customize Annotation Track Colors", expanded=False):
                st.caption(
                    "Assign specific colors per category (or a colormap for numeric tracks), per "
                    "annotation track — updates the column/row annotation bars on the heatmap below "
                    "dynamically as you change them."
                )
                MAX_VALUES_FOR_COLOR_UI = 12

                def _annotation_color_controls(col, series, axis_label):
                    kind = heatmap_module.classify_annotation_series(series)
                    if kind == "continuous":
                        st.markdown(f"**{col}** ({axis_label} track, continuous)")
                        return st.selectbox(
                            f"Colormap for {col}", heatmap_module.CONTINUOUS_ANNOT_CMAPS,
                            key=f"heatmap_annotcmap_{axis_label}_{col}"
                        )
                    values = list(series.dropna().unique())
                    if len(values) > MAX_VALUES_FOR_COLOR_UI:
                        st.caption(f"**{col}** ({axis_label} track) has {len(values)} distinct values — "
                                   "too many to customize individually; using auto-assigned colors.")
                        return None
                    st.markdown(f"**{col}** ({axis_label} track)")
                    value_cols = st.columns(min(4, len(values)) or 1)
                    track_colors = {}
                    for i, v in enumerate(values):
                        default_c = heatmap_module.GROUP_PALETTE[i % len(heatmap_module.GROUP_PALETTE)]
                        with value_cols[i % len(value_cols)]:
                            track_colors[str(v)] = st.color_picker(
                                str(v), default_c, key=f"heatmap_annotcolor_{axis_label}_{col}_{v}"
                            )
                    return track_colors

                for col in col_annot_cols:
                    result_colors = _annotation_color_controls(col, meta.loc[heat_cols, col], "col")
                    if result_colors:
                        annotation_colors[col] = result_colors
                for col in row_annot_cols:
                    result_colors = _annotation_color_controls(col, row_annotations[col], "row")
                    if result_colors:
                        annotation_colors[col] = result_colors

        st.subheader("Feature Selection")
        st.caption("Dynamically adjust which genes appear in the heatmap by FDR or p-value — "
                   "the figure below updates immediately as you change the metric or threshold.")
        c1, c2, c3 = st.columns([1.2, 1.5, 1])
        metric_options = [m for m in ["FDR", "p-value"] if m in result.columns]
        if not metric_options:
            st.error("The current Statistics result has no FDR or p-value column to filter on.")
            st.stop()
        sig_metric = c1.selectbox("Filter by", metric_options, key="heatmap_sig_metric")
        sig_threshold = c2.slider(
            f"{sig_metric} ≤", min_value=0.0, max_value=1.0,
            value=0.05, step=0.01, key="heatmap_sig_threshold"
        )
        n_preview = (result[sig_metric] <= sig_threshold).sum()
        with c3:
            st.metric("Genes at this cutoff", n_preview)
        if n_preview > heatmap_module.HIDE_ROW_LABELS_ABOVE:
            st.warning(
                f"{n_preview} genes selected — beyond {heatmap_module.HIDE_ROW_LABELS_ABOVE}, "
                f"row labels become unreadable and will be hidden automatically. Tighten the cutoff "
                f"above for a labeled figure, or proceed for an unlabeled overview heatmap."
            )

        st.subheader("Clustering Options")
        c3, c4, c5 = st.columns(3)
        cluster_mode = c3.selectbox("Clustering", ["Rows only", "Both rows and columns", "Columns only",
                                                     "No clustering"], key="heatmap_cluster_mode")
        distance = c4.selectbox("Distance metric", ["euclidean", "correlation"], key="heatmap_distance")
        linkage_m = c5.selectbox("Linkage method", ["ward", "average", "complete"], key="heatmap_linkage")
        cluster_rows = cluster_mode in ("Both rows and columns", "Rows only")
        cluster_cols = cluster_mode in ("Both rows and columns", "Columns only")

        st.subheader("Color Scale Customization")
        c6, c7 = st.columns(2)
        use_custom_gradient = c6.checkbox("Use a custom color gradient instead of a preset palette", value=False,
                                           key="heatmap_use_gradient")
        reverse_cmap = c7.checkbox("Reverse colormap", value=False, key="heatmap_reverse_cmap")

        if use_custom_gradient:
            gc1, gc2, gc3 = st.columns(3)
            with gc1:
                color_low = st.color_picker("Low color", "#2166AC", key="heatmap_grad_low")
            with gc2:
                color_mid = st.color_picker("Mid color", "#FFFFFF", key="heatmap_grad_mid")
            with gc3:
                color_high = st.color_picker("High color", "#B2182B", key="heatmap_grad_high")
            custom_colors = [color_low, color_mid, color_high]
            cmap_name = None
        else:
            cmap_name = st.selectbox("Predefined color palette", heatmap_module.PREDEFINED_PALETTES,
                                      key="heatmap_cmap_name")
            custom_colors = None

        c8, c9 = st.columns(2)
        vmin = c8.number_input("Z-score minimum", value=-2.5, step=0.1, key="heatmap_vmin")
        vmax = c9.number_input("Z-score maximum", value=2.5, step=0.1, key="heatmap_vmax")

        use_breakpoints = st.checkbox("Use custom discrete color breakpoints (instead of a continuous scale)",
                                       value=False, key="heatmap_use_breakpoints")
        breakpoints = None
        if use_breakpoints:
            bp_text = st.text_input("Comma-separated breakpoints (e.g. -3,-1,0,1,3)",
                                     value=f"{vmin},{vmin/2:.2g},0,{vmax/2:.2g},{vmax}", key="heatmap_bp_text")
            try:
                breakpoints = [float(x.strip()) for x in bp_text.split(",") if x.strip()]
            except ValueError:
                st.warning("Couldn't parse breakpoints — using a continuous scale instead.")
                breakpoints = None

        st.subheader("Figure Size")
        c_size1, c_size2, c_size3 = st.columns(3)
        use_custom_size = c_size1.checkbox("Customize heatmap panel size", value=False, key="heatmap_custom_size")
        heatmap_width_in = heatmap_height_in = None
        if use_custom_size:
            heatmap_width_in = c_size2.number_input("Width (inches)", value=8.0, min_value=2.0, max_value=30.0,
                                                      step=0.5, key="heatmap_width")
            heatmap_height_in = c_size3.number_input("Height (inches)", value=6.0, min_value=2.0, max_value=30.0,
                                                       step=0.5, key="heatmap_height")
            st.caption("All other elements (dendrograms, annotation bar, legend, colorbar, labels) "
                       "scale automatically with the size you set here.")

        st.subheader("Title & Font Size")
        ft1, ft2 = st.columns(2)
        font_family_ui = ft1.selectbox("Font family", ["sans-serif", "serif", "monospace"], key="heatmap_font_family")
        font_size_ui = ft2.slider("Base font size", 6, 18, 10, key="heatmap_font_size",
                                   help="Applied proportionally to the title, tick labels, legend text, "
                                        "colorbar label, and annotation track labels.")
        st.caption("Row labels still shrink automatically for very large gene counts, capped at the "
                   "size chosen above — the cap doesn't force oversized text when hundreds of rows are shown.")

        if len(heat_cols) < 2:
            st.warning("Select at least one group with 2+ samples to build a heatmap.")
        else:
            sig_feats = result[result[sig_metric] <= sig_threshold].index

            if len(sig_feats) < 2:
                st.warning(
                    f"Only {len(sig_feats)} significant gene(s) at this cutoff — "
                    "clustering needs at least 2. Relax the threshold above."
                )
            else:
                try:
                    title_suffix = f"{len(sig_feats)} genes, {sig_metric} ≤ {sig_threshold}"
                    fig_heat, z_ordered, notes = heatmap_module.clustered_heatmap(
                        log2_df, sig_feats, meta=meta, col_annot_cols=col_annot_cols,
                        cluster_rows=cluster_rows, cluster_cols=cluster_cols,
                        distance=distance, linkage_method=linkage_m,
                        cmap_name=cmap_name, custom_colors=custom_colors, reverse_cmap=reverse_cmap,
                        vmin=vmin, vmax=vmax, breakpoints=breakpoints,
                        heatmap_width_in=heatmap_width_in, heatmap_height_in=heatmap_height_in,
                        row_meta=row_annotations, row_annot_cols=row_annot_cols,
                        annotation_colors=annotation_colors,
                        font_family=font_family_ui, font_size=font_size_ui, title_suffix=title_suffix,
                    )
                    st.session_state.heatmap_fig = fig_heat
                    for note in notes:
                        st.info(f"ℹ️ {note}")
                    utils.render_figure(st, fig_heat)
                    st.caption(f"{len(sig_feats)} genes shown ({sig_metric} ≤ {sig_threshold}), "
                               f"{len(heat_cols)} samples across {len(selected_heat_groups)} selected group(s), "
                               f"row-scaled (z-score).")

                    heatmap_tag = "_".join(_safe_tag(g) for g in selected_heat_groups) or "AllGroups"


                    st.subheader("Export Heatmap")
                    ec1, ec2, ec3 = st.columns(3)
                    export_fmt = ec1.selectbox("Format", ["PNG", "PDF", "SVG", "JPEG", "TIFF"],
                                                key="heatmap_export_fmt")
                    export_dpi = ec2.selectbox("Resolution (DPI)", [150, 300, 600], index=1,
                                                key="heatmap_export_dpi")
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
                            file_name=heatmap_filename, mime=mime_map[export_fmt],
                            key="dl_heatmap_fig"
                        )

                    st.subheader("Z-score Data Table")
                    st.caption(
                        "The exact row-scaled Z-score values used to render the heatmap above "
                        "(rows/columns in the same clustered order shown)."
                    )
                    st.dataframe(utils.format_for_display(z_ordered), width='stretch', height=300)
                    st.download_button(
                        f"Download {heatmap_tag}_Zscore_Table.csv", utils.to_download_bytes_csv(z_ordered),
                        file_name=f"{heatmap_tag}_Zscore_Table.csv", mime="text/csv", key="dl_heatmap_zscore"
                    )
                except ValueError as e:
                    st.warning(str(e))

# ===========================================================================
# TAB 10 — BOXPLOT
# ===========================================================================
if ACTIVE_PAGE == 9:
    st.header("Boxplot of Genes")
    if st.session_state.log2_data is None:
        st.info("Run Normalization first.")
    else:
        log2_df_full = st.session_state.log2_data
        meta = st.session_state.meta

        st.subheader("Select Genes & Groups")
        selected_genes = st.multiselect(
            "Search and select one or more genes",
            log2_df_full.index.tolist(), key="boxplot_genes"
        )
        box_cat_cols = utils.get_categorical_metadata_columns(meta, log2_df_full.columns)
        box_group_col = st.selectbox(
            "Grouping variable:", box_cat_cols,
            index=box_cat_cols.index("Diagnosis") if "Diagnosis" in box_cat_cols
            else (box_cat_cols.index("Group") if "Group" in box_cat_cols else 0),
            key="boxplot_group_col"
        )
        groups_all_box = meta.loc[meta.index.intersection(log2_df_full.columns), box_group_col].dropna().unique().tolist()
        selected_box_groups = st.multiselect(
            f"{box_group_col} values to include in the comparison (choose any subset — 2, 3, 4, or more)",
            groups_all_box, default=groups_all_box, key="boxplot_group_select"
        )

        if selected_genes and len(selected_box_groups) >= 2:
            test_name = "Welch's t-test" if len(selected_box_groups) == 2 else "one-way ANOVA"
            st.caption(
                f"Statistics computed on the normalized expression matrix: {test_name} "
                f"across the {len(selected_box_groups)} selected group(s). FDR here is corrected across "
                f"only the {len(selected_genes)} gene(s) shown, not the full gene panel — "
                f"for a panel-wide FDR, use the Statistics page."
            )

        st.subheader("Customize Appearance")
        c1, c2, c3 = st.columns(3)
        font_family = c1.selectbox("Font family", ["sans-serif", "serif", "monospace"], key="boxplot_font_family")
        font_size = c2.slider("Font size", 6, 20, 10, 1, key="boxplot_font_size")
        ncols = c3.number_input("Panels per row", min_value=1, max_value=6, value=3, step=1, key="boxplot_ncols")

        c4, c5 = st.columns(2)
        use_custom_box_size = c4.checkbox("Customize figure size", value=False, key="boxplot_custom_size")
        fig_width_in = fig_height_in = None
        if use_custom_box_size:
            fig_width_in = c5.number_input("Width (inches)", value=10.0, min_value=3.0, max_value=30.0, step=0.5,
                                            key="boxplot_width")
            fig_height_in = st.number_input("Height (inches)", value=6.0, min_value=3.0, max_value=30.0, step=0.5,
                                             key="boxplot_height")

        c6, c7, c8 = st.columns(3)
        show_points = c6.checkbox("Show individual data points", value=True, key="boxplot_points")
        show_mean = c7.checkbox("Show mean (dashed line)", value=False, key="boxplot_mean")
        show_median = c8.checkbox("Show median (solid line)", value=True, key="boxplot_median")

        st.caption("Optional: assign a color per group (defaults to a standard palette).")
        group_colors = {}
        if selected_box_groups:
            color_cols = st.columns(min(4, len(selected_box_groups)) or 1)
            for i, g in enumerate(selected_box_groups):
                with color_cols[i % len(color_cols)]:
                    default_c = boxplot_module.GROUP_PALETTE[i % len(boxplot_module.GROUP_PALETTE)]
                    group_colors[g] = st.color_picker(f"Color: {g}", default_c, key=f"boxcolor_{g}")

        if not selected_genes:
            st.info("Select at least one gene above to generate a boxplot.")
        elif len(selected_box_groups) < 2:
            st.warning("Select at least 2 groups to compare.")
        else:
            box_cols = [c for c in log2_df_full.columns if meta.loc[c, box_group_col] in selected_box_groups]
            log2_df_box = log2_df_full[box_cols]

            stats_table = boxplot_module.compute_stats_for_genes(
                log2_df_box, meta, selected_genes, selected_box_groups, group_col=box_group_col
            )
            fig_box = boxplot_module.boxplot_genes(
                log2_df_box, meta, selected_genes, selected_box_groups,
                stats_table=stats_table, fig_width_in=fig_width_in, fig_height_in=fig_height_in,
                font_size=font_size, font_family=font_family, group_colors=group_colors,
                show_points=show_points, show_mean=show_mean, show_median=show_median, ncols=ncols,
                group_col=box_group_col,
            )
            st.session_state.boxplot_fig = fig_box
            st.session_state.boxplot_stats = stats_table
            utils.render_figure(st, fig_box)

            boxplot_tag = "_".join(_safe_tag(g) for g in selected_box_groups) or "AllGroups"
            if len(selected_genes) == 1:
                boxplot_tag = f"{_safe_tag(selected_genes[0])}_{boxplot_tag}"

            st.subheader("Statistical Results")
            st.dataframe(utils.format_for_display(stats_table), width='stretch')
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
                    file_name=boxplot_filename, mime=mime_map[export_fmt],
                    key="dl_boxplot_fig"
                )

# ===========================================================================
# TAB 11 — GSEA (Gene Set Enrichment Analysis)
# ===========================================================================
if ACTIVE_PAGE == 10:
    st.header("Gene Set Enrichment Analysis")
    st.caption(
        "Three statistically-distinct methods: **STRING Enrichment Analysis** (Szklarczyk et "
        "al. 2023) needs only a gene list and computes enrichment server-side against STRING's "
        "own database; **Over-Representation Analysis (ORA)** is a hypergeometric test computed "
        "locally against a chosen gene-set library, using your full comparison list against an "
        "explicit background/universe; **Gene Set Enrichment Analysis (GSEA)** (Subramanian et "
        "al. 2005) ranks every detected gene and tests where each gene set falls in that ranking "
        "— a genuinely different algorithm from ORA, never used as a substitute for it. "
        "**STRING and ORA/GSEA's Enrichr-backed libraries require internet access**; GSEA can "
        "run fully offline if you upload your own .gmt gene set file."
    )

    if st.session_state.stats_result is None:
        st.warning("Run a comparison in the Statistics page first — GSEA analyzes the resulting "
                   "Log2FC/p-value/FDR gene list.")
    elif str(st.session_state.get("stats_contrast_label") or "").startswith("ANOVA"):
        st.info("The current Statistics result is an ANOVA across three or more groups, which has no single "
                "fold change. Enrichment analysis needs a two-group comparison: on the Statistics page, run a "
                "two-group test or a DESeq2 contrast between two groups.")
    else:
        stats_df = st.session_state.stats_result
        groups_used = st.session_state.get("stats_result_groups")
        group_col_used = st.session_state.get("stats_result_group_col") or "Group"
        if groups_used:
            st.caption(f"Reflects the two {group_col_used} values compared in the Statistics page: "
                       f"**{groups_used[0]}** vs **{groups_used[1]}**.")
        st.caption(
            "RNA-seq features are already gene identifiers, so no separate identifier→gene mapping step is "
            "needed here — the gene names in the expression matrix are used directly for "
            "STRING/Enrichr lookups (works best when they are HGNC symbols rather than Ensembl IDs)."
        )

        # ---- Significant-gene selector (shared by STRING and ORA) ----
        c1, c2 = st.columns(2)
        gsea_metric_options = [m for m in ["FDR", "p-value"] if m in stats_df.columns]
        sig_metric = c1.selectbox("Filter by", gsea_metric_options, key="gsea_sig_metric")
        sig_threshold = c2.number_input(f"{sig_metric} <", min_value=0.0, max_value=1.0,
                                         value=0.05, step=0.01,
                                         key="gsea_sig_threshold")
        sig_genes = stats_df[stats_df[sig_metric] < sig_threshold].index.tolist()
        st.caption(f"{len(sig_genes)} of {len(stats_df)} genes pass {sig_metric} < {sig_threshold}.")

        # All genes tested in this study — the recommended ORA/GSEA background, since it
        # reflects genes that could actually have been detected here, rather than assuming
        # the whole genome was tested.
        detected_genes = list(st.session_state.log2_data.index)

        method = st.selectbox(
            "Analysis type",
            ["Gene Set Enrichment Analysis (GSEA)", "Over-Representation Analysis (ORA)",
             "STRING Enrichment Analysis"],
            key="gsea_method"
        )

        # ===================================================================
        # Method 1: STRING Enrichment Analysis (requires internet)
        # ===================================================================
        if method == "STRING Enrichment Analysis":
            species_label = st.selectbox("Organism", list(enrichment.STRING_SPECIES.keys()),
                                          key="gsea_string_species")
            string_ver = enrichment.get_string_version()
            if string_ver:
                st.caption(f"STRING database version: {string_ver}")
            if st.button("Run STRING Enrichment Analysis", type="primary", key="gsea_run_string"):
                if len(sig_genes) < 1:
                    st.warning("No genes pass the current cutoff — relax it above.")
                else:
                    try:
                        with st.spinner("Querying STRING (string-db.org)..."):
                            result = enrichment.run_string_enrichment(
                                sig_genes, species=enrichment.STRING_SPECIES[species_label]
                            )
                        st.session_state.gsea_result = result
                        st.session_state["gsea_plot_on"] = st.session_state["gsea_dotplot_on"] = False
                        st.session_state.gsea_ranked_scores = None
                        st.session_state.gsea_method_used = f"STRING Enrichment Analysis ({species_label})"
                        st.session_state.gsea_run_meta = {
                            "Analysis method": "STRING Enrichment Analysis",
                            "Organism": species_label,
                            "Database": "STRING (own bundled GO/KEGG/Reactome/Pfam/InterPro annotation)",
                            "Database version": str(string_ver) if string_ver else "unavailable (could not query live)",
                            "# input genes": len(sig_genes),
                            "Significance cutoff": f"{sig_metric} < {sig_threshold}",
                        }
                        st.success(f"{len(result)} enriched terms returned.")
                    except enrichment.EnrichmentError as e:
                        st.error(str(e))

        # ===================================================================
        # Method 2: Over-Representation Analysis (ORA) — local hypergeometric test
        # ===================================================================
        elif method == "Over-Representation Analysis (ORA)":
            c1, c2 = st.columns(2)
            organism = c1.selectbox(
                "Organism", [o for o in enrichment.ORGANISM_CATALOG if "plant" not in o] +
                [o for o in enrichment.ORGANISM_CATALOG if "plant" in o],
                key="ora_organism"
            )
            org_cfg = enrichment.ORGANISM_CATALOG[organism]
            db_options = enrichment.ORA_GSEA_DATABASE_COLLECTIONS
            db_labels = {
                c: enrichment.display_collection_label(c) + ("" if org_cfg["databases"][c]["available"] else " — unavailable")
                for c in db_options
            }
            collection = c2.selectbox("Gene-set database", db_options, key="ora_collection",
                                       format_func=lambda c: db_labels[c])
            db_info = org_cfg["databases"][collection]
            if not db_info["available"]:
                st.warning(f"**{enrichment.display_collection_label(collection)}** is not available for {organism}: {db_info.get('note', '')}")
            else:
                st.caption(f"Database version: **{db_info['version']}**  ·  "
                           f"Gene identifier type: **{org_cfg['gene_id_type']}**")

            universe_choice = st.radio(
                "Background / universe",
                ["All detected genes in this study (recommended)", "Entire gene-set library"],
                key="ora_universe",
                help="The universe is the set of genes considered 'testable'. Using all genes "
                     "actually detected and retained after filtering in this study (rather than "
                     "assuming the whole genome) is the more defensible default, since it "
                     "reflects what could actually have been observed in this experiment."
            )
            c3, c4 = st.columns(2)
            min_size_ora = c3.number_input("Min gene set size", min_value=1, value=1, step=1, key="ora_min_size")
            max_size_ora = c4.number_input("Max gene set size (0 = no limit)", min_value=0, value=0, step=10,
                                            key="ora_max_size")

            if st.button("Run Over-Representation Analysis", type="primary", key="gsea_run_ora",
                         disabled=not db_info["available"]):
                if len(sig_genes) < 1:
                    st.warning("No genes pass the current cutoff — relax it above.")
                else:
                    try:
                        with st.spinner(f"Fetching {enrichment.display_collection_label(collection)} for {organism}..."):
                            gene_sets, resolved_library, resolved_version = enrichment.fetch_gene_set_library(
                                organism, collection
                            )
                        universe = detected_genes if universe_choice.startswith("All detected") else None
                        with st.spinner(f"Running hypergeometric ORA against {len(gene_sets)} gene sets..."):
                            result, run_meta = enrichment.run_ora_hypergeometric(
                                sig_genes, gene_sets, universe_genes=universe,
                                min_set_size=min_size_ora,
                                max_set_size=(max_size_ora if max_size_ora > 0 else None),
                            )
                        st.session_state.gsea_result = result
                        st.session_state["gsea_plot_on"] = st.session_state["gsea_dotplot_on"] = False
                        st.session_state.gsea_ranked_scores = None
                        st.session_state.gsea_method_used = f"Over-Representation Analysis ({organism}, {enrichment.display_collection_label(collection)})"
                        run_meta.update({
                            "Organism": organism, "Gene identifier type": org_cfg["gene_id_type"],
                            "Database name": resolved_library, "Database version": resolved_version,
                            "Significance cutoff": f"{sig_metric} < {sig_threshold}",
                        })
                        st.session_state.gsea_run_meta = run_meta
                        if len(result) == 0:
                            st.warning("No gene sets could be tested — none of the comparison genes were "
                                       "found in the chosen universe. See the metadata panel below for "
                                       "how many genes were mapped.")
                        else:
                            st.success(f"{len(result)} gene sets tested.")
                    except enrichment.EnrichmentError as e:
                        st.error(str(e))

        # ===================================================================
        # Method 3: Gene Set Enrichment Analysis (GSEA), preranked
        # ===================================================================
        else:
            st.markdown("**Gene set source**")
            gs_source = st.radio(
                "Gene sets from:",
                ["Organism gene-set database", "Upload custom .gmt file"],
                horizontal=True, key="gsea_prerank_source"
            )

            if gs_source == "Organism gene-set database":
                c1, c2 = st.columns(2)
                organism_g = c1.selectbox(
                    "Organism", [o for o in enrichment.ORGANISM_CATALOG if "plant" not in o] +
                    [o for o in enrichment.ORGANISM_CATALOG if "plant" in o],
                    key="gsea_organism"
                )
                org_cfg_g = enrichment.ORGANISM_CATALOG[organism_g]
                db_labels_g = {
                    c: enrichment.display_collection_label(c) + ("" if org_cfg_g["databases"][c]["available"] else " — unavailable")
                    for c in enrichment.ORA_GSEA_DATABASE_COLLECTIONS
                }
                collection_g = c2.selectbox("Gene-set database", enrichment.ORA_GSEA_DATABASE_COLLECTIONS,
                                             key="gsea_collection", format_func=lambda c: db_labels_g[c])
                db_info_g = org_cfg_g["databases"][collection_g]
                if not db_info_g["available"]:
                    st.warning(f"**{enrichment.display_collection_label(collection_g)}** is not available for {organism_g}: {db_info_g.get('note', '')}")
                else:
                    st.caption(f"Database version: **{db_info_g['version']}**  ·  "
                               f"Gene identifier type: **{org_cfg_g['gene_id_type']}**")
                if st.button("Fetch gene set library", key="gsea_fetch_library", disabled=not db_info_g["available"]):
                    try:
                        with st.spinner(f"Fetching {enrichment.display_collection_label(collection_g)} for {organism_g}..."):
                            gene_sets, resolved_library, resolved_version = enrichment.fetch_gene_set_library(
                                organism_g, collection_g
                            )
                        st.session_state.gsea_gene_sets = gene_sets
                        st.session_state.gsea_library_meta = {
                            "Organism": organism_g, "Gene identifier type": org_cfg_g["gene_id_type"],
                            "Database name": resolved_library, "Database version": resolved_version,
                        }
                        st.success(f"Fetched {len(gene_sets)} gene sets from '{resolved_library}' "
                                   f"(version: {resolved_version}).")
                    except enrichment.EnrichmentError as e:
                        st.error(str(e))
            else:
                gmt_file = st.file_uploader("Upload a .gmt gene set file", type=["gmt", "txt"], key="gsea_gmt_upload")
                if gmt_file is not None:
                    st.session_state.gsea_gene_sets = enrichment.parse_gmt_text(
                        gmt_file.getvalue().decode("utf-8", errors="ignore")
                    )
                    st.session_state.gsea_library_meta = {
                        "Organism": "user-uploaded (not organism-verified)",
                        "Gene identifier type": "as provided in uploaded file",
                        "Database name": gmt_file.name, "Database version": "user-supplied file, unversioned",
                    }
                    st.success(f"Parsed {len(st.session_state.gsea_gene_sets)} gene sets from the uploaded file.")

            c1, c2, c3 = st.columns(3)
            rank_method = c1.selectbox("Ranking metric", ["signed_neglogp", "log2fc", "statistic"],
                                        key="gsea_rank_method",
                                        help="signed_neglogp = sign(Log2FC) × -log10(p-value), the standard "
                                             "choice combining direction and significance. 'statistic' uses the "
                                             "Wald statistic column when pyDESeq2 produced the result.")
            min_size = c2.number_input("Min gene set size", min_value=2, value=15, step=1, key="gsea_min_size")
            max_size = c3.number_input("Max gene set size", min_value=5, value=500, step=5, key="gsea_max_size")
            n_perm = st.slider("Permutations (per gene set)", 100, 2000, 500, 100, key="gsea_n_perm",
                                help="More permutations give finer-grained p-values but take longer — "
                                     "500-1000 is typical for exploratory use.")

            if st.session_state.gsea_gene_sets is None:
                st.info("Load, fetch, or upload a gene set library above first.")
            elif st.button("Run Gene Set Enrichment Analysis (GSEA)", type="primary", key="gsea_run_prerank"):
                try:
                    ranked_scores = enrichment.compute_ranking_score(stats_df, method=rank_method)
                    with st.spinner(f"Running GSEA against {len(st.session_state.gsea_gene_sets)} "
                                     f"gene sets ({n_perm} permutations each — this can take a while)..."):
                        result, run_meta = enrichment.run_prerank_gsea(
                            ranked_scores, st.session_state.gsea_gene_sets,
                            n_perm=n_perm, min_size=min_size, max_size=max_size
                        )
                    st.session_state.gsea_result = result
                    st.session_state["gsea_plot_on"] = st.session_state["gsea_dotplot_on"] = False
                    st.session_state.gsea_ranked_scores = ranked_scores
                    st.session_state.gsea_method_used = f"Gene Set Enrichment Analysis (GSEA) ({rank_method})"
                    run_meta["Ranking metric"] = rank_method
                    run_meta.update(st.session_state.get("gsea_library_meta") or {})
                    st.session_state.gsea_run_meta = run_meta
                    st.success(f"{len(result)} gene sets tested (after size filtering).")
                except enrichment.EnrichmentError as e:
                    st.error(str(e))

        # ===================================================================
        # Shared results display
        # ===================================================================
        if st.session_state.gsea_result is not None:
            st.subheader(f"Results — {st.session_state.get('gsea_method_used', 'Enrichment')}")
            st.dataframe(utils.format_for_display(st.session_state.gsea_result), width='stretch', height=350)
            st.download_button(
                "Download GSEA_Results.csv", utils.to_download_bytes_csv(st.session_state.gsea_result),
                file_name="GSEA_Results.csv", mime="text/csv", key="dl_gsea_results"
            )

            if st.session_state.get("gsea_run_meta"):
                with st.expander("📋 Analysis metadata (for reproducibility)"):
                    meta_df = pd.DataFrame(
                        [{"Field": k, "Value": v} for k, v in st.session_state.gsea_run_meta.items()]
                    )
                    st.dataframe(meta_df, width='stretch', hide_index=True)
                    st.download_button(
                        "Download GSEA_Analysis_Metadata.csv", utils.to_download_bytes_csv(meta_df),
                        file_name="GSEA_Analysis_Metadata.csv", mime="text/csv", key="dl_gsea_meta"
                    )

            # ---- Dot plot of the top enriched pathways ----
            st.subheader("Pathway Dot Plot")
            st.caption(
                "One row per gene set: x-position is the NES (GSEA), Fold Enrichment (ORA), or Gene "
                "Ratio (STRING), dot size is the gene count, and dot color is the significance — the "
                "same visual style used across all three enrichment methods."
            )
            dc1, dc2, dc3, dc4 = st.columns(4)
            top_n_dot = dc1.selectbox("Show top", [5, 10, 20, 30, 50], index=1, key="gsea_dot_topn")
            has_nes = "NES" in st.session_state.gsea_result.columns
            selection_opts = ["Most significant"] + (["Top up + down"] if has_nes else [])
            dot_selection = dc2.selectbox("Selection", selection_opts, key="gsea_dot_selection")
            # STRING's results use its own lowercase TSV column names (fdr / p_value)
            # rather than the FDR / P-value columns GSEA and ORA produce; offer the
            # same canonical labels regardless of which method produced the results.
            _dotplot_col_alts = {"FDR": ("FDR", "fdr"), "P-value": ("P-value", "p_value")}
            color_opts = [label for label, alts in _dotplot_col_alts.items()
                          if any(a in st.session_state.gsea_result.columns for a in alts)]
            dot_color_by = dc3.selectbox("Color by", color_opts, key="gsea_dot_colorby")
            dot_cmap = dc4.selectbox("Colormap", ["coolwarm", "RdBu_r", "viridis", "plasma", "YlOrRd"],
                                      key="gsea_dot_cmap")
            dc5, dc6 = st.columns(2)
            dot_width = dc5.slider("Figure width (in)", 4.0, 16.0, 9.5, 0.5, key="gsea_dot_width")
            dot_font = dc6.slider("Font size", 6.0, 16.0, 10.0, 0.5, key="gsea_dot_font")

            try:
                fig_dot, plotted = enrichment.enrichment_dotplot(
                    st.session_state.gsea_result, top_n=top_n_dot, sort_by=dot_color_by,
                    selection=dot_selection, color_by=dot_color_by, cmap_name=dot_cmap,
                    fig_width_in=dot_width, font_size=dot_font,
                )
                utils.render_figure(st, fig_dot)
                utils.render_figure_download(st, fig_dot, "GSEA_DotPlot", key_prefix="gsea_dotplot")
                st.download_button(
                    "Download GSEA_DotPlot_Data.csv", utils.to_download_bytes_csv(plotted),
                    file_name="GSEA_DotPlot_Data.csv", mime="text/csv", key="dl_gsea_dotplot_data"
                )
            except enrichment.EnrichmentError as e:
                st.warning(str(e))

            if (st.session_state.gsea_ranked_scores is not None and st.session_state.gsea_gene_sets is not None
                    and "Term" in st.session_state.gsea_result.columns):
                st.subheader("Enrichment Plot")
                term_choice = st.selectbox("Gene set to plot", st.session_state.gsea_result["Term"].tolist(),
                                            key="gsea_plot_term")
                if st.button("Generate Enrichment Plot", key="gsea_plot_btn"):
                    st.session_state["gsea_plot_on"] = True
                if st.session_state.get("gsea_plot_on"):
                    try:
                        # Pull this term's statistics from the results table so they can be
                        # annotated directly onto the figure.
                        res_df = st.session_state.gsea_result
                        row = res_df[res_df["Term"] == term_choice]
                        def _get(col):
                            if row.empty or col not in row.columns:
                                return None
                            try:
                                return float(row.iloc[0][col])
                            except (TypeError, ValueError):
                                return None
                        fig_gsea = enrichment.running_score_plot(
                            st.session_state.gsea_ranked_scores, st.session_state.gsea_gene_sets, term_choice,
                            nes=_get("NES"), pval=_get("P-value"), fdr=_get("FDR"),
                        )
                        utils.render_figure(st, fig_gsea)
                        utils.render_figure_download(st, fig_gsea, f"GSEA_{term_choice[:40]}",
                                                      key_prefix="gsea_plot")
                    except enrichment.EnrichmentError as e:
                        st.error(str(e))
