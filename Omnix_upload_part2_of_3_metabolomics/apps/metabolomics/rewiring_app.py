"""
Omnix Metabolomics — Correlation Rewiring Map (standalone app)

    streamlit run rewiring_app.py

Finds metabolite pairs whose relationship changes between two groups (couplings
lost, gained or flipped) and turns them into ranked, testable hypotheses. The same
analysis is also available inside the main app (app.py -> Network Analysis), where it
runs on the pipeline's normalized log2 data.
"""

import os
import sys

import pandas as pd
import streamlit as st

sys.path.insert(0, os.path.dirname(__file__))
from metabolomics_modules import normalization, rewiring_ui, utils  # noqa: E402

HERE = os.path.dirname(__file__)

st.set_page_config(page_title="Rewiring Map · Omnix Metabolomics", layout="wide", page_icon="🕸️")
st.markdown(
    """
    <style>
    [data-testid="stSidebar"] { background: linear-gradient(180deg, #f8fafc 0%, #e6f4f2 100%) !important;
                                border-right: 1px solid #e2e8f0; }
    [data-testid="stSidebar"] h1 { color: #00695c !important; }
    [data-testid="stButton"] button[kind="primary"] { background: #00695c !important; border: none !important; }
    [data-testid="stButton"] button[kind="primary"]:hover { background: #004d40 !important; }
    </style>
    """,
    unsafe_allow_html=True,
)

SOURCES = [
    rewiring_ui.DEMO_LABEL,
    "Omnix Metabolomics demo: Rich Clinical (4 groups × 9 samples)",
    "Upload my data",
]

st.sidebar.title("Rewiring Map")
st.sidebar.caption("Omnix Metabolomics · differential correlation")
source = st.sidebar.radio("Data source", SOURCES, key="rw_source")


def _load_upload():
    st.sidebar.markdown("**Peak area matrix** (rows = metabolites, columns = samples)")
    peak_file = st.sidebar.file_uploader("Peak area matrix", type=["csv", "xlsx", "xls"], key="rw_peak",
                                         label_visibility="collapsed")
    st.sidebar.markdown("**Sample metadata** (`Sample`, `Group`, optional `IsQC`, `Batch`, …)")
    meta_file = st.sidebar.file_uploader("Metadata", type=["csv", "xlsx", "xls"], key="rw_meta",
                                         label_visibility="collapsed")
    st.sidebar.markdown("**Metabolite annotations** (optional: `Metabolite`, `Pathway`, …)")
    ann_file = st.sidebar.file_uploader("Row annotations", type=["csv", "xlsx", "xls"], key="rw_ann",
                                        label_visibility="collapsed")
    if peak_file is None or meta_file is None:
        return None
    try:
        peak = utils.validate_peak_matrix(utils.load_table(peak_file))
        meta = utils.validate_metadata(utils.load_table(meta_file), peak.columns)
        ann = None
        if ann_file is not None:
            raw = utils.load_table(ann_file)
            first = raw.columns[0]
            ann = raw.rename(columns={first: "Metabolite"}).astype({"Metabolite": str}).set_index("Metabolite")
    except utils.DataValidationError as err:
        st.error(f"Couldn't use these files: {err}")
        return None
    except Exception as err:  # malformed spreadsheets, etc.
        st.error(f"Couldn't read these files: {err}")
        return None
    return peak, meta, ann, ", ".join(f.name for f in (peak_file, meta_file) if f is not None)


if source == SOURCES[0]:
    peak, meta, ann = rewiring_ui.load_simulated_demo()
    label = "Simulated demo"
    blurb = ("43 targeted metabolites with planted rewiring: glutamine decoupled from the TCA cycle, an "
             "IDO1-like tryptophan–kynurenine sign flip, a lactate–kynurenine coupling gained in Tumor, two balanced "
             "batches with a batch shift, and three outlier samples that create a decoy correlation.")
    with st.sidebar.expander("Use this demo in the main Omnix Metabolomics app"):
        st.caption("Download the files and load them on the Data Upload page (Targeted Metabolomics).")
        for name, data in rewiring_ui.demo_csv_bytes().items():
            st.download_button(name, data, file_name=name, mime="text/csv", key=f"rw_demo_{name}")
    median_default = False
elif source == SOURCES[1]:
    peak = utils.validate_peak_matrix(pd.read_csv(os.path.join(HERE, "sample_peak_area_matrix_richdemo.csv")))
    meta = utils.validate_metadata(pd.read_csv(os.path.join(HERE, "sample_metadata_richdemo.csv")), peak.columns)
    ann = pd.read_csv(os.path.join(HERE, "metabolite_row_annotations_richdemo.csv")).set_index("Metabolite")
    istds = normalization.detect_istds(peak.index, ann)          # internal standards are not biology
    peak, ann = normalization.drop_istds(peak, istds), ann.drop(index=[i for i in istds if i in ann.index])
    label = "Omnix Metabolomics Rich Clinical demo"
    blurb = ("The main app's Rich Clinical demo (200 metabolites, 9 samples per diagnosis group). With 9 samples "
             "per group only very large correlation changes can be detected, which makes it a useful example of "
             "what the small-sample warnings mean.")
    median_default = True
else:
    loaded = _load_upload()
    if loaded is None:
        st.title("Correlation Rewiring Map")
        st.info("Upload a peak area matrix and a sample metadata table in the sidebar. The formats are the same as "
                "Omnix Metabolomics' Data Upload page.")
        st.stop()
    peak, meta, ann, label = loaded
    blurb = f"{peak.shape[0]} features × {peak.shape[1]} samples from {label}."
    median_default = True

st.title("Correlation Rewiring Map")
st.markdown(
    "Finds metabolite pairs whose **relationship** changes between two groups, even when their levels don't. "
    "Edges show couplings that are lost, gained or flipped, linked to the enzymatic reactions that connect them "
    "(Human-GEM); pathways are tested as sets with the same libraries as MSEA; the hypothesis cards turn it all into "
    "ranked, testable ideas."
)
st.caption(blurb)

rewiring_ui.render_workspace(
    st, peak, meta, ann, key_prefix=f"rwapp_{SOURCES.index(source)}", source_label=label,
    data_is_log2=False, default_median_center=median_default,
)
