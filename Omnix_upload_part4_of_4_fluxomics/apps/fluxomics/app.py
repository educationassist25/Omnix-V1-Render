"""
Omnix™ Fluxomics — isotope-tracing analysis for unbiased LC-MS metabolomics.
Run:  streamlit run app.py

Inside the Omnix portal (omnix_app.py) this file runs as the Fluxomics platform. Its navigation then
follows the other Omnix platforms: the SIDEBAR lists the workflow groups, the MAIN PAGE shows the selected
group's steps as a row of buttons, and the selected step renders below it. The steps themselves (ui/,
core/) are unchanged.
"""
import re
import sys
from pathlib import Path

import streamlit as st

_HERE = str(Path(__file__).resolve().parent)
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

# Inside the Omnix portal the portal sets the page title and header; run on its own, the app sets its own.
OMNIX = globals().get("OMNIX", False)
if not OMNIX:
    st.set_page_config(page_title="Omnix Fluxomics", page_icon="🧪", layout="wide", initial_sidebar_state="expanded")

from core import __version__  # noqa: E402
from ui import common as C  # noqa: E402
from ui import (step1_upload, step2_design, step4_normalize, step5_incorporation,  # noqa: E402
                step6_boxplots, step7_stats, step8_heatmap, step9_pathways, step11_kinetics,
                step13_export)

APP_NAME = "Fluxomics"

CSS = """
<style>
h1,h2,h3{letter-spacing:-0.01em}
.step-head{border-bottom:1px solid #e6e5ee;padding-bottom:0.5rem;margin-bottom:1rem}
.step-kicker{font-size:0.78rem;letter-spacing:.08em;text-transform:uppercase;color:#000;font-weight:700}
.step-title{font-size:1.45rem;font-weight:700;margin:0.1rem 0}
.step-desc{color:#000;max-width:72ch;margin:0}
.direction{background:#f4f2fc;border:1px solid #dcd6f5;border-radius:8px;padding:8px 12px;margin:6px 0 12px;font-size:0.92rem}
.direction .pos{color:#000;font-weight:700}
.note{background:#fff8e6;border-left:3px solid #e0a400;padding:8px 12px;border-radius:4px;font-size:0.9rem;margin:6px 0}
.okbox{background:#effaf3;border-left:3px solid #1f9d55;padding:8px 12px;border-radius:4px;font-size:0.9rem;margin:6px 0}
.ds-card{background:#ffffff;border:1px solid #bfe3dc;border-radius:8px;padding:8px 10px;font-size:0.82rem;color:#000;margin:12px 0}
.st-key-intbox{background:#ffffff !important;border:1.5px solid #00897b !important;border-radius:8px;padding:8px 10px;margin:4px 0 12px}
.st-key-intbox label p{font-weight:600}
.int-head{font-size:0.8rem;font-weight:700;color:#000;line-height:1.3;padding-bottom:12px;display:block}
[data-testid="stMetricValue"]{font-variant-numeric:tabular-nums}
/* v2.4.2: all text solid black (captions, labels, help text, metrics, tabs, expanders, notes) */
.stApp, .stApp p, .stApp li, .stApp label, .stApp span, .stApp h1, .stApp h2, .stApp h3, .stApp h4, .stApp h5,
[data-testid="stCaptionContainer"], [data-testid="stCaptionContainer"] *, .stCaption, small,
[data-testid="stWidgetLabel"], [data-testid="stWidgetLabel"] *, [data-testid="stMetricLabel"], [data-testid="stMetricLabel"] *,
[data-testid="stMetricValue"], [data-testid="stExpander"] summary *, [data-baseweb="tab"] p, .note, .okbox, .direction,
[data-testid="stMarkdownContainer"], [data-testid="stMarkdownContainer"] *:not(code),
[data-testid="stTooltipContent"], [data-baseweb="select"] *, input, textarea {color:#000000 !important; opacity:1 !important}
.stButton button[kind="primary"] *, .stDownloadButton button[kind="primary"] *, .stFormSubmitButton button[kind="primary"] * {color:#ffffff !important}
[data-testid="stAlertContentError"] *, [data-testid="stAlertContentWarning"] *, [data-testid="stAlertContentInfo"] *,
[data-testid="stAlertContentSuccess"] * {color:#000000 !important}

/* ===== Omnix platform shell (same as Omnix Metabolomics): sidebar = workflow groups, main = step row ===== */
[data-testid="stSidebarContent"] { display: flex; flex-direction: column; }
[data-testid="stSidebarHeader"] { order: 0; }
[data-testid="stSidebarUserContent"] { order: 1; padding-bottom: 0px !important; }
[data-testid="stSidebar"] { background: linear-gradient(180deg, #f8fafc 0%, #e6f4f2 100%) !important;
  border-right: 1px solid #e2e8f0; }
[data-testid="stSidebar"] h1 { color: #00695c !important; }

div[data-testid="stVerticalBlock"][class*="st-key-nav_group_box"] { background: #e6f4f2; border: 1px solid #bfe3dc;
  border-radius: 10px; padding: 8px !important; gap: 4px !important; }
[class*="st-key-nav_group_box"] [data-testid="stButton"] button { background: transparent !important;
  border: 1px solid transparent !important; border-left: 3px solid #bfe3dc !important; box-shadow: none !important;
  border-radius: 6px !important; justify-content: flex-start !important; text-align: left !important;
  white-space: nowrap !important; overflow: visible !important; padding: 10px 12px !important; width: 100%;
  transform: none !important; min-height: 0 !important;
  transition: background-color 0.12s ease-in-out, color 0.12s ease-in-out, border-left-color 0.12s ease-in-out; }
[class*="st-key-nav_group_box"] [data-testid="stButton"] button p { color: #00695c !important; font-weight: 700 !important;
  font-size: 18px !important; line-height: 1.15 !important; letter-spacing: 0.1px; white-space: nowrap !important;
  word-break: keep-all !important; overflow-wrap: normal !important; hyphens: none !important; }
[class*="st-key-nav_group_box"] [data-testid="stButton"] button:hover { background: rgba(0, 105, 92, 0.10) !important;
  border-left-color: #00897b !important; }
[class*="st-key-nav_group_box"] [data-testid="stButton"] button:hover p,
[class*="st-key-nav_group_box"] [data-testid="stButton"] button:focus:not(:active) p { color: #00695c !important; }

[class*="st-key-nav_page_row"] [data-testid="stHorizontalBlock"] { display: flex; flex-wrap: wrap; row-gap: 8px; column-gap: 10px; }
[class*="st-key-nav_page_row"] [data-testid="stColumn"] { width: auto !important; flex: 0 0 auto !important; min-width: 0 !important; }
[class*="st-key-nav_page_row"] [data-testid="stButton"] button { background: #f8fafc !important; border: 1px solid #e2e8f0 !important;
  box-shadow: none !important; border-radius: 8px !important; white-space: nowrap !important; overflow: visible !important;
  padding: 10px 22px !important; transform: none !important;
  transition: background-color 0.12s ease-in-out, border-color 0.12s ease-in-out; }
[class*="st-key-nav_page_row"] [data-testid="stButton"] button p { color: #334155 !important; font-weight: 600 !important;
  font-size: 22px !important; line-height: 1.15 !important; white-space: nowrap !important; word-break: keep-all !important;
  overflow-wrap: normal !important; hyphens: none !important; }
[class*="st-key-nav_page_row"] [data-testid="stButton"] button:hover { background: #eef2f7 !important; border-color: #cbd5e1 !important; }
[class*="st-key-nav_page_row"] [data-testid="stButton"] button:focus:not(:active) p { color: #334155 !important; }

/* action buttons teal, download buttons orange (as in the other Omnix platforms) */
[data-testid="stButton"] button { background: linear-gradient(135deg, #00897b 0%, #00695c 100%) !important; border: none !important;
  border-radius: 8px !important; padding: 0.5rem 1.1rem !important; box-shadow: 0 2px 6px rgba(0, 105, 92, 0.35) !important;
  transition: transform 0.12s ease-in-out, box-shadow 0.12s ease-in-out, background 0.12s ease-in-out !important;
  height: auto !important; min-height: 2.5rem !important; white-space: normal !important; }
[data-testid="stButton"] button p { color: #ffffff !important; font-weight: 600 !important; white-space: normal !important;
  overflow-wrap: break-word !important; word-break: break-word !important; line-height: 1.25 !important; }
[data-testid="stButton"] button:hover { background: linear-gradient(135deg, #00695c 0%, #004d40 100%) !important;
  box-shadow: 0 4px 12px rgba(0, 105, 92, 0.45) !important; transform: translateY(-1px) !important; }
[data-testid="stButton"] button:hover p, [data-testid="stButton"] button:focus:not(:active) p { color: #ffffff !important; }
[data-testid="stButton"] button:active { transform: translateY(0) !important; }
[data-testid="stButton"] button:disabled { background: #cbd5e1 !important; box-shadow: none !important; transform: none !important; }
[data-testid="stButton"] button:disabled p { color: #64748b !important; }
[data-testid="stDownloadButton"] button { background: linear-gradient(135deg, #f57c00 0%, #e65100 100%) !important; border: none !important;
  border-radius: 8px !important; padding: 0.5rem 1.1rem !important; box-shadow: 0 2px 6px rgba(230, 81, 0, 0.35) !important;
  height: auto !important; min-height: 2.5rem !important; white-space: normal !important; }
[data-testid="stDownloadButton"] button p { color: #ffffff !important; font-weight: 600 !important; white-space: normal !important;
  line-height: 1.25 !important; }
[data-testid="stDownloadButton"] button:hover { background: linear-gradient(135deg, #e65100 0%, #bf360c 100%) !important;
  box-shadow: 0 4px 12px rgba(230, 81, 0, 0.45) !important; transform: translateY(-1px) !important; }
[data-testid="stDownloadButton"] button:hover p { color: #ffffff !important; }
</style>
"""
st.markdown(CSS, unsafe_allow_html=True)
if not OMNIX:   # on its own: room for Streamlit's toolbar and a comfortable reading width
    st.markdown("<style>.block-container{padding-top:3.4rem;max-width:1400px}</style>", unsafe_allow_html=True)
C.init_state()

PAGES = {"upload": step1_upload, "design": step2_design, "normalize": step4_normalize,
         "incorporation": step5_incorporation, "boxplots": step6_boxplots, "stats": step7_stats,
         "pathways": step9_pathways, "heatmap": step8_heatmap,
         "kinetics": step11_kinetics, "export": step13_export}

# ===========================================================================
# Grouped navigation (as Omnix Metabolomics). The workflow steps and their numbers are unchanged;
# st.session_state.step stays the one source of truth, so every "Continue to …" / "Go to …" button
# inside the steps keeps working and simply moves the sidebar group along.
# ===========================================================================
NAV_GROUPS = [
    ("Setup", ["upload", "design"]),
    ("Normalization", ["normalize"]),
    ("Isotope Labeling", ["incorporation", "kinetics"]),
    ("Statistical Analysis", ["stats"]),
    ("Exploratory Analysis", ["boxplots", "heatmap"]),
    ("Pathway Analysis", ["pathways"]),
    ("Export & Report", ["export"]),
]
GROUP_OF = {step: name for name, steps in NAV_GROUPS for step in steps}
STEP_LABEL = dict(C.STEPS)


def _slug(text):
    """CSS-class-safe key fragment for a group/step label."""
    return re.sub(r"[^a-z0-9]+", "-", text.lower()).strip("-")


if st.session_state.step not in GROUP_OF:
    st.session_state.step = NAV_GROUPS[0][1][0]
active_group = GROUP_OF[st.session_state.step]

# --- Sidebar: group names only ---
st.sidebar.title(APP_NAME)
st.sidebar.caption("Stable-isotope tracing metabolomics & flux analysis (LC-MS)")
with st.sidebar.container(key="nav_group_box"):
    for group_name, steps in NAV_GROUPS:
        if st.button(group_name, key=f"navgrp__{_slug(group_name)}", use_container_width=True):
            st.session_state.step = steps[0]
            st.rerun()

ds = st.session_state.ds
if ds is not None:
    n_lab = int((ds.S.SampleType == "13C labeled").sum())
    st.sidebar.markdown(f"<div class='ds-card'><b>Dataset</b><br>{'Compound Discoverer' if ds.source == 'cd' else 'Generic table'}"
                        f"<br>{ds.CT.Compound.nunique()} metabolites · {len(ds.CT)} isotopologue channels<br>"
                        f"{len(ds.S)} samples ({n_lab} labeled)</div>", unsafe_allow_html=True)
    with st.sidebar:
        try:
            box = st.container(border=True, key="intbox")       # Streamlit ≥ 1.39: keyed, styled via CSS
        except TypeError:
            box = st.container()
        with box:
            st.markdown("<div class='int-head'>Intensity dataset for downstream steps</div>",
                        unsafe_allow_html=True)
            st.radio("Intensity dataset for downstream steps", ["log2", "log2_iqr"], key="intensity_choice",
                     format_func={"log2_iqr": "Log2 + IQR-normalized intensity", "log2": "Log2 intensity"}.get,
                     label_visibility="collapsed")
with st.sidebar.expander("Session"):
    if st.button("Reset all steps", key="reset"):
        for k in list(st.session_state.keys()):
            if not str(k).startswith("omnix_"):        # inside Omnix: keep the portal's own state (page, access)
                del st.session_state[k]
        st.rerun()
    st.caption(f"Fluxomics version {__version__}")

# Highlight the active group: solid teal pill with an orange edge.
_g = _slug(active_group)
st.markdown(
    f"""<style>
    [class*="st-key-nav_group_box"] [class*="st-key-navgrp__{_g}"] [data-testid="stButton"] button,
    [class*="st-key-nav_group_box"] [class*="st-key-navgrp__{_g}"] [data-testid="stButton"] button:focus:not(:active) {{
        background: #00695c !important; border-left-color: #e65100 !important; }}
    [class*="st-key-nav_group_box"] [class*="st-key-navgrp__{_g}"] [data-testid="stButton"] button p,
    [class*="st-key-nav_group_box"] [class*="st-key-navgrp__{_g}"] [data-testid="stButton"] button:focus:not(:active) p {{
        color: #ffffff !important; }}
    [class*="st-key-nav_group_box"] [class*="st-key-navgrp__{_g}"] [data-testid="stButton"] button:hover {{
        background: #004d40 !important; }}
    [class*="st-key-nav_group_box"] [class*="st-key-navgrp__{_g}"] [data-testid="stButton"] button:hover p {{
        color: #ffffff !important; }}
    </style>""",
    unsafe_allow_html=True,
)

# --- Main page: the active group's steps, as a row of buttons (✓ = step completed) ---
steps = dict(NAV_GROUPS)[active_group]
with st.container(key="nav_page_row"):
    cols = st.columns(len(steps))
    for col, key in zip(cols, steps):
        with col:
            label = ("✓ " if C.done(key) else "") + STEP_LABEL[key]
            if st.button(label, key=f"navpg__{_g}__{key}", use_container_width=False):
                st.session_state.step = key
                st.rerun()

# The selected step: orange, as the page tabs of the other platforms.
_p = f"navpg__{_g}__{st.session_state.step}"
st.markdown(
    f"""<style>
    [class*="st-key-nav_page_row"] [class*="st-key-{_p}"] [data-testid="stButton"] button,
    [class*="st-key-nav_page_row"] [class*="st-key-{_p}"] [data-testid="stButton"] button:focus:not(:active) {{
        background: #e65100 !important; border-color: #e65100 !important; }}
    [class*="st-key-nav_page_row"] [class*="st-key-{_p}"] [data-testid="stButton"] button p,
    [class*="st-key-nav_page_row"] [class*="st-key-{_p}"] [data-testid="stButton"] button:focus:not(:active) p {{
        color: #ffffff !important; font-weight: 700 !important; }}
    [class*="st-key-nav_page_row"] [class*="st-key-{_p}"] [data-testid="stButton"] button:hover {{
        background: #bf360c !important; border-color: #bf360c !important; }}
    [class*="st-key-nav_page_row"] [class*="st-key-{_p}"] [data-testid="stButton"] button:hover p {{
        color: #ffffff !important; }}
    </style>""",
    unsafe_allow_html=True,
)

st.divider()

PAGES[st.session_state.step].render()
