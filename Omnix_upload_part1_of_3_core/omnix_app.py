"""
Omnix™ — Integrated Multi-Omics Data Analytics.

    streamlit run omnix_app.py

One site holding the three analysis platforms (apps/metabolomics, apps/proteomics,
apps/transcriptomics) behind the Omnix header, plus Subscription, Contact and Help pages.
Each platform also still runs on its own: streamlit run apps/<platform>/app.py
"""

import os
import sys

import streamlit as st

ROOT = os.path.dirname(os.path.abspath(__file__))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

_ICON = os.path.join(ROOT, "omnix_portal", "assets", "omnix_icon_192.png")
st.set_page_config(page_title="Omnix™ · Integrated Multi-Omics Data Analytics", layout="wide",
                   page_icon=_ICON if os.path.isfile(_ICON) else "🧬",
                   initial_sidebar_state="auto")

from omnix_portal import portal  # noqa: E402

portal.main()
