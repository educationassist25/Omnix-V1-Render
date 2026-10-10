"""
Full-pipeline integration test: drive the real MetaboAI Pro UI from data loading to the
Rewiring Map, with no hand-built session state.

    Data Upload (Multiple datasets -> Demo Data -> Load Demo Datasets)
      -> Auto-Process All Datasets
      -> Combined Normalization Data -> Generate Combined Normalized Data
      -> Network Analysis -> Correlation Rewiring Map -> Run

    python -m pytest tests/test_rewiring_pipeline.py -q
"""

import os

import numpy as np
from streamlit.testing.v1 import AppTest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MAIN = os.path.join(ROOT, "app.py")


def _button(at, text):
    return next(b for b in at.button if text in b.label)


def _goto(at, group, page):
    at.session_state["nav_group"] = group
    at.session_state["nav_page"] = page
    return at.run()


def test_multi_dataset_demo_through_full_pipeline_to_rewiring():
    at = AppTest.from_file(MAIN, default_timeout=300).run()
    assert not at.exception

    at.radio[0].set_value("Multiple datasets (combine methods/modes)").run()
    at.radio(key="multi_data_source").set_value("Demo Data").run()
    _button(at, "Load Demo Datasets").click().run()
    assert not at.exception
    assert len(at.session_state["datasets"]) == 4

    _button(at, "Auto-Process All Datasets").click().run()
    assert not at.exception

    _goto(at, "Normalization", "Combined Normalization Data (Multi-Dataset)")
    _button(at, "Generate Combined Normalized Data").click().run()
    assert not at.exception
    combined = at.session_state["log2_data"]
    assert combined is not None and combined.shape[0] > 100
    assert any("::" in str(i) for i in combined.index)       # multi-dataset feature names
    before = combined.copy()

    _goto(at, "Network Analysis", "Correlation Rewiring Map")
    assert not at.exception
    assert "combined normalized data" in " ".join(str(c.value) for c in at.caption)
    # 24 samples over 4 groups -> 6 per group: a note about small groups (not a warning)
    assert "6 samples in the smaller group" in " ".join(str(c.value) for c in at.caption)
    # the multimethod demo ships Method/Pathway-style row annotations; they must be offered
    pw_opts = at.selectbox(key="rewire_main_pw_col").options
    assert len(pw_opts) >= 2
    at.number_input(key="rewire_main_nperm").set_value(100)
    at.number_input(key="rewire_main_nboot").set_value(50)
    _button(at, "Run Rewiring Analysis").click().run()
    assert not at.exception
    stored = at.session_state["rewire_main_result"]
    assert stored is not None
    res = stored["result"]
    # combined data is already log2 -> must not be transformed again, and names are display names
    assert res.prep.logged is False
    assert not any("::" in f for f in res.prep.features)
    assert np.isfinite(res.pairs["z"]).all()
    assert at.get("iframe")
    # the page must not modify pipeline state
    import pandas as pd
    pd.testing.assert_frame_equal(at.session_state["log2_data"], before)
