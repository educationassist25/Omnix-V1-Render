"""
The two apps, Metabolomics (metabolomics_app.py) and Lipidomics (lipidomics_app.py), run the same
pipeline as app.py and differ only in the analysis types offered.

    python -m pytest tests/test_split_apps.py -q
"""

import os

import pytest
from streamlit.testing.v1 import AppTest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
APPS = {
    "Metabolomics": (os.path.join(ROOT, "metabolomics_app.py"), ["Untargeted Metabolomics", "Targeted Metabolomics"]),
    "Lipidomics": (os.path.join(ROOT, "lipidomics_app.py"), ["Untargeted Lipidomics", "Targeted Lipidomics"]),
}
NAV = ["Setup", "Preprocessing", "Normalization", "Statistical Analysis", "Exploratory Analysis",
       "Biomarker Discovery", "Pathway Analysis", "Network Analysis"]


def _button(at, text):
    return next(b for b in at.button if text in b.label)


def _goto(at, group, page):
    at.session_state["nav_group"] = group
    at.session_state["nav_page"] = page
    return at.run()


@pytest.mark.parametrize("family", list(APPS))
def test_app_offers_only_its_two_analysis_types_and_same_pages(family):
    path, types = APPS[family]
    at = AppTest.from_file(path, default_timeout=120).run()
    assert not at.exception
    assert at.sidebar.title[0].value == family
    assert at.selectbox[0].label == "Analysis type" and at.selectbox[0].options == types
    assert at.session_state["data_type"] == types[0]
    # every navigation group of the original app is present, in the same order
    groups = [b.label for b in at.sidebar.button]
    assert groups == NAV
    # multi-dataset upload: each dataset's analysis type offers the same two types
    at.radio[0].set_value("Multiple datasets (combine methods/modes)").run()
    at.radio(key="multi_data_source").set_value("Real Data").run()
    assert at.selectbox(key="ds_type_0").options == types


def test_original_app_still_offers_all_four_types():
    at = AppTest.from_file(os.path.join(ROOT, "app.py"), default_timeout=120).run()
    assert not at.exception
    assert at.sidebar.title[0].value == "Metabolomics"
    assert at.selectbox[0].options == APPS["Metabolomics"][1][:1] + ["Targeted Metabolomics", "Untargeted Lipidomics",
                                                                     "Targeted Lipidomics"]


@pytest.mark.parametrize("family", list(APPS))
def test_single_dataset_targeted_and_untargeted_run_through_statistics(family):
    path, types = APPS[family]
    for dtype in types:
        at = AppTest.from_file(path, default_timeout=300).run()
        at.selectbox[0].set_value(dtype).run()
        _button(at, "Load & Validate Data").click().run()
        assert not at.exception, dtype
        assert at.session_state["raw_peak_df"] is not None and at.session_state["data_type"] == dtype
        _goto(at, "Normalization", "Normalization")
        assert not at.exception, dtype


@pytest.mark.parametrize("family", list(APPS))
def test_multi_dataset_demo_through_full_pipeline(family):
    path, types = APPS[family]
    at = AppTest.from_file(path, default_timeout=300).run()
    at.radio[0].set_value("Multiple datasets (combine methods/modes)").run()
    at.radio(key="multi_data_source").set_value("Demo Data").run()
    _button(at, "Load Demo Datasets").click().run()
    assert not at.exception
    assert list(at.session_state["datasets"]) == types           # only this app's two demo datasets
    _button(at, "Auto-Process All Datasets").click().run()
    assert not at.exception
    _goto(at, "Normalization", "Combined Normalization Data (Multi-Dataset)")
    _button(at, "Generate Combined Normalized Data").click().run()
    assert not at.exception
    combined = at.session_state["log2_data"]
    assert combined is not None and len(combined) > 50
    prefixes = {str(i).split("::", 1)[0] for i in combined.index}
    assert prefixes == set(types)
    for group, page in [("Statistical Analysis", "Statistics"), ("Exploratory Analysis", "PCA"),
                        ("Exploratory Analysis", "Heatmap"), ("Pathway Analysis", "Metabolite-Set Enrichment Analysis (MSEA)"),
                        ("Network Analysis", "Correlation Rewiring Map")]:
        _goto(at, group, page)
        assert not at.exception, page
