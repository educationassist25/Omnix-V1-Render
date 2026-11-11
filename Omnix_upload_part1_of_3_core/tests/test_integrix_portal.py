"""
Integrix inside the Omnix portal: registration, plans, the module -> Integrix bridge and the session hand-over.

    python -m pytest tests/test_integrix_portal.py -q
"""
import hashlib
import io
import os
import sys
import warnings
import zipfile

import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest

warnings.filterwarnings("ignore")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PORTAL = os.path.join(ROOT, "omnix_app.py")
NORM = os.path.join(ROOT, "apps", "integrix", "demo_data", "normalized")
sys.path.insert(0, ROOT)

from omnix_portal import access, billing, config, integrix_bridge  # noqa: E402
from omnix_portal import portal  # noqa: E402

TEST_KEY = "OMX-TEST-KEYS-1234-ABCD"
TEST_HASH = hashlib.sha256(TEST_KEY.encode()).hexdigest()


@pytest.fixture(autouse=True)
def _reset_guess_counter():
    access._fails.clear()
    yield
    access._fails.clear()


def _portal(secrets=None, timeout=300):
    at = AppTest.from_file(PORTAL, default_timeout=timeout)
    if secrets is not None:
        at.secrets["omnix"] = secrets
    return at.run()


def _open(at, page):
    at.session_state["omnix_page"] = page
    at.run()
    if at.session_state["omnix_page"] != page:   # Streamlit 1.52 test tool: a value set right after a click can be lost
        at.session_state["omnix_page"] = page
        at.run()
    return at


def _activate(at, key):
    _open(at, "subscription")
    at.text_input(key="omnix_key_input").input(key)
    next(b for b in at.button if b.label == "Activate").click().run()
    return at


def _module_state(layer):
    """What a module keeps in the session after Normalization, built from Integrix's normalized demo tables."""
    df = pd.read_csv(os.path.join(NORM, f"{layer}.csv"))
    meta = pd.read_csv(os.path.join(NORM, "metadata.csv")).set_index("Sample")
    meta["IsQC"] = False
    state = {"meta": meta, "processing_notes": ["Missingness filter 50%", "kNN imputation (k=5)",
                                                 "Median normalization, log2"]}
    if layer == "metabolomics":
        X = df.set_index("Metabolite").drop(columns=["HMDB"])
        state["row_annotations"] = df.set_index("Metabolite")[["HMDB"]]
    elif layer == "proteomics":
        X = df.set_index("Protein").drop(columns=["Gene"])
        state["gene_map"] = dict(zip(df["Protein"], df["Gene"]))
    else:
        X = df.set_index("Gene")
    state["log2_data"] = X
    return state


# ---------------------------------------------------------------------------------------------- registration
def test_integrix_is_registered_and_found():
    assert "integrix" in config.PLATFORMS and config.PLATFORMS["integrix"]["integration"]
    assert config.INTEGRATION_PLATFORMS == ["integrix"]
    assert "integrix" not in config.DATA_PLATFORMS and len(config.DATA_PLATFORMS) == 5
    assert portal.find_app("integrix") and portal.find_app("integrix").endswith(os.path.join("integrix", "app.py"))
    assert "integrix" in portal.INTEGRATION_APPS and "integrix" not in portal.DATA_APPS


def test_plans_cover_integrix_only_with_all_five_platforms():
    assert "integrix" in access._platforms({})                                       # no list: everything
    assert "integrix" in access._platforms({"platforms": "all"})
    assert "integrix" in access._platforms({"platforms": list(config.DATA_PLATFORMS)})
    assert "integrix" not in access._platforms({"platforms": ["proteomics"]})
    assert "integrix" not in access._platforms({"platforms": ["proteomics", "metabolomics"]})
    assert "integrix" not in billing.CHOOSABLE                                       # not a single/duo choice


# ---------------------------------------------------------------------------------------------- bridge
def test_bridge_export_round_trips_through_integrix_reader():
    sys.path.insert(0, os.path.join(ROOT, "apps", "integrix"))
    from trace_core import omnix_io
    files = []
    for layer in ("transcriptomics", "proteomics", "metabolomics"):
        st_ = _module_state(layer)
        st_["meta"].loc["S01", "IsQC"] = True                                       # a QC injection is dropped
        matrix, meta, info = integrix_bridge.build(layer, st_)
        assert "S01" not in matrix.columns and "S01" not in set(meta["Sample"]) and "IsQC" not in meta.columns
        assert matrix.shape[1] - len(info["id_columns"]) == 39
        assert info["imputation"].startswith("kNN") and "Median normalization" in info["normalization"]
        files.append((f"{layer}_integrix_export.zip", integrix_bridge.to_zip(matrix, meta, info)))
    _, _, info_p = integrix_bridge.build("proteomics", _module_state("proteomics"))
    assert info_p["id_columns"] == ["Protein", "Gene"]
    _, _, info_m = integrix_bridge.build("metabolomics", _module_state("metabolomics"))
    assert info_m["id_columns"] == ["Metabolite", "HMDB"]
    raw, meta, infos, msgs = omnix_io.read_exports(files)
    assert set(raw) == {"rna", "prot", "met"} and meta is not None and len(meta) == 39
    assert all(i["module"].startswith("Omnix ") for i in infos.values())
    assert not [m for m in msgs if "raw" in m.lower()], msgs


def test_bridge_needs_three_non_qc_samples():
    st_ = _module_state("proteomics")
    st_["meta"]["IsQC"] = True
    st_["meta"].iloc[:2, st_["meta"].columns.get_loc("IsQC")] = False
    assert integrix_bridge.build("proteomics", st_) is None
    assert integrix_bridge.build("proteomics", {"meta": st_["meta"]}) is None


# ---------------------------------------------------------------------------------------------- portal pages
@pytest.mark.parametrize("page", ["integrix", "help-integrix", "data-preparation"])
def test_integrix_deep_links(page):
    at = AppTest.from_file(PORTAL, default_timeout=300)
    at.query_params["page"] = page
    at.run()
    assert not at.exception and at.session_state["omnix_page"] == page


def test_home_lists_integrix_and_demo_runs():
    at = _portal()
    labels = [b.label for b in at.button]
    assert "Open Integrix  →" in labels and "Integrix" in labels
    _open(at, "integrix")
    assert not at.exception
    assert "Run analysis" in [b.label for b in at.sidebar.button]
    assert "Subscribe or enter access key" in [b.label for b in at.button]           # demo banner
    next(b for b in at.sidebar.button if b.label == "Run analysis").click().run(timeout=600)
    assert not at.exception
    assert "trace" in at.session_state


def test_demo_users_cannot_use_own_data_in_integrix():
    at = _portal({"access_keys": {TEST_HASH: {"name": "Lab", "plan": "Standard"}}})
    _open(at, "integrix")
    at.radio(key="src").set_value("Upload my data").run()
    assert not at.exception
    ups = at.get("file_uploader")
    assert ups and all(u.disabled for u in ups)


def test_single_omics_key_keeps_integrix_in_demo():
    at = _portal({"access_keys": {TEST_HASH: {"name": "Sam", "plan": "Single-Omics", "platforms": ["proteomics"]}}})
    _activate(at, TEST_KEY)
    _open(at, "integrix")
    assert "Upgrade your plan" in [b.label for b in at.button]


def test_multi_omics_key_unlocks_integrix_uploads():
    at = _portal({"access_keys": {TEST_HASH: {"name": "Lab", "plan": "Multi-Omics", "expires": "2099-12-31"}}})
    _activate(at, TEST_KEY)
    _open(at, "integrix")
    at.radio(key="src").set_value("Upload my data").run()
    ups = at.get("file_uploader")
    assert ups and not any(u.disabled for u in ups)


def test_sent_data_is_used_from_this_omnix_session():
    at = _portal({"access_keys": {TEST_HASH: {"name": "Lab", "plan": "Multi-Omics", "expires": "2099-12-31"}}},
                 timeout=900)
    _activate(at, TEST_KEY)
    inbox = {}
    for layer in ("transcriptomics", "proteomics", "metabolomics"):
        m, md, info = integrix_bridge.build(layer, _module_state(layer))
        inbox[layer] = {"platform_name": layer.capitalize(), "layer": layer, "zip": integrix_bridge.to_zip(m, md, info),
                        "zip_name": f"{layer}_integrix_export.zip", "n_features": len(m), "n_samples": len(md),
                        "sent": "10:00", "default_use": True}
    at.session_state[integrix_bridge.INBOX] = inbox
    _open(at, "integrix")
    assert not at.exception
    assert at.session_state["src"] == "From this Omnix session"
    assert len([c for c in at.sidebar.checkbox if c.key and c.key.startswith("ib_use_")]) == 3
    next(b for b in at.sidebar.button if b.label == "Run analysis").click().run(timeout=900)
    assert not at.exception and "trace" in at.session_state
    _open(at, "proteomics")                                       # the inbox survives switching apps
    assert integrix_bridge.INBOX in at.session_state


def test_proteomics_send_to_integrix_after_normalization():
    at = _portal({"access_keys": {TEST_HASH: {"name": "Lab", "plan": "Multi-Omics", "expires": "2099-12-31"}}},
                 timeout=600)
    _activate(at, TEST_KEY)
    _open(at, "proteomics")
    assert "Send to Integrix" not in [b.label for b in at.sidebar.button]            # before Normalization

    def click(label):
        next(b for b in at.button if b.label.startswith(label) or label in b.label).click().run()

    def step(group, page):
        at.session_state["nav_group"] = group
        at.session_state["nav_page"] = page
        at.run()
    click("Load & Validate Data")
    step("Preprocessing", "Data Cleaning & Imputation")
    click("Apply Missingness Filter")
    click("Apply Imputation")
    step("Preprocessing", "QC Validation")
    click("Confirm QC & Proceed")
    step("Normalization", "Normalization")
    click("Apply Log2 Transformation")
    click("Apply Median Centering Normalization")
    assert not at.exception
    next(b for b in at.sidebar.button if b.label == "Send to Integrix").click().run()
    inbox = at.session_state[integrix_bridge.INBOX]
    assert "proteomics" in inbox and inbox["proteomics"]["n_samples"] >= 3
    with zipfile.ZipFile(io.BytesIO(inbox["proteomics"]["zip"])) as z:
        assert {"proteomics.csv", "metadata.csv", "proteomics_omnix.json"} <= set(z.namelist())
    assert "Open Integrix  →" in [b.label for b in at.sidebar.button]
