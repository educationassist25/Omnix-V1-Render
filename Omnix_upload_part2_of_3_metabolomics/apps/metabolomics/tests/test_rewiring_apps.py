"""
Headless UI tests (Streamlit AppTest) for the standalone Rewiring Map app and the
Network Analysis page inside the main MetaboAI Pro app.

    python -m pytest tests/test_rewiring_apps.py -q
"""

import os
import sys

import numpy as np
import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from metabolomics_modules import rewiring as rw, utils  # noqa: E402

STANDALONE = os.path.join(ROOT, "rewiring_app.py")
MAIN = os.path.join(ROOT, "app.py")


def _run_button(at):
    return next(b for b in at.button if b.label.startswith("Run Rewiring"))


def _texts(elements):
    return " ".join(str(e.value) for e in elements)


# ---------------------------------------------------------------------------
# standalone app
# ---------------------------------------------------------------------------
def test_standalone_demo_runs_end_to_end():
    at = AppTest.from_file(STANDALONE, default_timeout=180).run()
    assert not at.exception
    assert "Run Rewiring Analysis" in [b.label for b in at.button]
    assert "click **Run Rewiring Analysis**" in _texts(at.info)
    _run_button(at).click().run()
    assert not at.exception
    msg = _texts(at.success)
    assert "rewired pairs" in msg and "Global test p" in msg and "permutation FDR" in msg
    assert "Pathway level:" in msg and "KEGG pathways rewired" in msg
    assert at.get("iframe"), "interactive view should be embedded"
    labels = [d.label for d in at.get("download_button")]
    for want in ("Rewired pairs (CSV)", "All pairs (CSV)", "Hypothesis evidence (JSON)", "Interactive report (HTML)",
                 "Pathway-level tests (CSV)", "Reaction pairs (CSV)", "Metabolite ID / pathway mapping (CSV)"):
        assert want in labels
    # figures download from the buttons right under each image in the view, not from a separate section
    assert not [x for x in labels if x.startswith(("Download network", "Download correlation plot"))]


def test_standalone_settings_change_prompts_rerun():
    at = AppTest.from_file(STANDALONE, default_timeout=180).run()
    _run_button(at).click().run()
    at.radio(key="rwapp_0_method").set_value("Spearman").run()
    assert not at.exception
    assert "Settings changed since the last run" in _texts(at.caption)
    _run_button(at).click().run()
    assert "Settings changed" not in _texts(at.caption)


def test_standalone_rich_demo_warns_about_small_groups():
    at = AppTest.from_file(STANDALONE, default_timeout=180).run()
    at.sidebar.radio[0].set_value("Omnix Metabolomics demo: Rich Clinical (4 groups × 9 samples)").run()
    assert not at.exception
    assert "9 samples in the smaller group" in _texts(at.caption)
    assert not at.warning
    _run_button(at).click().run()
    assert not at.exception
    assert at.get("iframe")


def test_standalone_upload_mode_without_files_shows_instructions():
    at = AppTest.from_file(STANDALONE, default_timeout=60).run()
    at.sidebar.radio[0].set_value("Upload my data").run()
    assert not at.exception
    assert "Upload a peak area matrix" in _texts(at.info)


# ---------------------------------------------------------------------------
# main app integration
# ---------------------------------------------------------------------------
def _goto_rewiring(at):
    at.session_state["nav_group"] = "Network Analysis"
    at.session_state["nav_page"] = "Correlation Rewiring Map"
    return at


def test_main_app_home_page_still_loads():
    at = AppTest.from_file(MAIN, default_timeout=180).run()
    assert not at.exception
    assert "Network Analysis" in [b.label for b in at.sidebar.button]


def test_main_app_page_without_data_offers_demo():
    at = _goto_rewiring(AppTest.from_file(MAIN, default_timeout=180)).run()
    assert not at.exception
    assert "Complete **Normalization** first" in _texts(at.info)
    at.toggle(key="rewire_use_demo").set_value(True).run()
    assert not at.exception
    _run_button(at).click().run()
    assert not at.exception
    assert "rewired pairs" in _texts(at.success)


def _pipeline_like_state(prefix_names=False):
    """What Normalization leaves in session_state: log2 bio-sample matrix, validated
    metadata indexed by Sample, row annotations indexed by Metabolite."""
    d = rw.simulate_rewiring_demo()
    peak = d["peak_df"]
    meta = utils.validate_metadata(d["meta"].copy(), peak.columns)
    log2 = np.log2(peak)
    if prefix_names:
        log2.index = ["Targeted Metabolomics::" + n for n in log2.index]
    ann = d["annotations"].set_index("Metabolite")
    return log2, meta, ann


@pytest.mark.parametrize("multi", [False, True])
def test_main_app_page_uses_normalized_pipeline_data(multi):
    log2, meta, ann = _pipeline_like_state(prefix_names=multi)
    at = _goto_rewiring(AppTest.from_file(MAIN, default_timeout=180))
    at.session_state["log2_data"] = log2
    at.session_state["meta"] = meta
    at.session_state["row_annotations"] = ann
    at.session_state["data_mode"] = "multi" if multi else "single"
    at.run()
    assert not at.exception
    assert ("combined normalized data" if multi else "normalized log2 data") in _texts(at.caption)
    assert at.selectbox(key="rewire_main_pw_col").value == "KEGG (as in MSEA)"
    assert "rewire_main_batch_col" not in [s.key for s in at.selectbox]      # no batch-column option
    assert "Batch column (for batch check)" not in [s.label for s in at.selectbox]
    _run_button(at).click().run()
    assert not at.exception
    msg = _texts(at.success)
    assert "rewired pairs" in msg and "43 metabolites" in msg
    stored = at.session_state["rewire_main_result"]["result"]
    assert stored.prep.pathway_source == "library"             # MSEA KEGG library
    assert stored.prep.logged is False                          # already log2 -> not transformed again
    keys = {h["key"] for h in stored.hypotheses}
    assert "Alanine, aspartate and glutamate metabolism|Citrate cycle (TCA cycle)|lost" in keys
    assert not any("::" in f for f in stored.prep.features)
    # only real pathway libraries are offered (no annotation-column or clustering options)
    assert at.selectbox(key="rewire_main_pw_col").options == [
        "KEGG (as in MSEA)", "SMPDB (as in MSEA)", "Human-GEM subsystems", "LIPID MAPS (as in MSEA)"]


def test_main_app_other_state_untouched_by_rewiring_page():
    log2, meta, ann = _pipeline_like_state()
    at = _goto_rewiring(AppTest.from_file(MAIN, default_timeout=180))
    at.session_state["log2_data"] = log2
    at.session_state["meta"] = meta
    at.session_state["row_annotations"] = ann
    at.run()
    _run_button(at).click().run()
    assert at.session_state["stats_result"] is None
    pd.testing.assert_frame_equal(at.session_state["log2_data"], log2)


def test_standalone_reaction_scope_and_three_per_group():
    """Reaction-linked scope runs from the UI, and a 3 + 3 dataset runs with the exact permutation note."""
    at = AppTest.from_file(STANDALONE, default_timeout=180).run()
    at.radio(key="rwapp_0_scope").set_value("Reaction-linked pairs (≤ 2 steps)").run()
    _run_button(at).click().run()
    assert not at.exception
    res = at.session_state["rwapp_0_result"]["result"]
    assert res.cfg.pair_scope == "reaction"
    assert res.pairs["tested"].sum() == (res.pairs["rxn_steps"] > 0).sum() > 0
    assert "(reaction-linked)" in _texts(at.success)


def test_standalone_network_controls_figure_formats_and_fdr_cutoff():
    at = AppTest.from_file(STANDALONE, default_timeout=180).run()
    at.select_slider(key="rwapp_0_fdr").set_value(1.0).run()
    assert "A high FDR q-value cutoff" in _texts(at.caption)
    _run_button(at).click().run()
    assert not at.exception
    res = at.session_state["rwapp_0_result"]["result"]
    assert res.cfg.cutoff_on == "fdr" and res.cfg.fdr == 1.0
    assert "FDR q-value" in _texts(at.success) and "≤ 1." in _texts(at.success)
    # every tested pair with |dr| >= min_delta is called at cutoff 1
    t = res.pairs[res.pairs["tested"]]
    assert (t["delta_r"].abs() >= res.cfg.min_delta).sum() == t["class"].isin(("lost", "gained", "flipped")).sum()
    # network default: top 10 pathways by p-value <= 0.05
    assert at.selectbox(key="rwapp_0_net_top").value == "Top 10"
    assert at.radio(key="rwapp_0_net_metric").value == "p-value"
    assert at.select_slider(key="rwapp_0_net_cut").value == 0.05
    at.selectbox(key="rwapp_0_net_top").set_value("Top 5").run()
    at.radio(key="rwapp_0_net_metric").set_value("FDR q-value").run()
    assert not at.exception
    # figures download from the buttons under each image in the view; no separate figure-format section
    assert not [e for e in at.expander if "More figure formats" in e.label]
    assert "rwapp_0_fig_fmt" not in [s.key for s in at.selectbox]


def test_network_page_writes_p_value_and_fdr_q_value_in_full():
    """No bare 'q' anywhere on the page: options, summary, tables, downloads and the embedded view;
    the pair cutoff is always on the FDR q-value (no 'Rewired-pair cutoff on' choice)."""
    at = AppTest.from_file(STANDALONE, default_timeout=180).run()
    assert "rwapp_0_cut_on" not in [r.key for r in at.radio]
    assert "Rewired-pair cutoff on" not in [r.label for r in at.radio]
    assert at.select_slider(key="rwapp_0_fdr").label.startswith("Rewired-pair FDR q-value cutoff")
    _run_button(at).click().run()
    assert "FDR q-value" in _texts(at.success) and "Global test p-value" in _texts(at.success)
    assert at.radio(key="rwapp_0_net_metric").options == ["p-value", "FDR q-value"]
    res = at.session_state["rwapp_0_result"]["result"]
    from metabolomics_modules import rewiring_ui
    cols = list(rewiring_ui.labelled(res.rewired_pairs()).columns)
    assert "FDR q-value (used)" in cols and "q_value" not in cols and "q_bh" not in cols
    pcols = list(rewiring_ui.labelled(res.pathway_table_display(), "pathways", ("Control", "Tumor")).columns)
    assert "FDR q-value (BH)" in pcols and "p-value (permutation)" in pcols
    assert not [c for c in pcols if "_" in c or c.lower() in ("q", "p")], pcols      # every column written out
    assert "Mean |r| Control" in pcols
    html = rewiring_ui.render_view_html(res)
    assert '"pw_metric": "p"' in html and '"pw_cutoff": 0.05' in html and '"pw_top": null' in html
    assert "FDR q-value" in html and ">FDR q<" not in html and '<th class="n">q</th>' not in html
