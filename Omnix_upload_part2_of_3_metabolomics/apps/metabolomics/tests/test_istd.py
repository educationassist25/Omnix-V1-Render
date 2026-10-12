"""
ISTD normalization for Targeted Metabolomics and Targeted Lipidomics: one internal standard per
dataset, linked to every feature in the row annotations' ISTD column.

    python -m pytest tests/test_istd.py -q
"""

import os
import sys

import numpy as np
import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from metabolomics_modules import dataset_manager, normalization as nz, utils  # noqa: E402

TARGETED_FILES = [
    # matrix, annotation file, Method filter, the dataset's ISTD, number of analytes
    ("sample_peak_area_matrix_targeted.csv", "metabolite_row_annotations_targeted.csv", None, "ISTD_D4-Alanine", 59),
    ("sample_peak_area_matrix_multimethod_targeted_metabolomics.csv", "metabolite_row_annotations_multimethod.csv",
     "Targeted Metabolomics", "ISTD_D4-Alanine", 49),
    ("sample_peak_area_matrix_multimethod_targeted_lipidomics.csv", "metabolite_row_annotations_multimethod.csv",
     "Targeted Lipidomics", "ISTD_PC(15:0/18:1-d7)", 39),
    ("sample_peak_area_matrix_richdemo.csv", "metabolite_row_annotations_richdemo.csv", None, "ISTD_D4-Alanine", 200),
]


def _read(name):
    return pd.read_csv(os.path.join(ROOT, name))


# ---------------------------------------------------------------------------
# Demo data: one ISTD per targeted dataset, every feature linked to it
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("matrix, annot, method, istd, n_analytes", TARGETED_FILES)
def test_demo_has_one_istd_and_every_feature_is_linked_to_it(matrix, annot, method, istd, n_analytes):
    peak = utils.validate_peak_matrix(_read(matrix))
    ann = _read(annot)
    if method:
        ann = ann[ann["Method"] == method]
    ann = ann.set_index("Metabolite")
    assert set(ann.index) == set(peak.index)
    assert nz.detect_istds(peak.index, ann) == [istd] and len(peak) - 1 == n_analytes
    assert (ann["ISTD"] == istd).all()                         # every feature (and the ISTD row) -> the ISTD
    assert nz.default_istd(peak.index, ann, method, method) == istd
    vals = peak.loc[istd]
    assert (vals > 0).all() and vals.std() / vals.mean() * 100 < 10


def test_demo_lipids_use_a_lipid_standard():
    peak = _read("sample_peak_area_matrix_multimethod_targeted_lipidomics.csv")
    assert "ISTD_D4-Alanine" not in set(peak["Metabolite"])
    assert peak["Metabolite"].iloc[0] == "ISTD_PC(15:0/18:1-d7)"


def test_demo_untargeted_files_have_no_istds():
    for m in ["sample_peak_area_matrix_untargeted.csv", "sample_peak_area_matrix_multimethod_untargeted_metabolomics.csv",
              "sample_peak_area_matrix_multimethod_untargeted_lipidomics.csv"]:
        assert not nz.detect_istds(_read(m)["Metabolite"])


# ---------------------------------------------------------------------------
# Mapping and normalization math
# ---------------------------------------------------------------------------
def _toy():
    peak = pd.DataFrame(
        {"S1": [100., 200., 300., 400., 10., 20.], "S2": [120., 180., 330., 500., 20., 40.]},
        index=["A", "B", "C", "D", "ISTD_x", "ISTD_y"])
    ann = pd.DataFrame({"Class": ["c1", "c2", "c3", "c1", "c1", "c2"],
                        "ISTD": ["ISTD_y", "", "", "", "ISTD_x", "ISTD_y"]}, index=peak.index)
    return peak, ann


def test_mapping_priority_column_then_class_then_pooled():
    peak, ann = _toy()
    m = nz.build_istd_mapping(peak.index, ["ISTD_x", "ISTD_y"], ann)
    assert m.at["A", "ISTD"] == "ISTD_y" and m.at["A", "Matched by"] == nz.MATCH_COLUMN   # column beats class
    assert m.at["B", "ISTD"] == "ISTD_y" and m.at["B", "Matched by"] == nz.MATCH_CLASS
    assert m.at["D", "ISTD"] == "ISTD_x" and m.at["D", "Matched by"] == nz.MATCH_CLASS
    assert m.at["C", "ISTD"] == nz.POOLED_ISTD
    assert "ISTD_x" not in m.index
    single = nz.build_istd_mapping(peak.index, ["ISTD_x"], ann)
    assert single.at["A", "ISTD"] == "ISTD_x" and single.at["A", "Matched by"] == nz.MATCH_CLASS
    assert single.at["A", "Annotated ISTD"] == "ISTD_y"                # asked for, but not selected
    assert single.at["B", "ISTD"] == "ISTD_x" and single.at["B", "Matched by"] == nz.MATCH_SINGLE


def test_normalized_values_are_ratios_to_the_matched_istd():
    peak, ann = _toy()
    out, m = nz.istd_normalize_mapped(peak, ["ISTD_x", "ISTD_y"], ann)
    assert list(out.index) == ["A", "B", "C", "D"]
    np.testing.assert_allclose(out.loc["A"], peak.loc["A"] / peak.loc["ISTD_y"])
    np.testing.assert_allclose(out.loc["D"], peak.loc["D"] / peak.loc["ISTD_x"])
    pooled = np.sqrt(peak.loc["ISTD_x"] * peak.loc["ISTD_y"])
    np.testing.assert_allclose(out.loc["C"], peak.loc["C"] / pooled)
    # a zero ISTD value gives NaN only for its own analytes, in that sample
    peak.loc["ISTD_x", "S1"] = 0
    out, _ = nz.istd_normalize_mapped(peak, ["ISTD_x", "ISTD_y"], ann)
    assert np.isnan(out.at["D", "S1"]) and np.isnan(out.at["C", "S1"]) and not np.isnan(out.at["A", "S1"])
    # the legacy single-ISTD function and the wrapper agree with the mapped version
    old = nz.istd_normalize(peak, "ISTD_y")
    new, _ = nz.istd_normalize_mapped(peak, ["ISTD_y"])
    pd.testing.assert_frame_equal(old, new, check_names=False)
    _, unmatched = nz.istd_normalize_multi(peak, ["ISTD_x", "ISTD_y"], ann)
    assert unmatched == ["C"]


def test_unselected_istds_are_dropped_and_lookups_handle_prefixes_and_duplicates():
    peak, ann = _toy()
    out, _ = nz.istd_normalize_mapped(peak, ["ISTD_x"], ann, drop_features=["ISTD_y"])
    assert "ISTD_y" not in out.index
    # multi-dataset annotation files may use "<label>::<feature>" names, and one plain name may occur
    # for several methods: the row of the dataset's own Method wins
    ann2 = pd.DataFrame({"ISTD": ["ISTD_x", "ISTD_y", "ISTD_y"],
                         "Method": ["Targeted Metabolomics", "Targeted Lipidomics", "Targeted Lipidomics"]},
                        index=["A", "A", "LipidRun::B"])
    m = nz.build_istd_mapping(["A", "B"], ["ISTD_x", "ISTD_y"], ann2, label="LipidRun", method="Targeted Lipidomics")
    assert m.at["A", "ISTD"] == "ISTD_y" and m.at["B", "ISTD"] == "ISTD_y"
    m = nz.build_istd_mapping(["A"], ["ISTD_x", "ISTD_y"], ann2, method="Targeted Metabolomics")
    assert m.at["A", "ISTD"] == "ISTD_x"
    # detection: by name, and by being named in an 'Internal Standard' column
    ann3 = pd.DataFrame({"Internal Standard": ["PC 15:0-18:1(d7)", ""]}, index=["PC(16:0/18:1)", "Glc"])
    assert nz.detect_istds(["PC 15:0-18:1(d7)", "PC(16:0/18:1)", "IS_Leu-d3", "Glc", "Isoleucine"], ann3) == \
        ["PC 15:0-18:1(d7)", "IS_Leu-d3"]


def test_keep_istds_restores_filtered_standards():
    peak, _ = _toy()
    filtered = peak.loc[["A", "C"]]
    kept = nz.keep_istds(filtered, peak, ["ISTD_x", "ISTD_y"])
    assert list(kept.index) == ["A", "C", "ISTD_x", "ISTD_y"]


# ---------------------------------------------------------------------------
# Auto-Process (multi-dataset Quick Start)
# ---------------------------------------------------------------------------
def _multi_datasets():
    meta = utils.validate_metadata(_read("sample_metadata_multimethod.csv"),
                                   utils.validate_peak_matrix(_read("sample_peak_area_matrix_targeted.csv")).columns)
    ann = _read("metabolite_row_annotations_multimethod.csv").set_index("Metabolite")
    qc_cols, sample_cols = utils.split_qc_and_samples(None, meta)
    datasets = {}
    for dtype in ["Untargeted Metabolomics", "Targeted Metabolomics", "Untargeted Lipidomics", "Targeted Lipidomics"]:
        e = dataset_manager.make_dataset_entry(dtype, dtype)
        e["raw_df"] = utils.validate_peak_matrix(
            _read(f"sample_peak_area_matrix_multimethod_{dtype.lower().replace(' ', '_')}.csv"))
        e["qc_cols"], e["sample_cols"] = qc_cols, sample_cols
        datasets[dtype] = e
    return datasets, ann


def test_auto_process_normalizes_each_targeted_dataset_with_its_one_istd():
    datasets, ann = _multi_datasets()
    # give the lipid ISTD a QC CV above 20%: the CV filter must still keep it
    lip = datasets["Targeted Lipidomics"]
    lip["raw_df"].loc["ISTD_PC(15:0/18:1-d7)", lip["qc_cols"]] = [1e4, 5e4, 1e4, 5e4, 1e4, 5e4]
    notes = dataset_manager.auto_process_all_datasets(datasets, row_annotations=ann)
    for dtype, istd in [("Targeted Metabolomics", "ISTD_D4-Alanine"), ("Targeted Lipidomics", "ISTD_PC(15:0/18:1-d7)")]:
        ds = datasets[dtype]
        assert istd in ds["raw_peak_df_qc"].index
        assert set(ds["istd_mapping"]["ISTD"]) == {istd}
        assert not any("ISTD" in f for f in ds["log2_data"].index)
        base = ds["raw_peak_df_qc"]
        feats = ds["log2_data"].index
        expected = np.log2(base.loc[feats].astype(float).div(base.loc[istd].astype(float), axis=1))
        np.testing.assert_allclose(ds["log2_data"].values.astype(float), expected.values, equal_nan=True)
        assert any(f"using '{istd}'" in n for n in notes)


def test_auto_process_untargeted_drops_istd_rows_and_targeted_without_istd_fails_clearly():
    peak = utils.validate_peak_matrix(_read("sample_peak_area_matrix_richdemo.csv"))
    meta = utils.validate_metadata(_read("sample_metadata_richdemo.csv"), peak.columns)
    ann = _read("metabolite_row_annotations_richdemo.csv").set_index("Metabolite")
    qc_cols, sample_cols = utils.split_qc_and_samples(peak, meta)
    e = dataset_manager.make_dataset_entry("Rich", "Untargeted Metabolomics")
    e.update(raw_df=peak, qc_cols=qc_cols, sample_cols=sample_cols)
    dataset_manager.auto_process_dataset(e, row_annotations=ann)
    assert not any("ISTD" in f for f in e["log2_data"].index) and len(e["log2_data"]) <= 200
    e2 = dataset_manager.make_dataset_entry("Rich", "Targeted Metabolomics")
    e2.update(raw_df=peak.drop(index=nz.detect_istds(peak.index, ann)), qc_cols=qc_cols, sample_cols=sample_cols)
    with pytest.raises(ValueError, match="no internal standard"):
        dataset_manager.auto_process_dataset(e2, row_annotations=None)


# ---------------------------------------------------------------------------
# The app: Normalization page, single dataset
# ---------------------------------------------------------------------------
def _goto(at, group, page):
    at.session_state["nav_group"] = group
    at.session_state["nav_page"] = page
    return at.run()


def _button(at, text):
    return next(b for b in at.button if text in b.label)


@pytest.mark.parametrize("dtype, demo, istd, n_analytes", [
    ("Targeted Metabolomics", "Standard", "ISTD_D4-Alanine", 59),
    ("Targeted Lipidomics", "Standard", "ISTD_PC(15:0/18:1-d7)", 39),
    ("Targeted Metabolomics", "Rich Clinical Demo", "ISTD_D4-Alanine", 200),
])
def test_app_targeted_istd_step_preselects_the_istd(dtype, demo, istd, n_analytes):
    at = AppTest.from_file(os.path.join(ROOT, "app.py"), default_timeout=300).run()
    at.selectbox[0].set_value(dtype).run()
    at.radio(key="single_demo_choice").set_value(demo).run()
    _button(at, "Load & Validate Data").click().run()
    assert not at.exception
    _goto(at, "Normalization", "Normalization")
    assert not at.exception
    assert at.selectbox(key="single_istd_choice").value == istd
    _button(at, "Apply ISTD Normalization").click().run()
    assert not at.exception
    norm = at.session_state["istd_normalized"]
    assert norm.shape[0] == n_analytes and istd not in norm.index
    assert set(at.session_state["istd_mapping"]["ISTD"]) == {istd}
    raw = at.session_state["raw_peak_df"]
    expected = raw.loc[norm.index, norm.columns].astype(float).div(raw.loc[istd, norm.columns].astype(float), axis=1)
    np.testing.assert_allclose(norm.values.astype(float), expected.values)
    _button(at, "Apply Strict Log2 Transformation").click().run()
    assert not at.exception and at.session_state["log2_data"].shape[0] == n_analytes


def test_app_untargeted_leaves_istd_rows_out():
    at = AppTest.from_file(os.path.join(ROOT, "app.py"), default_timeout=300).run()
    at.selectbox[0].set_value("Untargeted Metabolomics").run()
    at.radio(key="single_demo_choice").set_value("Rich Clinical Demo").run()
    _button(at, "Load & Validate Data").click().run()
    _goto(at, "Normalization", "Normalization")
    _button(at, "Apply Strict Log2 Transformation").click().run()
    _button(at, "Apply Median-IQR Normalization").click().run()
    assert not at.exception
    out = at.session_state["log2_data"]
    assert out.shape[0] == 200 and not any("ISTD" in f for f in out.index)
