"""Automated checks of the core engine. Run:  python -m pytest -q tests"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from core import correction as cr, stats as stt, kinetics as kn, labeling as lb  # noqa: E402
from core import dataset as dsm, normalization as nm, pathmap as pmap, pathway_templates as ptpl  # noqa: E402
from core import design as dz, pathways as pw, interpret as ip  # noqa: E402
from core.io_reorganize import parse_filenames, process_files, TYPE_LAB, TYPE_REF, TYPE_BLANK  # noqa: E402



# ---------------- statistics ----------------
def test_bh_matches_statsmodels():
    from statsmodels.stats.multitest import multipletests
    p = np.random.default_rng(0).uniform(size=300)
    assert np.allclose(stt.bh_fdr(p), multipletests(p, method="fdr_bh")[1])


def test_bh_keeps_nan():
    out = stt.bh_fdr([0.01, np.nan, 0.04])
    assert np.isnan(out[1]) and np.isfinite(out[[0, 2]]).all()


def test_welch_matches_scipy_and_direction():
    from scipy import stats
    rng = np.random.default_rng(1)
    X = pd.DataFrame(rng.normal(10, 1, (20, 6)), columns=list("abcdef"))
    X.loc[:, ["d", "e", "f"]] += 2  # B higher
    T = stt.welch_table(X, ["a", "b", "c"], ["d", "e", "f"], "A", "B")
    for k in range(5):
        t, p = stats.ttest_ind(X.loc[k, ["d", "e", "f"]], X.loc[k, ["a", "b", "c"]], equal_var=False)
        assert np.isclose(T.p[k], p) and np.isclose(T.t[k], t)
    assert (T.log2FC > 0).mean() > 0.9
    assert np.allclose(T.FC, 2 ** T.log2FC)


def test_welch_warnings():
    X = pd.DataFrame([[1, 1, 1, 2, 3, 4], [1, np.nan, np.nan, 2, 3, 4]], columns=list("abcdef"))
    T = stt.welch_table(X, ["a", "b", "c"], ["d", "e", "f"])
    assert "zero variance" in T.Warnings[0]
    assert "not tested" in T.Warnings[1] and np.isnan(T.p[1])


def test_tukey_matches_statsmodels():
    from statsmodels.stats.multicomp import pairwise_tukeyhsd
    rng = np.random.default_rng(2)
    g = {"x": rng.normal(0, 1, 4), "y": rng.normal(1, 1, 5), "z": rng.normal(2, 1, 3)}
    mine = stt.tukey_kramer(g)
    vals = np.concatenate(list(g.values())); lab = np.repeat(list(g), [len(v) for v in g.values()])
    ref = pairwise_tukeyhsd(vals, lab).summary().data[1:]
    assert np.allclose(mine["p (Tukey)"], [r[3] for r in ref], atol=1e-3)
    assert np.allclose(mine["Difference (B−A)"], [r[2] for r in ref], atol=1e-4)


def test_anova_matches_scipy():
    from scipy import stats
    rng = np.random.default_rng(3)
    X = pd.DataFrame(rng.normal(0, 1, (5, 9)), columns=[f"s{i}" for i in range(9)])
    groups = {"A": ["s0", "s1", "s2"], "B": ["s3", "s4", "s5"], "C": ["s6", "s7", "s8"]}
    T, TK = stt.anova_table(X, groups)
    for k in range(5):
        F, p = stats.f_oneway(*[X.loc[k, s] for s in groups.values()])
        assert np.isclose(T["p (ANOVA)"][k], p)
    assert len(TK) == 5 * 3


# ---------------- natural abundance correction ----------------
def test_correction_round_trip():
    true = np.array([0.4, 0.05, 0.3, 0.1, 0.05, 0.06, 0.04])
    A = cr.correction_matrix("C6H8O7", 7, purity=0.99)
    assert np.allclose(cr.correct_mid(A @ true, "C6H8O7", purity=0.99), true, atol=1e-9)


def test_correction_unlabeled_becomes_m0():
    nat = cr.natural_mid("C6H12O6", 7)
    c = cr.correct_mid(nat, "C6H12O6", purity=1.0)
    assert c[0] > 0.999


def test_natural_m1_citrate():
    assert abs(cr.natural_mid("C6H8O7", 3)[1] - 0.0608) < 0.002


# ---------------- labeling metrics ----------------
def _fake_rows():
    cols = {f"Relative Area {i}": 0.0 for i in range(26)}
    r = dict(Name="X", Mode="RP_Pos", Formula="C3 H6 O3", nC=3, Sample="s1", **cols)
    r["Relative Area 0"], r["Relative Area 1"], r["Relative Area 2"], r["Relative Area 3"] = 20, 10, 40, 30
    return pd.DataFrame([r])


def test_fraction_example_from_spec():
    F, S = lb.compute_labeling(_fake_rows(), ["s1"])
    assert np.allclose(F.Fraction.values, [0.2, 0.1, 0.4, 0.3])
    assert np.isclose(S.Incorporation[0], 80.0)
    assert np.isclose(S.MeanLabels[0], 0.1 + 0.8 + 0.9)
    assert np.isclose(S.Enrichment[0], 1.8 / 3 * 100)
    assert np.isclose(S.SumCheck[0], 1.0)


# ---------------- v2: imputation and normalization ----------------
def _matrix(seed=4, miss=0.1):
    rng = np.random.default_rng(seed)
    base = rng.normal(20, 2, (60, 1))
    L = base + rng.normal(0, 0.3, (60, 9))
    X = pd.DataFrame(2 ** L, index=[f"m{i} | M+{i % 3}" for i in range(60)], columns=[f"s{j}" for j in range(9)])
    M = rng.uniform(size=X.shape) < miss
    return X.mask(M), M


def test_robust_iqr_equalizes_samples_and_keeps_log2_units():
    X, _ = _matrix(miss=0.0)
    X["s1"] *= 4                                                   # +2 log2 offset
    L, bad = nm.log2_transform(X)
    N, info = nm.robust_iqr(L)
    assert bad == 0
    assert N.median().std() < 1e-9 and (N.quantile(.75) - N.quantile(.25)).std() < 1e-9
    assert abs(float(np.median(info["IQR (log2)"])) - float((N.quantile(.75) - N.quantile(.25)).iloc[0])) < 1e-9


# ---------------- v2: statistics ----------------
def test_p_adjust_matches_statsmodels():
    from statsmodels.stats.multitest import multipletests
    p = np.random.default_rng(0).uniform(size=200) ** 3
    for m, sm in (("bh", "fdr_bh"), ("by", "fdr_by"), ("holm", "holm"), ("bonferroni", "bonferroni")):
        assert np.allclose(stt.p_adjust(p, m), multipletests(p, method=sm)[1])


def test_two_way_anova_matches_statsmodels():
    import statsmodels.api as sm
    import statsmodels.formula.api as smf
    rng = np.random.default_rng(1)
    D = pd.DataFrame({"Sample": [f"s{i}" for i in range(18)], "Group": np.repeat(["a", "b"], 9),
                      "TimeLabel": np.tile(np.repeat(["12", "18", "40"], 3), 2)})
    X = pd.DataFrame(rng.normal(0, 1, (1, 18)), columns=D.Sample)
    X.iloc[0, :9] += 1
    X.iloc[0, 3] = np.nan
    T = stt.two_way_anova(X, D)
    d = D.assign(y=X.iloc[0].values).dropna()
    ref = sm.stats.anova_lm(smf.ols("y ~ C(Group) * C(TimeLabel)", d).fit(), typ=2)
    assert np.isclose(T["p Group"][0], ref["PR(>F)"].iloc[0])
    assert np.isclose(T["p TimeLabel"][0], ref["PR(>F)"].iloc[1])
    assert np.isclose(T["p Group×TimeLabel"][0], ref["PR(>F)"].iloc[2])


def test_games_howell_and_welch_anova_reference_values():
    # reference values from pingouin 0.5 (pairwise_gameshowell, welch_anova) for these data
    g = {"x": np.array([0.1, -0.5, 0.9, 0.3]), "y": np.array([1.8, 0.2, 2.9, 0.4, 1.1]),
         "z": np.array([2.2, 1.7, 2.6])}
    GH = stt.games_howell(g)
    assert np.allclose(GH["Difference (B−A)"], [g["y"].mean() - g["x"].mean(), g["z"].mean() - g["x"].mean(),
                                                g["z"].mean() - g["y"].mean()])
    try:
        import pingouin as pg
    except Exception:
        pytest.skip("pingouin not installed")
    df = pd.DataFrame([(k, v) for k, vs in g.items() for v in vs], columns=["g", "v"])
    ref = pg.pairwise_gameshowell(df, dv="v", between="g")
    assert np.allclose(GH["p (Games-Howell)"], ref.pval, atol=1e-6)
    assert np.isclose(stt._welch_anova(list(g.values()))[1], pg.welch_anova(df, dv="v", between="g")["p_unc"][0])


def test_contrast_builder():
    D = dz.finalize_design(pd.DataFrame({
        "Sample": [f"s{i}" for i in range(12)], "SampleType": TYPE_LAB, "Group": np.repeat(["A", "B"], 6),
        "Treatment": "", "TimeValue": np.tile(np.repeat([1.0, 2.0], 3), 2), "TimeUnit": "h", "Include": True,
        "Replicate": 1}))
    c = dz.build_contrasts(D, "arms_per_time", "A", "B")
    assert [x["name"] for x in c] == ["B vs A · 1 h", "B vs A · 2 h"] and c[0]["samples_b"] == ["s6", "s7", "s8"]
    assert len(dz.build_contrasts(D, "time_vs_first")) == 2
    D2 = D.copy(); D2.loc[D2.Sample == "s0", "Treatment"] = "Drug"
    assert "A + Drug" in dz.groups_of(dz.finalize_design(D2))


# ---------------- v2: data model, no summing ----------------
def test_generic_table_keeps_every_isotopologue_and_never_sums():
    import io as _io
    T = pd.DataFrame({"Metabolite": ["Lac"] * 3, "Formula": ["C3H6O3"] * 3, "Isotopologue": ["M+0", "M+1", "M+3"],
                      "s1": [100.0, 5.0, 40.0], "s2": [90.0, 4.0, 60.0]})
    d = dsm.read_generic(_io.BytesIO(T.to_csv(index=False).encode()), "t.csv", min_area=0)
    assert list(d.CT.i) == [0, 1, 3]                               # one row per isotopologue, never summed
    assert d.CT.loc[d.CT.i == 3, "s1"].iloc[0] == 40.0
    dup = pd.concat([T, T.iloc[[0]]])
    with pytest.raises(ValueError):
        dsm.read_generic(_io.BytesIO(dup.to_csv(index=False).encode()), "t.csv")
    F, S = lb.compute_from_ct(d.CT, ["s1", "s2"])
    assert np.allclose(F[F.Sample == "s1"].Fraction, [100 / 145, 5 / 145, 40 / 145])


def test_vectorized_labeling_matches_reference_implementation():
    F, S = lb.compute_labeling(_fake_rows(), ["s1"])
    CT = pd.DataFrame({"Compound": "X", "Formula": "C3 H6 O3", "Method": "RP_Pos", "n": 3, "i": [0, 1, 2, 3],
                       "Channel": [f"X | M+{i}" for i in range(4)], "s1": [20.0, 10, 40, 30]})
    F2, S2 = lb.compute_from_ct(CT, ["s1"])
    assert np.allclose(F2.Fraction, F.Fraction) and np.isclose(S2.Enrichment[0], S.Enrichment[0])


# ---------------- v2.5: curated pathway maps ----------------
def _blocks_for(binding, cols, n=5):
    rng = np.random.default_rng(0)
    mk = lambda: pd.DataFrame(rng.random((n, len(cols))) * 50, index=range(n),
                              columns=pd.MultiIndex.from_tuples(cols, names=["Group", "Hours"]))
    comps = set(binding.values())
    return {"iso": {c: mk() for c in comps}, "z": {c: mk() - 25 for c in comps}}, {c: list(range(n)) for c in comps}


@pytest.mark.parametrize("key", list(ptpl.T))
def test_pathway_template_layout_has_no_overlaps(key):
    tpl = ptpl.T[key]
    for e in tpl.edges:
        assert e["a"] in tpl.nodes and e["b"] in tpl.nodes
    binding = {nid: f"cmp_{nid}" for k, (nid, n) in enumerate(tpl.nodes.items()) if n["kind"] == "met" and k % 2 == 0}
    cols = [("A", 1.0), ("A", 4.0), ("B", 1.0), ("B", 4.0)]
    blocks, rows = _blocks_for(binding, cols)
    pm_ = pmap.layout(tpl, binding, rows, cols, ["iso", "z"])
    R = pm_.nodes
    x0, x1 = (R.cx - R.w / 2).to_numpy(), (R.cx + R.w / 2).to_numpy()
    y0, y1 = (R.cy - R.h / 2).to_numpy(), (R.cy + R.h / 2).to_numpy()
    for i in range(len(R)):
        for j in range(i + 1, len(R)):
            assert not (x0[i] < x1[j] and x0[j] < x1[i] and y0[i] < y1[j] and y0[j] < y1[i]), (R.id[i], R.id[j])
    px0, py0, px1, py1 = pm_.panel                 # everything inside the panel
    assert x0.min() >= px0 and x1.max() <= px1 and y0.min() >= py0 and y1.max() <= py1
    assert len(pm_.edges) >= len([e for e in tpl.edges if "arc" not in e])


def test_pathway_binding_by_kegg_and_generic_nodes():
    A = pd.DataFrame({"Compound": ["Citric acid", "Pyruvic acid", "Some acyl-CoA"],
                      "KEGG": ["C00158", "C00022", "C01144"], "Confidence": ["High", "High", "High"],
                      "Label": ["Citrate", "Pyruvate", ""], "StandardName": ["Citrate", "Pyruvate", ""],
                      "Pathways": ["TCA cycle", "Glycolysis", "Fatty acid metabolism"]})
    b = pmap.bind(ptpl.T["tca"], A, A.Compound)
    assert b["cit"] == "Citric acid" and b["pyr"] == "Pyruvic acid"
    # generic class nodes (3-hydroxyacyl-CoA) bind only to user-confirmed IDs
    assert "hacyl" not in pmap.bind(ptpl.T["fatty_acid"], A, A.Compound)
    A.loc[2, "Confidence"] = "User-confirmed"
    assert pmap.bind(ptpl.T["fatty_acid"], A, A.Compound)["hacyl"] == "Some acyl-CoA"
    # detected members of the annotation pathway that are not on the drawing are listed below it
    assert pmap.extras(ptpl.T["glycolysis"], A, A.Compound, pmap.bind(ptpl.T["glycolysis"], A, A.Compound)) == []
    A.loc[0, "Pathways"] = "Glycolysis"
    assert pmap.extras(ptpl.T["glycolysis"], A, A.Compound, {}) == ["Citric acid", "Pyruvic acid"]


def test_pathway_figures_render():
    tpl = ptpl.T["tca"]
    binding = {"cit": "Citrate", "akg": "aKG", "mal": "Malate"}
    cols = [("A", 1.0), ("B", 1.0)]
    blocks, rows = _blocks_for(binding, cols)
    pm_ = pmap.layout(tpl, binding, rows, cols, ["iso", "z"], extra=[])
    fig = pmap.render_mpl(pm_, blocks, rows, width=10)
    assert abs(fig.get_size_inches()[0] - 10) < 0.05
    import matplotlib.pyplot as plt
    plt.close(fig)
    svg = pmap.render_svg(pm_, blocks, rows)
    assert svg.startswith("<svg") and svg.count("<title>") >= 3 * 2 * 5 * 2       # hover text on every cell
    assert "pmv-stage" in pmap.viewer_html(svg)


def test_m0_only_metabolites_excluded():
    from core.io_reorganize import exclude_m0_only
    W = pd.DataFrame({"Method": ["RP"] * 3, "Selection": [""] * 3, "Compound": ["a", "b", "b"],
                      "Peak": ["1 of 1"] * 3, "Isotopologue": ["M+0", "M+0", "M+1"], "s1": [1.0, 2.0, 3.0]})
    W2, Wr, kept, names = exclude_m0_only(W, W.iloc[0:0], {("RP", "a", "1 of 1"), ("RP", "b", "1 of 1")})
    assert names == ["a"] and list(W2.Compound) == ["b", "b"] and ("RP", "a", "1 of 1") not in kept
    assert len(Wr) == 1 and Wr.Selection.iloc[0].startswith("Removed: only M+0")


def test_low_value_rule_in_generic_route():
    CT = pd.DataFrame({"Compound": ["A"] * 3 + ["B"] * 2, "Formula": "C3H6O3", "Method": "x", "n": 3,
                       "i": [0, 1, 2, 0, 1], "Channel": ["A | M+0", "A | M+1", "A | M+2", "B | M+0", "B | M+1"],
                       "s1": [5e4, 500.0, 3e3, 800.0, 2e3], "s2": [6e4, 2e3, 4e3, 9e4, 3e3]})
    out, log = dsm.apply_zero_rule(CT, ["s1", "s2"], 1000)
    assert list(out.Channel) == ["A | M+0", "A | M+2"] and log


def test_reference_sample_naming():
    from core.io_reorganize import parse_filename, reference_recognises
    assert parse_filename("08302026_RP_Pos_CellLine_GeneKD_1001_12Hr.raw")["Sample"] == "GeneKD_12h_1001"
    assert parse_filename("08302026_RP_Pos_CellLine_shControl_1004_18Hr.raw")["Sample"] == "shCtrl_18h_1004"
    assert parse_filename("08302026_RP_Pos_CellLine_shControl_1009_12C.raw")["Sample"] == "12C_Ref_1009"
    assert reference_recognises(["x_GeneKD_1_12Hr.raw", "x_Blank_3.raw"])
    assert not reference_recognises(["0101_RP_Pos_HeLa_WT_101_30min.raw"])


# ---------------- kinetics ----------------
def test_first_order_recovers_parameters():
    rng = np.random.default_rng(5)
    t = np.repeat([2, 4, 8, 16, 32], 3).astype(float)
    y = 60 * (1 - np.exp(-0.15 * t)) + rng.normal(0, 1, len(t))
    means = pd.Series(y).groupby(t).mean()
    r = kn.fit_one(t, y, means)
    assert r["converged"] and abs(r["k"] - 0.15) < 0.03 and abs(r["E_inf"] - 60) < 4
    assert r["t_half_lo"] < np.log(2) / 0.15 < r["t_half_hi"]


def test_kinetics_refuses_two_timepoints():
    t = np.array([12, 12, 12, 18, 18, 18.0]); y = np.array([10, 11, 9, 20, 21, 19.0])
    r = kn.fit_one(t, y, pd.Series(y).groupby(t).mean())
    assert r["status"] == "Not fitted" and "≥ 3" in r["reason"]


def test_kinetics_refuses_falling_labeling():
    t = np.repeat([12, 18, 40.0], 3); y = np.array([30, 31, 29, 60, 61, 59, 40, 41, 39.0])
    r = kn.fit_one(t, y, pd.Series(y).groupby(t).mean())
    assert r["status"] == "Not fitted" and "falls" in r["reason"]


def test_ranking_requires_supported_fits():
    fits = pd.DataFrame([dict(Compound="a", Group="g", status="Fitted with warnings", t_half=3.0, t_half_lo=1,
                              t_half_hi=9, E_inf=50, R2=0.9)])
    R, msg = kn.speed_ranking(fits)
    assert R.empty and "No ranking" in msg


# ---------------- names, mapping, interpretation ----------------
def test_generic_filename_parsing():
    S = parse_filenames(["0101_RP_Pos_HeLa_WT_101_30min.raw", "0101_RP_Pos_HeLa_KO_102_2h.raw",
                         "0101_RP_Pos_HeLa_WT_103_12C.raw", "0101_RP_Pos_Blank_1.raw"])
    d = S.set_index("Filename")
    assert d.loc["0101_RP_Pos_HeLa_WT_101_30min.raw", "TimeUnit"] == "min"
    assert d.loc["0101_RP_Pos_HeLa_KO_102_2h.raw", "Group"] == "KO"
    assert d.loc["0101_RP_Pos_HeLa_WT_103_12C.raw", "SampleType"] == TYPE_REF
    assert d.loc["0101_RP_Pos_Blank_1.raw", "SampleType"] == TYPE_BLANK


def test_mapping_confidence_rules():
    data = pd.DataFrame({"Compound": ["Citric acid", "D-Glutamine", "3-Oxoglutaric acid", "Palmitoylcarnitine",
                                      "Octenoylcarnitine", "Totally unknown"],
                         "Formula": ["C6 H8 O7", "C5 H10 N2 O3", "C5 H6 O5", "C23 H45 N O4", "C15 H27 N O4", "C99"]})
    A = pw.map_metabolites(data).set_index("Compound")
    assert A.loc["Citric acid", "Confidence"] == "High" and "TCA cycle" in A.loc["Citric acid", "Pathways"]
    assert A.loc["D-Glutamine", "Confidence"] == "Medium"
    assert A.loc["3-Oxoglutaric acid", "Confidence"] == "Low" and A.loc["3-Oxoglutaric acid", "Pathways"] == ""
    assert A.loc["Octenoylcarnitine", "Confidence"] == "Class rule" and A.loc["Octenoylcarnitine", "Pathways"] == ""
    assert A.loc["Totally unknown", "Confidence"] == "Unmapped"


def test_lint_flags_causal_language():
    hits = ip.lint("glycolysis was activated and showed increased flux")
    assert "activated" in hits and "increased flux" in hits
    assert ip.lint("consistent with a greater relative contribution") == []


# ---------------- full pipeline on the synthetic example data ----------------
GEN = ROOT / "example_data" / "generic"


@pytest.mark.skipif(not (GEN / "synthetic_isotopologue_table.csv").exists(), reason="example data not present")
def test_pipeline_on_synthetic_example():
    t, m = GEN / "synthetic_isotopologue_table.csv", GEN / "synthetic_sample_metadata.csv"
    out = dsm.read_generic(str(t), t.name, str(m), m.name)
    d = out[0] if isinstance(out, tuple) else out
    D = dz.finalize_design(dz.default_design(d.S))
    assert dz.groups_of(D) == ["shControl", "GeneKD"] and dz.times_of(D) == [2.0, 8.0, 24.0]
    s = dz.analysed(D).Sample.tolist()
    assert len(s) == 18 and dz.reference_samples(D)
    F, S = lb.compute_from_ct(d.CT, s)
    assert np.allclose(S.SumCheck, 1.0)
    obs, forms = lb.reference_m1_ct(d.CT, dz.reference_samples(D)[0])
    assert cr.check_correction_status(obs, forms)["verdict"] == "corrected"
    lab = [x for x in s if x in d.CT.columns]
    assert not (d.CT[lab] == 0).any().any()
    # synthetic: no sample, group or file name of a real study
    text = t.read_text() + m.read_text()
    assert "HER2" not in text and "RT112" not in text


def test_figure_export_all_formats():
    import matplotlib.pyplot as plt
    from core.style import fig_bytes, EXPORT_FORMATS
    fig, ax = plt.subplots(figsize=(2, 2)); ax.plot([0, 1], [0, 1])
    sigs = {"PNG": b"\x89PNG", "TIFF": (b"II*\x00", b"MM\x00*"), "PDF": b"%PDF", "SVG": b"<?xml", "JPEG": b"\xff\xd8"}
    for f in EXPORT_FORMATS:
        data = fig_bytes(fig, f, 300)
        sig = sigs[f]
        assert any(data.startswith(s) for s in (sig if isinstance(sig, tuple) else (sig,))), f
    plt.close(fig)


def test_zero_rule_deletes_rows_and_compounds():
    CT = pd.DataFrame({"Compound": ["A"] * 3 + ["B"] * 2, "Formula": "C3H6O3", "Method": "x", "n": 3,
                       "i": [0, 1, 2, 0, 1], "Channel": "", "s1": [10.0, 0.0, 3.0, 0.0, 4.0], "s2": [9.0, 2.0, 1.0, 5, 3],
                       "ref": [0.0, 0.0, 0.0, 0.0, 0.0]})
    out, log = dsm.apply_zero_rule(CT, ["s1", "s2"])
    assert list(zip(out.Compound, out.i)) == [("A", 0), ("A", 2)]       # A M+1 row removed, B removed (M+0 = 0)
    assert any("M+0 = 0" in x for x in log)


def test_two_group_tests_match_scipy():
    from scipy import stats
    rng = np.random.default_rng(3)
    X = pd.DataFrame(rng.normal(0, 1, (4, 6)), columns=list("abcdef"))
    T = stt.two_group(X, ["a", "b", "c"], ["d", "e", "f"], "student")
    assert np.isclose(T.p[0], stats.ttest_ind(X.loc[0, ["d", "e", "f"]], X.loc[0, ["a", "b", "c"]]).pvalue)
    T = stt.two_group(X, [], [], "paired", pairs=[("a", "d"), ("b", "e"), ("c", "f")])
    assert np.isclose(T.p[1], stats.ttest_rel(X.loc[1, ["d", "e", "f"]], X.loc[1, ["a", "b", "c"]]).pvalue)
    assert np.allclose(T.FDR, stt.bh_fdr(T.p), equal_nan=True)
    A, _ = stt.anova_method(X, {"g1": ["a", "b"], "g2": ["c", "d"], "g3": ["e", "f"]}, "kruskal", "none")
    assert np.isclose(A.p[0], stats.kruskal(X.loc[0, ["a", "b"]], X.loc[0, ["c", "d"]], X.loc[0, ["e", "f"]]).pvalue)


# ---------------- v2.4.2: correction methods, pathway cards ----------------
@pytest.mark.parametrize("method", list(cr.METHODS))
def test_correction_methods_round_trip(method):
    true = np.array([0.4, 0.05, 0.3, 0.05, 0.05, 0.05, 0.1])
    kw = dict(resolving_power=7e4, instrument="orbitrap") if method == "accucor" else {}
    A = cr.correction_matrix("C6H8O7", 7, purity=0.99, method=method, **kw)
    assert np.allclose(cr.correct_mid(A @ true, "C6H8O7", purity=0.99, method=method, **kw), true, atol=1e-8)


def test_accucor_resolution_dependence():
    # 15N in ATP is unresolved from 13C at 70k but resolved at very high resolving power
    lo = cr.unresolved_isotopes("C10H16N5O13P3", resolving_power=7e4)
    hi = cr.unresolved_isotopes("C10H16N5O13P3", resolving_power=1e6)
    assert lo["N"] and not hi["N"]
    # at very high resolution AccuCor gives the same corrected MID as tracer-only correction (C/H/N compound)
    m = cr.correction_matrix("C5H5N5", 6, method="isocor_lr") @ np.array([0.5, 0.1, 0.1, 0.1, 0.1, 0.1])
    a = cr.correct_mid(m, "C5H5N5", method="accucor", resolving_power=1e9)
    b = cr.correct_mid(m, "C5H5N5", method="isocor_hr")
    assert np.allclose(a, b, atol=1e-9)


# ---------------- v2.5.3: Compound Discoverer "Exchange Rate in Percent" denominator ----------------
def test_fraction_denominator_uses_every_isotopologue(tmp_path):
    """Fraction = Area(M+i) / Σ Area(M+0…M+n) over ALL isotopologues — also rows removed by the minimum-value rule."""
    areas = [155650000, 46659402, 30864767, 9727657, 2479805, 783031]
    cd = [63.23003408, 18.95454932, 12.53826061, 3.951686897, 1.007376696, 0.318092404]   # CD output
    T = pd.DataFrame({"Metabolite": "X", "Formula": "C5H9NO4", "Method": "HILIC_Neg",
                      "Isotopologue": [f"M+{i}" for i in range(6)], "s1": areas, "s2": areas})
    T.loc[5, "s2"] = 500.0                              # below the minimum value in one labeled sample
    M = pd.DataFrame({"Sample": ["s1", "s2"], "Type": "Labeled", "Group": "A", "Time": 1, "TimeUnit": "h"})
    (tmp_path / "t.csv").write_text(T.to_csv(index=False)); (tmp_path / "m.csv").write_text(M.to_csv(index=False))
    out = dsm.read_generic(str(tmp_path / "t.csv"), "t.csv", str(tmp_path / "m.csv"), "m.csv", min_area=1000)
    d = out[0] if isinstance(out, tuple) else out
    assert 5 not in set(d.CT.i)                         # M+5 removed from the intensity analyses …
    assert sorted(d.labeling_table.i) == list(range(6))  # … but kept for labeling
    F, S = lb.compute_from_ct(d.labeling_table, ["s1"])
    assert np.allclose(F.sort_values("i").Fraction * 100, cd, atol=1e-6)
    s1 = S.iloc[0]
    assert np.isclose(s1.Incorporation, 100 - cd[0], atol=1e-6)
    assert np.isclose(s1.MeanLabels, sum(i * f / 100 for i, f in enumerate(cd)), atol=1e-6)
    assert np.isclose(s1.Enrichment, s1.MeanLabels / 5 * 100)


# ---------------------------------------------------------------------------
# v2.5.7 heatmap: grouped row labels, proportional cells
# ---------------------------------------------------------------------------
def test_heatmap_grouped_rows_and_balanced_cells():
    from core import heatmap as hm
    rows = [f"{m} | M+{i}" for m in ["B acid", "A acid", "C acid"] for i in (2, 0, 1)]
    Z = pd.DataFrame(np.random.default_rng(0).normal(size=(9, 6)), index=rows)
    order = hm.grouped_order(Z, cluster=False)
    assert order[:3] == ["A acid | M+0", "A acid | M+1", "A acid | M+2"]
    order_c = hm.grouped_order(Z, cluster=True)
    comps = [r.split(" | ")[0] for r in order_c]
    assert all(comps[i] == comps[i + 1] for i in (0, 1, 3, 4, 6, 7))       # blocks stay contiguous
    assert hm.row_runs(comps) == [(comps[0], 0, 2), (comps[3], 3, 5), (comps[6], 6, 8)]
    txt = hm.channel_text(order, ns={"A acid | M+1"})
    assert txt["A acid | M+1"] == "M+1 (n.s.)" and txt["A acid | M+0"] == "M+0"
    for nr, nc in [(5, 6), (20, 12), (400, 36)]:
        cw, ch = hm.balanced_cells(nr, nc)
        assert 0.11 <= ch <= 0.22 and 0.14 <= cw <= 0.45
        assert cw * nc <= 7.0 + 1e-9 or cw == 0.14
    f, r_, _ = hm.draw_heatmap(Z.loc[order], pd.DataFrame({"Group": list("aabbcc")}, index=Z.columns),
                               cluster_rows=False, row_text=txt, row_groups={r: r.split(" | ")[0] for r in order},
                               annotations=["Group"], cell_w=0.2, cell_h_in=0.15)
    assert r_ == order
    import matplotlib.pyplot as plt
    plt.close(f)


def test_verify_online_subset_offline(tmp_path, monkeypatch):
    import requests
    from core import pathways as pw
    A = pw.map_metabolites(pd.DataFrame({"Compound": ["Citrate", "Lactate"], "Formula": ["C6H8O7", "C3H6O3"]}))
    def boom(*a, **k):
        raise requests.ConnectionError("offline")
    monkeypatch.setattr(requests, "get", boom)
    cache = tmp_path / "c.json"
    out = pw.verify_online(A, cache_file=cache, compounds=["Citrate"])
    v = out.set_index("Compound").Verification
    assert "lookup failed" in v["Citrate"] and not str(v.get("Lactate", "") or "")
    assert cache.read_text().strip() == "{}"                 # failures are not cached (retried later)
    assert pw.online_available() is False


def test_heatmap_pathway_sections():
    from core import heatmap as hm
    secs = ["A", "A", "B", "B", "B", "C"]
    assert hm.row_positions(secs, 0.5) == [0, 1, 2.5, 3.5, 4.5, 6.0]
    # same metabolite in two adjacent pathways -> two separate runs
    assert hm.keyed_runs(["Glu", "Glu", "Glu", "Gln", "Gln", "Cit"], secs) == \
        [("Glu", 0, 1), ("Glu", 2, 2), ("Gln", 3, 4), ("Cit", 5, 5)]
    fam = hm._families("Arial")
    lines, fs = hm.wrap_to_length("Glutamate metabolism", 8, fam, 0.6)
    assert len(lines) == 2 and fs <= 8
    assert hm.wrap_to_length("TCA", 8, fam, 2.0) == (["TCA"], 8)
    rows = [f"{p} ‖ X{i} | M+{j}" for i, p in enumerate(["P1", "P2"]) for j in range(3)]
    Z = pd.DataFrame(np.random.default_rng(1).normal(size=(6, 4)), index=rows)
    f, order, _ = hm.draw_heatmap(Z, pd.DataFrame({"Group": list("aabb")}, index=Z.columns), cluster_rows=False,
                                  row_sections={r: r.split(" ‖ ")[0] for r in rows}, section_gap=1.0,
                                  row_groups={r: r.split(" ‖ ")[1].split(" | ")[0] for r in rows},
                                  annotations=["Group"], cell_w=0.2, cell_h_in=0.15, title="t")
    assert order == rows
    import matplotlib.pyplot as plt
    plt.close(f)


# ---------------------------------------------------------------------------
# v2.5.11 templates
# ---------------------------------------------------------------------------
def test_templates_load():
    from pathlib import Path
    from core import io_reorganize as io, dataset as dsm
    T = Path(__file__).resolve().parent.parent / "example_data" / "templates"
    r = io.process_files([T / "cd_sil" / "CD_SIL_Template_RP_Pos.xlsx", T / "cd_sil" / "CD_SIL_Template_HILIC_Neg.xlsx"])
    assert r.naming == "reference" and set(m for m, _ in r.sources) == {"RP_Pos", "HILIC_Neg"}
    assert (r.S.SampleType == io.TYPE_LAB).sum() == 18 and (r.S.SampleType == io.TYPE_REF).sum() == 1
    assert r.W.Compound.nunique() > 50
    g = T / "generic"
    meta = g / "Generic_Sample_Metadata.csv"
    d1 = dsm.read_generic(str(g / "Generic_Isotopologue_Table_with_12C_Ref.csv"), "t.csv", str(meta), "m.csv")
    d2 = dsm.read_generic(str(g / "Generic_Isotopologue_Table_without_12C_Ref.csv"), "t.csv", str(meta), "m.csv")
    assert len(d1.S) == 19 and len(d2.S) == 18 and d1.CT.Compound.nunique() == 63
