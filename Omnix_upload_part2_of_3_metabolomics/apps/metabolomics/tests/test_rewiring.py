"""
Tests for metabolomics_modules/rewiring.py (Correlation Rewiring Map engine).

    python -m pytest tests/test_rewiring.py -q

Groups:
  * statistical correctness against numpy / scipy / statsmodels
  * invariances the method must respect
  * recovery of planted rewiring and false-discovery control (simulation)
  * calibration under the null, including small groups
  * hypothesis-card checks (outlier decoy, imputation artifact, batch)
  * input handling and serialization
"""

import json
import math
import os
import sys

import numpy as np
import pandas as pd
import pytest
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from metabolomics_modules import rewiring as rw  # noqa: E402


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _demo_prep(**sim_kwargs):
    d = rw.simulate_rewiring_demo(**sim_kwargs)
    prep = rw.prepare_data(d["peak_df"], d["meta"], "Group", "Control", "Tumor",
                           pathways=d["annotations"].set_index("Metabolite")["Pathway"], batch_col="Batch")
    return d, prep


def _found(res):
    return {frozenset((a, b)): c for a, b, c in zip(res.pairs.feature_a, res.pairs.feature_b, res.pairs["class"])
            if c in rw.REWIRED}


def _gaussian_frame(rng, m, n_a, n_b, loading=0.7, noise=0.6):
    L = rng.standard_normal((m, 5)) * loading
    X = L @ rng.standard_normal((5, n_a + n_b)) + rng.standard_normal((m, n_a + n_b)) * noise + 15
    df = pd.DataFrame(X, index=[f"f{i}" for i in range(m)], columns=[f"s{j}" for j in range(n_a + n_b)])
    meta = pd.DataFrame({"Group": ["A"] * n_a + ["B"] * n_b}, index=df.columns)
    return df, meta


@pytest.fixture(scope="module")
def demo_result():
    d, prep = _demo_prep()
    res = rw.run_rewiring(prep, rw.RewiringConfig(n_permutations=200, n_bootstrap=200))
    return d, prep, res


# ---------------------------------------------------------------------------
# statistical correctness
# ---------------------------------------------------------------------------
def test_correlations_match_numpy_and_scipy():
    rng = np.random.default_rng(1)
    X = rng.standard_normal((12, 25))
    idx = np.arange(25)
    np.testing.assert_allclose(rw.correlation_matrix(X, idx, "pearson"), np.corrcoef(X), atol=1e-10)
    rho = stats.spearmanr(X.T).statistic
    np.testing.assert_allclose(rw.correlation_matrix(X, idx, "spearman"), rho, atol=1e-10)


def test_fisher_z_and_bh_match_reference(demo_result):
    _, prep, res = demo_result
    row = res.pairs.iloc[123]
    xa = prep.X[row.i, prep.idx_a], prep.X[row.j, prep.idx_a]
    xb = prep.X[row.i, prep.idx_b], prep.X[row.j, prep.idx_b]
    ra, rb = stats.pearsonr(*xa).statistic, stats.pearsonr(*xb).statistic
    z = (np.arctanh(rb) - np.arctanh(ra)) / math.sqrt(1 / (prep.n_a - 3) + 1 / (prep.n_b - 3))
    assert row.r_a == pytest.approx(ra, abs=1e-10)
    assert row.r_b == pytest.approx(rb, abs=1e-10)
    assert row.z == pytest.approx(z, abs=1e-8)
    assert row.p_value == pytest.approx(2 * stats.norm.sf(abs(z)), rel=1e-8)
    q_ref = multipletests(res.pairs["p_value"].values, method="fdr_bh")[1]
    np.testing.assert_allclose(res.pairs["q_bh"].values, q_ref, rtol=1e-12)


def test_spearman_uses_inflated_variance():
    assert rw.fisher_se(30, 30, "spearman") == pytest.approx(math.sqrt(1.06) * rw.fisher_se(30, 30, "pearson"))


def test_welch_fold_change_matches_scipy(demo_result):
    _, prep, res = demo_result
    k = 5
    t, p = stats.ttest_ind(prep.X[k, prep.idx_b], prep.X[k, prep.idx_a], equal_var=False)
    assert res.features["t"].iat[k] == pytest.approx(t)
    assert res.features["p_value"].iat[k] == pytest.approx(p)
    assert res.features["log2fc"].iat[k] == pytest.approx(prep.X[k, prep.idx_b].mean() - prep.X[k, prep.idx_a].mean())


def test_permutation_fdr_matches_bruteforce_definition():
    rng = np.random.default_rng(3)
    obs = np.abs(rng.standard_normal(60)) * 1.3
    obs[:5] += 4
    n_perm = 7
    null = np.sort(np.abs(rng.standard_normal(60 * n_perm)).astype(np.float32))
    q = rw.permutation_fdr(obs, null, n_perm)
    for i in range(len(obs)):
        best = 1.0
        for t in obs[obs <= obs[i]]:
            fdr = min(1.0, (max(np.sum(null >= np.float32(t)), 1) / n_perm) / np.sum(obs >= t))
            best = min(best, fdr)
        assert q[i] == pytest.approx(best, abs=1e-12)
    order = np.argsort(obs)
    assert np.all(np.diff(q[order]) <= 1e-12), "q must not increase with |z|"
    assert q.min() > 0


def test_classification_rules():
    cfg = rw.RewiringConfig(fdr=0.05, min_delta=0.45)
    pairs = pd.DataFrame({
        "r_a": [0.8, 0.0, 0.7, 0.8, 0.7, 0.2],
        "r_b": [0.0, 0.8, -0.7, 0.75, 0.1, 0.9],
        "q_value": [0.01, 0.01, 0.01, 0.9, 0.2, 0.01],
    })
    pairs["delta_r"] = pairs["r_b"] - pairs["r_a"]
    out = [None if pd.isna(c) else c for c in rw.classify_pairs(pairs, cfg)["class"]]
    assert out == ["lost", "gained", "flipped", "intact", None, "gained"]


# ---------------------------------------------------------------------------
# invariances
# ---------------------------------------------------------------------------
def test_invariant_to_feature_and_sample_order():
    d = rw.simulate_rewiring_demo(n_per_group=20, seed=4)
    peak, meta = d["peak_df"], d["meta"]
    pw = d["annotations"].set_index("Metabolite")["Pathway"]
    r1 = rw.run_rewiring(rw.prepare_data(peak, meta, "Group", "Control", "Tumor", pathways=pw),
                         rw.RewiringConfig(fdr_method="bh", n_permutations=0, n_bootstrap=20))
    rng = np.random.default_rng(0)
    peak2 = peak.iloc[rng.permutation(len(peak)), rng.permutation(peak.shape[1])]
    r2 = rw.run_rewiring(rw.prepare_data(peak2, meta, "Group", "Control", "Tumor", pathways=pw),
                         rw.RewiringConfig(fdr_method="bh", n_permutations=0, n_bootstrap=20))
    key = lambda df: df.assign(k=[frozenset(x) for x in zip(df.feature_a, df.feature_b)]).set_index("k")  # noqa: E731
    a, b = key(r1.pairs), key(r2.pairs)
    b = b.loc[a.index]
    flip = a.feature_a.values != b.feature_a.values
    np.testing.assert_allclose(a.r_a, b.r_a, atol=1e-10)
    np.testing.assert_allclose(a.q_bh, b.q_bh, atol=1e-10)
    assert (a["class"].fillna("") == b["class"].fillna("")).all()
    assert flip.any()   # the shuffle really changed pair orientation


def test_group_mean_shift_does_not_change_correlations():
    rng = np.random.default_rng(2)
    df, meta = _gaussian_frame(rng, 20, 25, 25)
    prep1 = rw.prepare_data(df, meta, "Group", "A", "B", log_transform="no")
    shifted = df.copy()
    shifted.loc[:, meta.Group == "B"] += rng.normal(0, 3, size=(20, 1))    # per-feature fold changes in B
    prep2 = rw.prepare_data(shifted, meta, "Group", "A", "B", log_transform="no")
    p1, p2 = rw.differential_correlation(prep1), rw.differential_correlation(prep2)
    np.testing.assert_allclose(p1.z, p2.z, atol=1e-9)


# ---------------------------------------------------------------------------
# recovery of planted rewiring
# ---------------------------------------------------------------------------
def test_demo_recovers_planted_rewiring(demo_result):
    d, _, res = demo_result
    t, found = d["truth"], _found(res)
    recall = {k: np.mean([p in found for p in t[k]]) for k in ("lost", "gained", "flipped")}
    assert recall["flipped"] == 1.0
    assert recall["lost"] >= 0.7
    assert recall["gained"] >= 0.6
    for k in ("lost", "gained", "flipped"):
        assert all(found[p] == k for p in t[k] if p in found), f"misclassified {k} pair"
    allowed = t["lost"] | t["gained"] | t["flipped"] | t["ambiguous"] | t["decoy"]
    assert [tuple(p) for p in found if p not in allowed] == []
    assert res.permutation["p_value"] < 0.01


def test_false_discovery_proportion_across_seeds():
    fdps, recalls = [], []
    for seed in range(10):
        d, prep = _demo_prep(seed=seed)
        res = rw.run_rewiring(prep, rw.RewiringConfig(n_permutations=100, n_bootstrap=20, seed=seed))
        t, found = d["truth"], _found(res)
        true_any = t["lost"] | t["gained"] | t["flipped"]
        disc = [p for p in found if p not in t["ambiguous"] and p not in t["decoy"]]
        fdps.append(np.mean([p not in true_any for p in disc]) if disc else 0.0)
        recalls.append(np.mean([p in found for p in true_any]))
    assert np.mean(fdps) <= 0.08, fdps
    assert np.mean(recalls) >= 0.7, recalls


# ---------------------------------------------------------------------------
# calibration under the null
# ---------------------------------------------------------------------------
def test_null_n30_no_discoveries_and_uniformish_global_p():
    any_disc, pvals = [], []
    for seed in range(100, 115):
        d, prep = _demo_prep(seed=seed, rewire=False, outliers=False)
        res = rw.run_rewiring(prep, rw.RewiringConfig(n_permutations=100, n_bootstrap=10, seed=seed))
        any_disc.append(sum(res.counts()[k] for k in rw.REWIRED) > 0)
        pvals.append(res.permutation["p_value"])
    assert np.mean(any_disc) <= 0.2
    assert np.mean(np.array(pvals) < 0.05) <= 0.2
    assert np.median(pvals) > 0.15


def test_small_groups_permutation_fdr_is_calibrated_and_beats_bh():
    rng = np.random.default_rng(0)
    hits = {"bh": [], "permutation": []}
    for rep in range(20):
        df, meta = _gaussian_frame(rng, 120, 9, 9)
        prep = rw.prepare_data(df, meta, "Group", "A", "B", log_transform="no")
        for fm in hits:
            res = rw.run_rewiring(prep, rw.RewiringConfig(fdr_method=fm, n_permutations=100, n_bootstrap=5, seed=rep))
            hits[fm].append(sum(res.counts()[k] for k in rw.REWIRED) > 0)
    rate_perm, rate_bh = np.mean(hits["permutation"]), np.mean(hits["bh"])
    assert rate_perm <= 0.15, rate_perm
    assert rate_perm <= rate_bh


# ---------------------------------------------------------------------------
# hypothesis cards
# ---------------------------------------------------------------------------
def test_true_modules_are_robust_with_curated_text(demo_result):
    _, _, res = demo_result
    by_key = {h["key"]: h for h in res.hypotheses}
    gln = by_key["Glutamine / Amino Acid Metabolism|TCA Cycle|lost"]
    trp = by_key["Tryptophan–Kynurenine Metabolism|Tryptophan–Kynurenine Metabolism|flipped"]
    for h in (gln, trp):
        assert h["verdict"] == "robust"
        assert h["curated_interpretation"]
        assert h["stability"] >= 0.9
        assert h["batch_check"]["ok"] is True and len(h["batch_check"]["levels"]) == 2
    assert "Citrate" in gln["hidden_features"] or "Aspartate" in gln["hidden_features"]
    assert res.hypotheses[0]["rank"] == 1
    assert all(res.hypotheses[i]["score"] >= res.hypotheses[i + 1]["score"] for i in range(len(res.hypotheses) - 1))


def test_outlier_decoy_is_flagged(demo_result):
    d, _, res = demo_result
    decoy = [h for h in res.hypotheses if any(frozenset(p) in d["truth"]["decoy"] for p in h["pair_names"])]
    assert decoy, "decoy should be called at the pair level in the default demo"
    h = decoy[0]
    assert h["verdict"] != "robust"
    assert h["influence_check"]["ok"] is False
    assert set(h["influence_check"]["dropped"]) >= {"TUM_05", "TUM_18", "TUM_26"}


def test_decoy_flagged_across_seeds():
    flagged = []
    for seed in range(8):
        d, prep = _demo_prep(seed=seed)
        res = rw.run_rewiring(prep, rw.RewiringConfig(n_permutations=100, n_bootstrap=50, seed=seed))
        for h in res.hypotheses:
            if any(frozenset(p) in d["truth"]["decoy"] for p in h["pair_names"]):
                flagged.append(h["verdict"] != "robust")
    assert flagged and all(flagged)


def test_imputation_artifact_is_flagged():
    """Two independent metabolites that are 'not detected' in the same 10 samples of one
    group get the same imputed floor value there, which fakes a strong correlation."""
    rng = np.random.default_rng(7)
    X = rng.standard_normal((12, 50)) + 15
    df = pd.DataFrame(2.0 ** X, index=[f"f{i}" for i in range(12)], columns=[f"s{j}" for j in range(50)])
    meta = pd.DataFrame({"Group": ["A"] * 25 + ["B"] * 25}, index=df.columns)
    nd = rng.choice(meta.index[meta.Group == "B"], 10, replace=False)
    df.loc[["f0", "f1"], nd] = 0.0
    prep = rw.prepare_data(df, meta, "Group", "A", "B", pathways=pd.Series("Pathway X", index=df.index),
                           max_missing_frac=0.5)
    assert prep.imputed.sum() == 20
    res = rw.run_rewiring(prep, rw.RewiringConfig(fdr_method="bh", n_permutations=0, n_bootstrap=50, min_delta=0.3))
    pair = res.pairs[(res.pairs.feature_a == "f0") & (res.pairs.feature_b == "f1")].iloc[0]
    assert pair.r_b > 0.75 and pair["class"] == "gained"
    h = next(h for h in res.hypotheses if ["f0", "f1"] in h["pair_names"])
    assert h["imputation_check"]["ok"] is False and h["imputation_check"]["n_imputed"] == 20
    assert h["verdict"] != "robust"
    assert "imputation artifact" in h["interpretation"]


def test_median_centering_removes_loading_artifact():
    rng = np.random.default_rng(11)
    m, n = 30, 30
    X = rng.standard_normal((m, 2 * n)) * 0.5 + 15
    X[:, n:] += rng.normal(0, 1.0, size=(1, n))      # group B: strong per-sample dilution differences
    df = pd.DataFrame(X, index=[f"f{i}" for i in range(m)], columns=[f"s{j}" for j in range(2 * n)])
    meta = pd.DataFrame({"Group": ["A"] * n + ["B"] * n}, index=df.columns)
    cfg = rw.RewiringConfig(fdr_method="bh", n_permutations=0, n_bootstrap=10)
    raw = rw.run_rewiring(rw.prepare_data(df, meta, "Group", "A", "B", log_transform="no"), cfg)
    cen = rw.run_rewiring(rw.prepare_data(df, meta, "Group", "A", "B", log_transform="no", sample_median_center=True), cfg)
    assert raw.counts()["gained"] > 100
    assert cen.counts()["gained"] + cen.counts()["lost"] <= 5


# ---------------------------------------------------------------------------
# input handling
# ---------------------------------------------------------------------------
def test_missing_values_zero_variance_and_qc_handling():
    d = rw.simulate_rewiring_demo(n_per_group=20, seed=1)
    peak, meta = d["peak_df"].copy(), d["meta"].copy()
    peak.iloc[0, :15] = np.nan                         # >30% missing in Control -> dropped
    peak.iloc[1, [2, 25]] = np.nan                     # imputed
    peak.iloc[2, :] = 1234.5                           # constant -> dropped
    peak["QC_1"] = peak.mean(axis=1)
    meta = pd.concat([meta, pd.DataFrame({"Sample": ["QC_1"], "Group": ["Control"], "IsQC": [True], "Batch": ["B1"]})])
    prep = rw.prepare_data(peak, meta, "Group", "Control", "Tumor")
    assert "QC_1" not in prep.samples
    assert peak.index[0] not in prep.features and peak.index[2] not in prep.features
    assert peak.index[1] in prep.features and prep.imputed[prep.features.index(peak.index[1])].sum() == 2
    assert not np.isnan(prep.X).any()
    assert prep.logged is True
    joined = " ".join(prep.notes)
    assert "dropped for >30% missing" in joined and "imputed" in joined and "no variation" in joined


def test_input_errors_are_user_facing():
    d = rw.simulate_rewiring_demo(n_per_group=10, seed=1)
    peak, meta = d["peak_df"], d["meta"]
    with pytest.raises(rw.RewiringInputError, match="two different groups"):
        rw.prepare_data(peak, meta, "Group", "Control", "Control")
    with pytest.raises(rw.RewiringInputError, match="No samples found for group 'Nope'"):
        rw.prepare_data(peak, meta, "Group", "Control", "Nope")
    with pytest.raises(rw.RewiringInputError, match="no column"):
        rw.prepare_data(peak, meta, "Diagnosis", "Control", "Tumor")
    with pytest.raises(rw.RewiringInputError, match="at least 12 samples"):
        rw.prepare_data(peak, meta, "Group", "Control", "Tumor", min_per_group=12)
    renamed = peak.rename(columns=lambda c: "x" + c)
    with pytest.raises(rw.RewiringInputError, match="None of the data columns"):
        rw.prepare_data(renamed, meta, "Group", "Control", "Tumor")
    with pytest.raises(rw.RewiringInputError):
        rw.RewiringConfig(method="kendall")


def test_feature_cap_keeps_most_variable():
    rng = np.random.default_rng(5)
    df, meta = _gaussian_frame(rng, 50, 15, 15)
    df.iloc[:10] *= 5                                    # 10 high-variance features
    prep = rw.prepare_data(df, meta, "Group", "A", "B", log_transform="no", max_features=10)
    assert prep.features == [f"f{i}" for i in range(10)]
    assert any("10 most variable" in n for n in prep.notes)


def test_pathway_matching_prefix_case_and_clustering():
    d = rw.simulate_rewiring_demo(n_per_group=15, seed=2)
    peak = d["peak_df"].copy()
    peak.index = ["Targeted Metabolomics::" + n for n in peak.index]      # multi-dataset names
    ann = d["annotations"].set_index("Metabolite")["Pathway"]
    ann.index = [n.upper() if i % 2 else n for i, n in enumerate(ann.index)]
    prep = rw.prepare_data(peak, d["meta"], "Group", "Control", "Tumor", pathways=ann)
    assert prep.pathway_source == "annotation" and rw.UNANNOTATED not in prep.pathways
    partial = ann.iloc[:20]
    prep2 = rw.prepare_data(d["peak_df"], d["meta"], "Group", "Control", "Tumor", pathways=partial)
    assert prep2.pathway_source == "mixed" and prep2.pathways.count(rw.UNANNOTATED) == len(prep2.features) - 20
    prep3 = rw.prepare_data(d["peak_df"], d["meta"], "Group", "Control", "Tumor", pathways=None)
    assert prep3.pathway_source == "cluster" and all(p.startswith("Cluster ") for p in prep3.pathways)
    assert len(set(prep3.pathways)) >= 2


def test_auto_log_detection():
    d = rw.simulate_rewiring_demo(n_per_group=10, seed=3)
    assert rw.prepare_data(d["peak_df"], d["meta"], "Group", "Control", "Tumor").logged is True
    logged = np.log2(d["peak_df"])
    assert rw.prepare_data(logged, d["meta"], "Group", "Control", "Tumor").logged is False


def test_pathway_tags_and_interpretation_fallback():
    assert rw.pathway_tag("Citrate cycle (TCA cycle)") == "tca"
    assert rw.pathway_tag("Alanine, aspartate and glutamate metabolism") == "glutamine"
    assert rw.pathway_tag("Tryptophan Metabolism") == "tryptophan"
    assert rw.pathway_tag("Glycolysis / Gluconeogenesis") == "glycolysis"
    assert rw.pathway_tag("Cluster 3") is None
    text, exp, curated = rw._interpret("Cluster 1", "Cluster 2", "gained")
    assert not curated and "independent cohort" in exp
    text, exp, curated = rw._interpret("Purine Metabolism", "Cluster 2", "lost")
    assert not curated and "¹⁵N" in exp


# ---------------------------------------------------------------------------
# serialization
# ---------------------------------------------------------------------------
def test_payload_is_strict_json_and_complete(demo_result):
    _, prep, res = demo_result
    payload = res.to_payload()
    text = json.dumps(payload, allow_nan=False)
    back = json.loads(text)
    assert len(back["features"]) == len(prep.features)
    rew = [p for p in back["pairs"] if p["cls"] != "intact"]
    assert len(rew) == sum(res.counts()[k] for k in rw.REWIRED)
    need = {str(p["i"]) for p in rew} | {str(p["j"]) for p in rew}
    assert need <= set(back["values"])
    assert len(back["hypotheses"]) <= 6
    assert any(h["verdict"] != "robust" for h in back["hypotheses"]), "a flagged card is always surfaced"
    json.loads(res.evidence_json())


def test_view_html_escapes_hostile_names():
    from metabolomics_modules import rewiring_ui
    d = rw.simulate_rewiring_demo(n_per_group=15, seed=2)
    peak = d["peak_df"].copy()
    peak.index = ["</script><img src=x onerror=alert(1)>" if i == 0 else n for i, n in enumerate(peak.index)]
    prep = rw.prepare_data(peak, d["meta"], "Group", "Control", "Tumor")
    res = rw.run_rewiring(prep, rw.RewiringConfig(n_permutations=20, n_bootstrap=10))
    html = rewiring_ui.render_view_html(res)
    body = html.split('<script id="payload" type="application/json">', 1)[1].split("</script>", 1)[0]
    assert "</script" not in body.lower()
    data = json.loads(body)
    assert data["features"][0]["name"].startswith("</script>")
