"""
Tests for the reaction library, MSEA-library mapping, reaction-scoped testing, the
pathway-level (set) test and small-sample (3 per group) behaviour of the Rewiring Map.

    python -m pytest tests/test_rewiring_reactions.py -q
"""

import itertools
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
from metabolomics_modules import rewiring as rw, rewiring_reactions as rr  # noqa: E402
import validate_reaction_library as vrl  # noqa: E402


@pytest.fixture(scope="module")
def lib():
    return rr.reaction_library()


def _gem(names):
    fm = rr.map_features(names, None, library=rr.LIB_CLUSTER)
    return list(fm.table["GEM_ID"])


# ---------------------------------------------------------------------------
# reaction library
# ---------------------------------------------------------------------------
def test_library_passes_gold_standard_and_negative_controls(lib):
    gold, unresolved, _ = vrl.check_gold(lib)
    assert not unresolved
    assert (gold.result == "ok").all(), gold[gold.result != "ok"].to_string()
    neg = vrl.check_negatives(lib)
    assert (neg.result == "ok").all(), neg[neg.result != "ok"].to_string()
    bad, ctx = vrl.check_currency(lib)
    assert bad == []
    assert all(r <= 45 for _, _, r in ctx), "context cofactors must not become two-step hubs"


def test_library_provenance(lib):
    info = lib.info
    assert "Human-GEM" in info["source"] and "CC BY 4.0" in info["source"]
    assert info["n_pairs"] > 4000 and info["version"]


def test_direction_and_two_step_route(lib):
    g = _gem(["Glutamine", "Glutamate", "Tryptophan", "Kynurenine", "Aspartate", "Oxaloacetate"])
    L = lib.links(g, max_steps=2)
    gl = L[(0, 1)]
    assert gl["steps"] == 1 and gl["direction"] in ("i>j", "both")
    tk = L[(2, 3)]
    assert tk["steps"] == 2 and tk["direction"] == "i>j" and "formylkynurenine" in tk["via_name"].lower()
    assert L[(4, 5)]["steps"] == 1                     # transamination main pair
    assert (4, 1) not in L and (1, 4) not in L          # aspartate-glutamate: not a substrate/product pair


def test_peptides_do_not_create_amino_acid_shortcuts(lib):
    g = _gem(["Leucine", "Tryptophan", "Valine", "Histidine", "Phenylalanine"])
    assert lib.links(g, max_steps=2) == {}


def test_name_resolution_layers(lib):
    fm = rr.map_features(["3-Hydroxykynurenine", "Quinolinate", "Glutathione (GSH)", "GSSG", "GABA",
                          "5-Hydroxy-L-tryptophan", "N-Formylkynurenine"], None, library="KEGG")
    t = fm.table
    assert (t["GEM_ID"] != "").all()
    assert (t["Match_Method"] != "Unmapped").all()     # IDs carried back for MSEA-library membership
    assert fm.pathway[0] == "Tryptophan metabolism" and fm.pathway[1] == "Tryptophan metabolism"


def test_kegg_anomer_groups_place_sugar_phosphates_in_glycolysis():
    fm = rr.map_features(["Glucose", "Glucose-6-phosphate", "Fructose-6-phosphate", "Fructose-1,6-bisphosphate",
                          "Pyruvate", "Lactate"], None, library="KEGG")
    gly = fm.sets["Glycolysis / Gluconeogenesis"]
    assert {0, 1, 2, 3, 4, 5} <= set(gly)


def test_prefixed_and_duplicated_names_map():
    fm = rr.map_features(["Targeted Metabolomics::Glutamine", "Untargeted Metabolomics::L-Glutamine"], None,
                         library="KEGG")
    assert fm.table["GEM_ID"].iloc[0] == fm.table["GEM_ID"].iloc[1] != ""


# ---------------------------------------------------------------------------
# engine: reaction scope, ratio shift, pathway-level test
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def demo():
    d = rw.simulate_rewiring_demo()
    prep = rw.prepare_data(d["peak_df"], d["meta"], "Group", "Control", "Tumor", batch_col="Batch",
                           pathway_library="KEGG")
    return d, prep


def test_reaction_scope_restricts_testing_and_fdr(demo):
    _, prep = demo
    res = rw.run_rewiring(prep, rw.RewiringConfig(pair_scope="reaction", n_permutations=100, n_bootstrap=30,
                                                  fdr_method="bh"))
    p = res.pairs
    tested = p["tested"].values
    assert tested.sum() == ((p["rxn_steps"] >= 1) & (p["rxn_steps"] <= 2)).sum() > 0
    assert p.loc[~tested, "q_value"].isna().all()
    expected = multipletests(p.loc[tested, "p_value"], method="fdr_bh")[1]
    assert np.allclose(p.loc[tested, "q_bh"], expected)
    assert set(p.loc[p["class"].isin(rw.REWIRED), "rxn_steps"]) <= {1, 2}


def test_direct_scope(demo):
    _, prep = demo
    res = rw.run_rewiring(prep, rw.RewiringConfig(pair_scope="direct", n_permutations=50, n_bootstrap=20))
    assert (res.pairs.loc[res.pairs["tested"], "rxn_steps"] == 1).all()


def test_ratio_shift_matches_manual(demo):
    _, prep = demo
    res = rw.run_rewiring(prep, rw.RewiringConfig(n_permutations=50, n_bootstrap=20))
    rt = res.reactions
    row = rt[(rt.substrate == "Glutamine") & (rt["product"] == "Glutamate")].iloc[0]
    i, j = prep.features.index("Glutamine"), prep.features.index("Glutamate")
    ratio = prep.X[j] - prep.X[i]
    a, b = ratio[prep.idx_a], ratio[prep.idx_b]
    assert row.ratio_shift == pytest.approx(b.mean() - a.mean())
    assert row.ratio_p == pytest.approx(stats.ttest_ind(b, a, equal_var=False).pvalue)
    # planted: glutamine -0.8, glutamate +0.5 log2 in tumor -> ratio shift ~ +1.3
    assert 0.9 < row.ratio_shift < 1.7


def test_pathway_test_finds_planted_pathways_and_not_intact_ones(demo):
    _, prep = demo
    res = rw.run_rewiring(prep, rw.RewiringConfig(n_permutations=200, n_bootstrap=20))
    pt = res.pathway_tests.set_index("pathway")
    assert pt.loc["Tryptophan metabolism", "q_value"] < 0.05
    assert pt.loc["Alanine, aspartate and glutamate metabolism", "q_value"] < 0.05
    assert pt.loc["Citrate cycle (TCA cycle)", "q_value"] > 0.2     # TCA-internal couplings are intact
    assert (pt["n_measured"] >= 3).all()


def test_pathway_test_p_values_match_bruteforce(demo):
    _, prep = demo
    sets = {k: v for k, v in prep.sets.items() if len(v) >= 3}
    perm = rw.permutation_test(prep, "pearson", 60, seed=3, sets=sets)
    se = rw.pair_scale(prep.n_a, prep.n_b, "pearson")
    m = len(prep.features)
    iu, ju = np.triu_indices(m, k=1)
    arr, exact = rw.label_arrangements(prep, 60, seed=3)
    Xc = rw._group_centered(prep)
    name = "Tryptophan metabolism"
    idx = np.array(sorted(sets[name]))
    sel = np.isin(iu, idx) & np.isin(ju, idx)
    z0 = rw._pair_z(prep.X, prep.idx_a, prep.idx_b, "pearson", se, iu, ju)
    obs = np.mean(z0[sel] ** 2)
    allidx = np.arange(len(prep.group))
    ge = sum(np.mean(rw._pair_z(Xc, ia, np.setdiff1d(allidx, ia), "pearson", se, iu, ju)[sel] ** 2) >= obs - 1e-12
             for ia in arr)
    assert perm["_set_p"][name] == pytest.approx((1 + ge) / (1 + len(arr)))


def test_null_pathway_test_is_calibrated():
    ps = []
    for seed in range(12):
        d = rw.simulate_rewiring_demo(n_per_group=10, seed=seed, rewire=False, outliers=False, batch_effect=False)
        prep = rw.prepare_data(d["peak_df"], d["meta"], "Group", "Control", "Tumor", pathway_library="KEGG")
        perm = rw.permutation_test(prep, "pearson", 200, seed=seed,
                                   sets={k: v for k, v in prep.sets.items() if len(v) >= 3})
        ps += list(perm["_set_p"].values())
    ps = np.array(ps)
    assert 0.0 <= np.mean(ps < 0.05) <= 0.12
    assert 0.35 < np.median(ps) < 0.65


# ---------------------------------------------------------------------------
# small samples
# ---------------------------------------------------------------------------
def test_three_per_group_runs_exactly_and_honestly():
    d = rw.simulate_rewiring_demo(n_per_group=3, seed=11)
    prep = rw.prepare_data(d["peak_df"], d["meta"], "Group", "Control", "Tumor", pathway_library="KEGG")
    assert any("20 possible label splits" in n and "below 0.1" in n for n in prep.notes)
    res = rw.run_rewiring(prep, rw.RewiringConfig(n_permutations=500, n_bootstrap=50))
    perm = res.permutation
    assert perm["exact"] and perm["n_permutations"] == math.comb(6, 3) == 20 and perm["p_floor"] == 0.1
    assert res.pairs["p_value"].isna().all() and res.pairs["q_bh"].isna().all()   # analytic test undefined
    assert (res.pairs["q_perm"] >= 0.05 - 1e-12).all()
    assert (res.pathway_tests["p_value"] >= 0.1 - 1e-12).all()      # split + mirror image always tie
    assert res.permutation["p_value"] >= 0.1 - 1e-12
    assert all(h["verdict"] in ("exploratory", "fragile") for h in res.hypotheses)
    payload = res.to_payload()
    json.dumps(payload, allow_nan=False)
    assert len(payload["top_pairs"]) > 0


def test_three_per_group_requires_permutations():
    d = rw.simulate_rewiring_demo(n_per_group=3, seed=11)
    prep = rw.prepare_data(d["peak_df"], d["meta"], "Group", "Control", "Tumor")
    with pytest.raises(rw.RewiringInputError, match="permutations"):
        rw.run_rewiring(prep, rw.RewiringConfig(n_permutations=0))


def test_two_per_group_is_refused():
    d = rw.simulate_rewiring_demo(n_per_group=3, seed=1)
    meta = d["meta"].copy()
    meta.loc[meta.Sample == "TUM_03", "Group"] = "Other"
    with pytest.raises(rw.RewiringInputError, match="at least 3"):
        rw.prepare_data(d["peak_df"], meta, "Group", "Control", "Tumor")


def test_exact_global_p_matches_bruteforce_enumeration():
    d = rw.simulate_rewiring_demo(n_per_group=4, seed=2)
    prep = rw.prepare_data(d["peak_df"], d["meta"], "Group", "Control", "Tumor")
    perm = rw.permutation_test(prep, "pearson", 1000, seed=0)
    assert perm["exact"] and perm["n_permutations"] == 70
    m = len(prep.features)
    iu, ju = np.triu_indices(m, k=1)
    se = rw.pair_scale(4, 4, "pearson")
    obs = np.mean(rw._pair_z(prep.X, prep.idx_a, prep.idx_b, "pearson", se, iu, ju) ** 2)
    Xc = rw._group_centered(prep)
    stats_all = [np.mean(rw._pair_z(Xc, np.array(c), np.setdiff1d(np.arange(8), c), "pearson", se, iu, ju) ** 2)
                 for c in itertools.combinations(range(8), 4)]
    assert perm["p_value"] == pytest.approx(np.mean(np.array(stats_all) >= obs - 1e-12))


def test_small_group_hypothesis_checks_are_marked_not_assessable():
    d = rw.simulate_rewiring_demo(n_per_group=4, seed=4)
    prep = rw.prepare_data(d["peak_df"], d["meta"], "Group", "Control", "Tumor", pathway_library="KEGG")
    pairs = rw.annotate_reactions(prep, rw.differential_correlation(prep), 2)
    pairs["tested"] = True
    pairs["q_value"] = 0.001                      # force calls to exercise the card builder
    pairs = rw.classify_pairs(pairs, rw.RewiringConfig())
    cards = rw.build_hypotheses(prep, pairs, rw.feature_statistics(prep), rw.RewiringConfig(n_bootstrap=10))
    assert cards and all(c["stability"] is None and c["influence_check"]["ok"] is None for c in cards)
    assert {c["verdict"] for c in cards} <= {"exploratory", "fragile"}


def test_null_three_per_group_never_calls_pathways():
    """Regression: with |r| -> 1 in 3-sample groups, the observed statistic must be computed exactly
    like the null (it once came out larger than its own split, giving p = 0)."""
    for seed in range(8):
        d0 = rw.simulate_rewiring_demo(n_per_group=3, seed=5000 + seed, rewire=False, outliers=False,
                                       batch_effect=False)
        prep0 = rw.prepare_data(d0["peak_df"], d0["meta"], "Group", "Control", "Tumor", pathway_library="KEGG")
        res0 = rw.run_rewiring(prep0, rw.RewiringConfig(n_permutations=200, n_bootstrap=10, seed=seed))
        assert (res0.pathway_tests["p_value"] >= 0.1 - 1e-12).all()
        assert not (res0.pathway_tests["q_value"] < 0.05).any()


# ---------------------------------------------------------------------------
# cutoffs, network selection, figures
# ---------------------------------------------------------------------------
def test_p_value_cutoff_and_cutoff_of_one(demo):
    _, prep = demo
    res_q = rw.run_rewiring(prep, rw.RewiringConfig(n_permutations=200, n_bootstrap=10))
    res_p = rw.run_rewiring(prep, rw.RewiringConfig(n_permutations=200, n_bootstrap=10, cutoff_on="p", fdr=0.05))
    called = res_p.pairs["class"].isin(rw.REWIRED)
    assert (res_p.pairs.loc[called, "p_perm"] <= 0.05).all()
    assert called.sum() >= res_q.pairs["class"].isin(rw.REWIRED).sum()      # p <= 0.05 is looser than q <= 0.05
    pp = res_p.pairs["p_perm"].dropna()
    assert ((pp > 0) & (pp <= 1)).all()
    res_1 = rw.run_rewiring(prep, rw.RewiringConfig(n_permutations=50, n_bootstrap=10, fdr=1.0))
    t = res_1.pairs[res_1.pairs["tested"]]
    assert t["class"].isin(rw.REWIRED).sum() == (t["delta_r"].abs() >= res_1.cfg.min_delta).sum()
    with pytest.raises(rw.RewiringInputError):
        rw.RewiringConfig(fdr=1.5)


def test_network_groups_top_n_and_cutoff(demo):
    from metabolomics_modules import rewiring_figures as rf
    _, prep = demo
    res = rw.run_rewiring(prep, rw.RewiringConfig(n_permutations=200, n_bootstrap=10))
    g_all, order_all, _ = rf.network_groups(res, 0)
    assert all(x is not None for x in g_all)
    for top in (5, 10, 15, 20):
        g, order, sel = rf.network_groups(res, top, "q", 1.0)
        assert len(sel) == min(top, len(res.pathway_tests)) and len(order) <= len(sel)
        shown = {i for i, x in enumerate(g) if x is not None}
        # every shown metabolite belongs to its pathway; every metabolite of a selected pathway is shown
        members = {r.pathway: set(r.feature_idx) for r in sel.itertuples()}
        assert all(i in members[g[i]] for i in shown)
        assert set().union(*members.values()) == shown
    g, order, sel = rf.network_groups(res, 10, "p", 1e-9)
    assert sel.empty and all(x is None for x in g)


def test_figures_render_in_every_format(demo):
    from metabolomics_modules import rewiring_figures as rf
    _, prep = demo
    res = rw.run_rewiring(prep, rw.RewiringConfig(n_permutations=100, n_bootstrap=10))
    i, j = prep.features.index("Tryptophan"), prep.features.index("Kynurenine")
    magic = {"PNG": b"\x89PNG", "JPEG": b"\xff\xd8", "TIFF": (b"II*\x00", b"MM\x00*"), "SVG": b"<?xml",
             "PDF": b"%PDF", "EPS": b"%!PS"}
    for fmt, head in magic.items():
        heads = head if isinstance(head, tuple) else (head,)
        net = rf.figure_bytes(rf.network_figure(res), fmt, 100)
        pair = rf.figure_bytes(rf.pair_figure(res, i, j), fmt, 100)
        assert any(net.startswith(h) for h in heads) and any(pair.startswith(h) for h in heads), fmt
    svg = rf.figure_bytes(rf.pair_figure(res, i, j), "SVG").decode("utf-8")
    assert "Relative Abundance (Tryptophan)" in svg and "Relative Abundance (Kynurenine)" in svg


def _polys_overlap(A, B):
    for P in (A, B):
        for k in range(len(P)):
            (x1, y1), (x2, y2) = P[k], P[(k + 1) % len(P)]
            ax, ay = -(y2 - y1), x2 - x1
            pa = [x * ax + y * ay for x, y in A]
            pb = [x * ax + y * ay for x, y in B]
            if max(pa) <= min(pb) or max(pb) <= min(pa):
                return False
    return True


def _check_geometry(rf, geo):
    """Names never touch each other, a dot, the pathway arcs or the pathway labels; pathway labels sit
    outside the arcs and never touch each other; leader lines stay outside the arcs."""
    import math
    pos, rn, r_arc = geo["pos"], geo["rn"], geo["r_arc"]
    polys = {i: rf.name_polygon(geo, i) for i in pos} if geo["labels"] else {}
    ids = sorted(polys)
    for a in range(len(ids)):
        for b in range(a + 1, len(ids)):
            assert not _polys_overlap(polys[ids[a]], polys[ids[b]]), (ids[a], ids[b])
    for i, P in polys.items():
        assert max(math.hypot(x, y) for x, y in P) < r_arc - 2          # inside the arcs
        for k, (_, x, y) in pos.items():                                # clear of every dot
            assert min(math.hypot(px - x, py - y) for px, py in P) > rf.NODE_R, (i, k)
    boxes = [it["box"] for it in geo["items"]]
    for k, (x0, x1, y0, y1) in enumerate(boxes):
        dx = x0 if x0 > 0 else (-x1 if x1 < 0 else 0.0)
        dy = y0 if y0 > 0 else (-y1 if y1 < 0 else 0.0)
        assert math.hypot(dx, dy) > r_arc + 2, geo["items"][k]["pw"]       # outside the ring
        for b in boxes[k + 1:]:
            assert not (x0 < b[1] and x1 > b[0] and y0 < b[3] and y1 > b[2])
    for it in geo["items"]:
        if it["lead"]:
            (x1, y1), (x2, y2) = it["lead"]
            for t in range(21):
                assert math.hypot(x1 + (x2 - x1) * t / 20, y1 + (y2 - y1) * t / 20) > r_arc + 1
            assert not rf._seg_hits(it["lead"], [b for b in boxes if b is not it["box"]], 0)


def test_network_figure_names_never_overlap(demo):
    """Full names on a ring that grows with the number of metabolites; the same geometry is drawn."""
    from metabolomics_modules import rewiring_figures as rf
    _, prep = demo
    res = rw.run_rewiring(prep, rw.RewiringConfig(n_permutations=100, n_bootstrap=10))
    for top, metric, cut in ((10, "p", 0.05), (5, "q", 1.0), (20, "q", 1.0), (0, "q", 1.0)):
        geo = rf.network_geometry(res, top, metric, cut)
        assert geo["labels"] and len(geo["pos"]) > 0
        _check_geometry(rf, geo)
    # dense ring: 150 metabolites with long names in 15 pathways (some with long names)
    rng = np.random.default_rng(5)
    names = [f"{'Phosphatidylethanolamine-like ' if j % 3 == 0 else 'M'}{k}-{j}" for k in range(15) for j in range(10)]
    pw = [f"Pathway {k} of a long lipid and amino acid metabolism" if k % 2 else f"P{k}" for k in range(15) for _ in range(10)]
    X = rng.standard_normal((150, 24)) + 10
    df = pd.DataFrame(X, index=names, columns=[f"s{i}" for i in range(24)])
    meta = pd.DataFrame({"Group": ["A"] * 12 + ["B"] * 12}, index=df.columns)
    p2 = rw.prepare_data(df, meta, "Group", "A", "B", pathways=pd.Series(pw, index=names), log_transform="no",
                         max_features=1000)
    r2 = rw.run_rewiring(p2, rw.RewiringConfig(n_permutations=30, n_bootstrap=5))
    geo = rf.network_geometry(r2, 0)
    assert len(geo["pos"]) == 150 and geo["rn"] > 215
    _check_geometry(rf, geo)
    fig = rf.network_figure(r2, 0)
    fig.canvas.draw()
    import matplotlib.pyplot as plt
    plt.close(fig)


NAME_CASES = [  # (name as it appears in data files, accepted HMDB IDs)
    ("β-Alanine", {"HMDB0000056"}), ("beta-Alanine", {"HMDB0000056"}), ("L-Alanine", {"HMDB0000161"}),
    ("α-Ketoglutarate", {"HMDB0000208"}), ("γ-Aminobutyric acid", {"HMDB0000112"}), ("α-Tocopherol", {"HMDB0001893"}),
    ("β-Hydroxybutyrate", {"HMDB0000357", "HMDB0000011"}), ("N,N-Dimethylglycine", {"HMDB0000092"}),
    ("Glucose 6\u2011phosphate", {"HMDB0001401"}), ("NAD\u207a", {"HMDB0000902"}),
    ("Adenosine 5\u2032-monophosphate", {"HMDB0000045"}), ("3\u2032,5\u2032-Cyclic AMP", {"HMDB0000058"}),
    ("2\u2032-Deoxyguanosine", {"HMDB0000085"}), ("L\u2010Lactic acid", {"HMDB0000190"}),
    ("L\u2013Lactic acid", {"HMDB0000190"}), ("L-(+)-Lactic acid", {"HMDB0000190"}), ("D-(+)-Glucose", {"HMDB0000122"}),
    ("Uric\u00a0acid", {"HMDB0000289"}), ("Taurine ", {"HMDB0000251"}), ("Kynurenine*", {"HMDB0000684"}),
    ("Pantothenic acid (Vitamin B5)", {"HMDB0000210"}), ("Glutamine_[M+H]+", {"HMDB0000641"}),
    ("Tryptophan [M+H]+", {"HMDB0000929"}), ("Tyrosine [M-H]-", {"HMDB0000158"}), ("Serine_RT2.45", {"HMDB0000187"}),
    ("Leucine/Isoleucine", {"HMDB0000687", "HMDB0000172"}), ("Citrate; Isocitrate", {"HMDB0000094", "HMDB0000193"}),
]


def test_special_characters_in_names_are_matched_correctly():
    """Greek letters, Unicode dashes / primes / spaces / superscripts, instrument tags, footnote marks
    and isomer groups: the right metabolite, never a different one."""
    fm = rr.map_features([n for n, _ in NAME_CASES], None, "KEGG")
    got = dict(zip([n for n, _ in NAME_CASES], fm.table["HMDB_ID"]))
    wrong = {n: got[n] for n, want in NAME_CASES if got[n] not in want}
    assert not wrong, wrong
    assert list(fm.table["Metabolite"]) == [n for n, _ in NAME_CASES]          # display names untouched
    assert "Greek letters spelled out" in fm.table["Name_Cleaning"].iat[0]


def test_name_cleaning_does_not_create_false_matches():
    names = ["3-Hydroxybutyrate", "2-Hydroxybutyrate", "N-Acetylaspartate", "PC(16:0/18:1)", "TG 16:0/18:1/18:2",
             "Unknown_123 [M+H]+", "Alanine", "Ω-Unknownium"]
    fm = rr.map_features(names, None, "KEGG")
    h = dict(zip(names, fm.table["HMDB_ID"]))
    assert h["3-Hydroxybutyrate"] != h["2-Hydroxybutyrate"]                    # locants are kept
    assert h["Alanine"] == "HMDB0000161"
    assert h["Unknown_123 [M+H]+"] == "" and h["Ω-Unknownium"] == ""
    assert "isomer group" not in fm.table["Name_Cleaning"].iat[3]               # slash inside a lipid name
    assert "isomer group" not in fm.table["Name_Cleaning"].iat[4]
    t, alts, _ = rr.clean_query_name("β-Alanine")
    assert t == "beta-Alanine" and alts == ["beta-Alanine"]


def test_annotation_ids_win_over_name_cleaning():
    ann = pd.DataFrame({"HMDB": ["HMDB0000056"]}, index=["β-Ala (my label)"])
    fm = rr.map_features(["β-Ala (my label)"], ann, "KEGG")
    assert fm.table["HMDB_ID"].iat[0] == "HMDB0000056"


def test_co_products_of_one_reaction_are_not_a_route(lib):
    """Glutamate and glycine are both released when a glutathione conjugate is hydrolysed; that
    single reaction must not make them 'two steps apart'. A real two-reaction route stays."""
    g = _gem(["Glutamate", "Glycine", "Cysteine", "Glutathione"])
    L = lib.links(g, max_steps=2)
    assert (0, 1) not in L
    for (a, b), rec in L.items():
        if rec["steps"] == 2:
            first, second = rec["reactions"].split(" → ")
            assert not set(first.split(";")) & set(second.split(";"))
    cys_gsh = L.get((2, 3))
    assert cys_gsh is not None and cys_gsh["steps"] == 2 and cys_gsh["direction"] in ("i>j", "both")
