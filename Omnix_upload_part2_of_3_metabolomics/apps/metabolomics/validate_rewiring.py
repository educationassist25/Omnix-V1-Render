"""
validate_rewiring.py - Simulation study for the Correlation Rewiring Map.

    python validate_rewiring.py            # ~1-2 minutes
    python validate_rewiring.py --quick    # fewer repetitions

Reports, for the engine's default settings (Pearson, permutation FDR < 0.05, |Δr| >= 0.45):
  1. Recovery of planted rewiring (recall per class, class accuracy, empirical false
     discovery proportion) over independent simulated datasets, n = 30 per group.
  2. Calibration under the global null (no rewiring), n = 30 per group.
  3. Small-group calibration (n = 9 per group): permutation FDR vs Benjamini-Hochberg,
     on Gaussian data and on the Rich Clinical demo with shuffled labels.
  4. Hypothesis-card checks: the outlier decoy is flagged; true modules are not.
  5. Sample-size study (3, 4, 6, 10, 20 per group): pair-level vs reaction-scoped vs
     pathway-level (MSEA-style) power, and calibration under the null.
"""

import argparse
import os
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from metabolomics_modules import normalization, rewiring as rw, utils  # noqa: E402


def demo_prep(seed, **kw):
    d = rw.simulate_rewiring_demo(seed=seed, **kw)
    prep = rw.prepare_data(d["peak_df"], d["meta"], "Group", "Control", "Tumor",
                           pathways=d["annotations"].set_index("Metabolite")["Pathway"], batch_col="Batch")
    return d, prep


def found_pairs(res):
    return {frozenset((a, b)): c for a, b, c in zip(res.pairs.feature_a, res.pairs.feature_b, res.pairs["class"])
            if c in rw.REWIRED}


def recovery(n_seeds, n_perm):
    rec = {k: [] for k in ("lost", "gained", "flipped")}
    acc, fdp, decoy, true_flagged = [], [], [], 0
    for seed in range(n_seeds):
        d, prep = demo_prep(seed)
        res = rw.run_rewiring(prep, rw.RewiringConfig(n_permutations=n_perm, n_bootstrap=100, seed=seed))
        t, f = d["truth"], found_pairs(res)
        true_any = t["lost"] | t["gained"] | t["flipped"]
        for k in rec:
            rec[k].append(np.mean([p in f for p in t[k]]))
        hits = [(p, k) for k in rec for p in t[k] if p in f]
        acc.append(np.mean([f[p] == k for p, k in hits]))
        disc = [p for p in f if p not in t["ambiguous"] and p not in t["decoy"]]
        fdp.append(np.mean([p not in true_any for p in disc]) if disc else 0.0)
        for h in res.hypotheses:
            pairs = [frozenset(x) for x in h["pair_names"]]
            if any(p in t["decoy"] for p in pairs):
                decoy.append(h["verdict"])
            elif pairs and all(p in true_any for p in pairs) and h["verdict"] != "robust":
                true_flagged += 1
    return rec, acc, fdp, decoy, true_flagged


def null_n30(n_reps, n_perm):
    any_disc, pvals = [], []
    for seed in range(1000, 1000 + n_reps):
        _, prep = demo_prep(seed, rewire=False, outliers=False)
        res = rw.run_rewiring(prep, rw.RewiringConfig(n_permutations=n_perm, n_bootstrap=10, seed=seed))
        any_disc.append(sum(res.counts()[k] for k in rw.REWIRED) > 0)
        pvals.append(res.permutation["p_value"])
    return np.array(any_disc), np.array(pvals)


def small_n(n_reps, n_perm):
    rng = np.random.default_rng(0)
    out = {"gaussian": {"bh": [], "permutation": []}, "richdemo": {"bh": [], "permutation": []}}
    for rep in range(n_reps):
        L = rng.standard_normal((147, 5)) * 0.7
        X = L @ rng.standard_normal((5, 18)) + rng.standard_normal((147, 18)) * 0.6 + 15
        df = pd.DataFrame(X, index=[f"f{i}" for i in range(147)], columns=[f"s{j}" for j in range(18)])
        meta = pd.DataFrame({"Group": ["A"] * 9 + ["B"] * 9}, index=df.columns)
        prep = rw.prepare_data(df, meta, "Group", "A", "B", log_transform="no")
        for fm in ("bh", "permutation"):
            r = rw.run_rewiring(prep, rw.RewiringConfig(fdr_method=fm, n_permutations=n_perm, n_bootstrap=5, seed=rep))
            out["gaussian"][fm].append(sum(r.counts()[k] for k in rw.REWIRED))
    peak = utils.validate_peak_matrix(pd.read_csv(os.path.join(HERE, "sample_peak_area_matrix_richdemo.csv")))
    peak = normalization.drop_istds(peak, normalization.detect_istds(peak.index))
    meta = utils.validate_metadata(pd.read_csv(os.path.join(HERE, "sample_metadata_richdemo.csv")), peak.columns)
    bio = meta[meta.Group.isin(["Healthy Control", "Prediabetic"])].copy()
    r2 = np.random.default_rng(5)
    for rep in range(n_reps):
        b = bio.copy()
        b["Group"] = r2.permutation(b["Group"].values)
        prep = rw.prepare_data(peak[b.index], b, "Group", "Healthy Control", "Prediabetic", sample_median_center=True)
        for fm in ("bh", "permutation"):
            r = rw.run_rewiring(prep, rw.RewiringConfig(fdr_method=fm, n_permutations=n_perm, n_bootstrap=5, seed=rep))
            out["richdemo"][fm].append(sum(r.counts()[k] for k in rw.REWIRED))
    return out


PLANTED_PATHWAYS = ("Tryptophan metabolism", "Alanine, aspartate and glutamate metabolism")


def sample_size_study(n_reps, n_perm, sizes=(3, 4, 6, 10, 20)):
    """Planted and null datasets at several group sizes, KEGG (MSEA) library, no outliers.
    For each size: pair-level recall over all planted pairs and over the planted pairs that are
    linked by a reaction (testing all pairs vs. reaction-linked pairs only), any-call rate and
    false discovery proportion, pathway-level detection of the two planted pathways, and under the
    null the share of datasets with any pair / any pathway called and the share of pathway p < 0.05."""
    out = {}
    for n in sizes:
        r = {"rec_all": [], "rec_rx_all": [], "rec_rx_scope": [], "fdp": [], "pw_hit": {k: [] for k in PLANTED_PATHWAYS},
             "null_any_pair": [], "null_any_pw": [], "null_set_p": [], "null_any_pair_rx": []}
        for seed in range(n_reps):
            d = rw.simulate_rewiring_demo(n_per_group=n, seed=1000 + seed, outliers=False)
            prep = rw.prepare_data(d["peak_df"], d["meta"], "Group", "Control", "Tumor", pathway_library="KEGG",
                                   batch_col="Batch" if n >= 10 else None)
            t = d["truth"]
            true_any = t["lost"] | t["gained"] | t["flipped"]
            res = rw.run_rewiring(prep, rw.RewiringConfig(n_permutations=n_perm, n_bootstrap=10, seed=seed))
            f = found_pairs(res)
            linked = {frozenset((a, b)) for a, b, s in zip(res.pairs.feature_a, res.pairs.feature_b, res.pairs.rxn_steps) if s}
            true_rx = true_any & linked
            r["rec_all"].append(np.mean([p in f for p in true_any]))
            r["rec_rx_all"].append(np.mean([p in f for p in true_rx]) if true_rx else np.nan)
            disc = [p for p in f if p not in t["ambiguous"]]
            r["fdp"].append(np.mean([p not in true_any for p in disc]) if disc else 0.0)
            pt = res.pathway_tests.set_index("pathway") if len(res.pathway_tests) else pd.DataFrame()
            for k in PLANTED_PATHWAYS:
                r["pw_hit"][k].append(bool(k in pt.index and pt.loc[k, "q_value"] < 0.05))
            res_rx = rw.run_rewiring(prep, rw.RewiringConfig(n_permutations=n_perm, n_bootstrap=10, seed=seed,
                                                             pair_scope="reaction"))
            f_rx = found_pairs(res_rx)
            r["rec_rx_scope"].append(np.mean([p in f_rx for p in true_rx]) if true_rx else np.nan)
            # null
            d0 = rw.simulate_rewiring_demo(n_per_group=n, seed=5000 + seed, rewire=False, outliers=False,
                                           batch_effect=False)
            prep0 = rw.prepare_data(d0["peak_df"], d0["meta"], "Group", "Control", "Tumor", pathway_library="KEGG")
            res0 = rw.run_rewiring(prep0, rw.RewiringConfig(n_permutations=n_perm, n_bootstrap=10, seed=seed))
            r["null_any_pair"].append(sum(res0.counts()[c] for c in rw.REWIRED) > 0)
            pt0 = res0.pathway_tests
            r["null_any_pw"].append(bool(len(pt0) and (pt0["q_value"] < 0.05).any()))
            r["null_set_p"] += list(pt0["p_value"]) if len(pt0) else []
            res0rx = rw.run_rewiring(prep0, rw.RewiringConfig(n_permutations=n_perm, n_bootstrap=10, seed=seed,
                                                              pair_scope="reaction"))
            r["null_any_pair_rx"].append(sum(res0rx.counts()[c] for c in rw.REWIRED) > 0)
        out[n] = r
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args()
    n_seeds, n_null, n_small, n_perm = (8, 10, 10, 100) if a.quick else (30, 30, 30, 200)
    t0 = time.time()
    rec, acc, fdp, decoy, true_flagged = recovery(n_seeds, n_perm)
    anyd, pv = null_n30(n_null, n_perm)
    sn = small_n(n_small, 100)
    ss = sample_size_study(8 if a.quick else 20, 200)
    print(f"\n## Rewiring Map validation ({n_seeds} planted datasets, {n_null} null datasets, "
          f"{n_small} small-n datasets per scenario; {time.time() - t0:.0f}s)\n")
    print("### 1. Recovery of planted rewiring (n = 30 per group)\n")
    print("| Metric | Mean | Min |\n|---|---|---|")
    for k in rec:
        print(f"| Recall, {k} pairs | {np.mean(rec[k]):.2f} | {np.min(rec[k]):.2f} |")
    print(f"| Class accuracy of recovered pairs | {np.mean(acc):.3f} | {np.min(acc):.3f} |")
    print(f"| False discovery proportion (target ≤ 0.05) | {np.mean(fdp):.3f} | max {np.max(fdp):.3f} |")
    print("\n### 2. Global null, n = 30 per group\n")
    print(f"- Datasets with any rewired pair called: {anyd.mean():.2f} (target ≤ 0.05)")
    print(f"- Global permutation p < 0.05: {np.mean(pv < 0.05):.2f} (target ≈ 0.05); median p = {np.median(pv):.2f}")
    print("\n### 3. Small groups (n = 9 per group), share of null datasets with any call\n")
    print("| Scenario | Benjamini–Hochberg | Permutation FDR |\n|---|---|---|")
    for sc, lab in (("gaussian", "Gaussian, 147 correlated features"), ("richdemo", "Rich Clinical demo, labels shuffled")):
        bh = np.mean(np.array(sn[sc]["bh"]) > 0)
        pm = np.mean(np.array(sn[sc]["permutation"]) > 0)
        print(f"| {lab} | {bh:.2f} (mean {np.mean(sn[sc]['bh']):.1f} pairs) | {pm:.2f} (mean {np.mean(sn[sc]['permutation']):.1f} pairs) |")
    print("\n### 4. Hypothesis-card checks\n")
    called = len(decoy)
    flagged = sum(v != "robust" for v in decoy)
    print(f"- Outlier decoy called as a pair in {called}/{n_seeds} datasets; card flagged (mixed/fragile) in {flagged}/{called}")
    print(f"- Modules made only of true rewired pairs that were flagged non-robust: {true_flagged}")
    n_ss = len(next(iter(ss.values()))["rec_all"])
    print(f"\n### 5. Sample size: pair-level vs reaction-scoped vs pathway-level ({n_ss} planted + {n_ss} null datasets per size)\n")
    print("Planted rewiring as in section 1 (no outlier decoy), KEGG library from the MSEA page. Reaction-linked recall "
          "is over the planted pairs that are linked by a Human-GEM reaction (direct or 2 steps).\n")
    print("| n per group | Label splits | Recall, all pairs | Recall of reaction-linked pairs: all-pairs test | "
          "…: reaction-scoped test | FDP | Tryptophan metabolism detected | Ala/Asp/Glu metabolism detected |")
    print("|---|---|---|---|---|---|---|---|")
    import math as _m
    for n, r in ss.items():
        splits = _m.comb(2 * n, n)
        print(f"| {n} | {splits if splits <= 200 else '> 200 (random)'} | {np.mean(r['rec_all']):.2f} | "
              f"{np.nanmean(r['rec_rx_all']):.2f} | {np.nanmean(r['rec_rx_scope']):.2f} | {np.mean(r['fdp']):.3f} | "
              f"{np.mean(r['pw_hit'][PLANTED_PATHWAYS[0]]):.2f} | {np.mean(r['pw_hit'][PLANTED_PATHWAYS[1]]):.2f} |")
    print("\nUnder the null (no rewiring):\n")
    print("| n per group | Datasets with any pair called (all / reaction-scoped) | Datasets with any pathway q < 0.05 | "
          "Pathway p < 0.05 (target ≈ 0.05) |")
    print("|---|---|---|---|")
    for n, r in ss.items():
        sp = np.array(r["null_set_p"])
        print(f"| {n} | {np.mean(r['null_any_pair']):.2f} / {np.mean(r['null_any_pair_rx']):.2f} | "
              f"{np.mean(r['null_any_pw']):.2f} | {np.mean(sp < 0.05):.3f} |")


if __name__ == "__main__":
    main()
