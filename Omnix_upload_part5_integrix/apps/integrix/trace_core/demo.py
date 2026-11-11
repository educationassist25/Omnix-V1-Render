"""
demo.py - synthetic paired tri-omics dataset with planted, known biology (ground truth for validation).

24 samples (12 Control, 12 Tumor), two processing batches. Planted effects (log2 units):
  1. Glycolysis up: enzymes up at mRNA and protein; glucose down, pyruvate and lactate up
     -> expected "Throughput-consistent (up)" reactions, Glycolysis shared-concordant.
  2. ASS1 loss (urea cycle): ASS1 down at mRNA and protein; citrulline up, argininosuccinate down
     -> expected "Candidate bottleneck" at L-citrulline:L-aspartate ligase.
  3. Fatty acid oxidation: proteins up without mRNA change -> "Protein-only (post-transcriptional)".
  4. Glutaminase step: GLS unchanged, glutamine up, glutamate down -> "Non-abundance regulation".
  5. Cholesterol biosynthesis: mRNA up, protein unchanged -> "mRNA-only (buffered)".
"""
from __future__ import annotations

import json
import os

import numpy as np
import pandas as pd

PLANT = {
    "glycolysis_genes": ["HK1", "HK2", "GPI", "PFKP", "PFKL", "ALDOA", "TPI1", "GAPDH", "PGK1", "PGAM1",
                         "ENO1", "PKM", "LDHA"],
    "glycolysis_mets": {"glucose": -0.7, "pyruvate": 0.8, "L-lactate": 1.1, "PEP": 0.6,
                        "fructose-1,6-bisphosphate": 0.6, "glucose-6-phosphate": 0.5},
    "bottleneck_gene": "ASS1",
    "bottleneck_mets": {"citrulline": 1.1, "argininosuccinate": -0.9},
    "fao_genes": ["CPT1A", "CPT2", "ACADM", "ACADVL", "HADHA", "HADHB", "ECHS1", "HADH", "ACAA2", "ACADS"],
    "nonabundance_gene": "GLS",
    "nonabundance_mets": {"glutamine": 0.9, "glutamate": -0.9},
    "buffered_genes": ["HMGCR", "HMGCS1", "MVK", "PMVK", "MVD", "FDPS", "FDFT1", "SQLE", "LSS", "CYP51A1",
                       "DHCR7", "IDI1"],
}


def make_demo(kb, out_dir: str, seed: int = 7, n_per_group: int = 12, effects: bool = True) -> dict:
    """effects=False writes a global-null dataset (no planted biology) for calibration tests."""
    rng = np.random.default_rng(seed)
    os.makedirs(out_dir, exist_ok=True)
    n = 2 * n_per_group
    samples = [f"S{i + 1:02d}" for i in range(n)]
    group = np.array(["Control"] * n_per_group + ["Tumor"] * n_per_group)
    batch = np.array(["A", "B"] * (n // 2))
    tumor = (group == "Tumor").astype(float) * (1.0 if effects else 0.0)
    meta = pd.DataFrame({"Sample": samples, "Group": group, "Batch": batch,
                         "Age": rng.integers(45, 80, n), "Sex": rng.choice(["F", "M"], n)})

    genes = kb.genes[kb.genes["symbol"].ne("")].drop_duplicates("symbol").reset_index(drop=True)
    G = len(genes)
    sym = genes["symbol"].values
    latent = rng.normal(size=(3, n))
    load = rng.normal(scale=0.25, size=(G, 3))
    mu = rng.normal(8.0, 1.6, G)
    rna = mu[:, None] + load @ latent + rng.normal(scale=0.35, size=(G, n))
    eff_r = np.zeros(G)
    eff_p = np.zeros(G)
    gi = {s: i for i, s in enumerate(sym)}
    for s in PLANT["glycolysis_genes"]:
        if s in gi:
            eff_r[gi[s]] += 0.9
            eff_p[gi[s]] += 0.9
    if PLANT["bottleneck_gene"] in gi:
        eff_r[gi[PLANT["bottleneck_gene"]]] -= 1.3
        eff_p[gi[PLANT["bottleneck_gene"]]] -= 1.3
    for s in PLANT["fao_genes"]:
        if s in gi:
            eff_p[gi[s]] += 0.9
    for s in PLANT["buffered_genes"]:
        if s in gi:
            eff_r[gi[s]] += 0.9
            eff_p[gi[s]] -= 0.0
    rna_eff = rna + eff_r[:, None] * tumor[None, :]
    # protein: coupled to the mRNA deviation, with its own noise and planted effects
    mu_p = rng.normal(20.0, 2.0, G)
    prot = (mu_p[:, None] + 0.55 * (rna - rna.mean(1, keepdims=True)) + rng.normal(scale=0.3, size=(G, n))
            + eff_p[:, None] * tumor[None, :])
    # buffered genes: protein does not follow their mRNA effect (only baseline coupling)
    # proteome coverage: ~65% of metabolic genes + every planted gene
    planted = set(PLANT["glycolysis_genes"]) | set(PLANT["fao_genes"]) | set(PLANT["buffered_genes"]) | \
        {PLANT["bottleneck_gene"], PLANT["nonabundance_gene"]}
    keep_p = rng.random(G) < 0.65
    keep_p |= np.isin(sym, list(planted))

    # metabolites: non-currency, with reactions, prefer those with HMDB IDs
    m = kb.metabolites
    cur = set(m.loc[m["currency"], "met_base"])
    rx_ok = kb.reactions.loc[~(kb.reactions["is_transport"] | kb.reactions["is_exchange"]), "rxn"]
    deg = kb.stoich[kb.stoich["rxn"].isin(rx_ok) & ~kb.stoich["met_base"].isin(cur)].groupby("met_base").size()
    cand = m[m["met_base"].isin(deg.index) & ~m["currency"]].copy()
    cand["deg"] = cand["met_base"].map(deg)
    cand = cand[(cand["hmdb"] != "") | (cand["kegg"] != "")]
    named = set(PLANT["glycolysis_mets"]) | set(PLANT["bottleneck_mets"]) | set(PLANT["nonabundance_mets"])
    chosen = pd.concat([cand[cand["name"].isin(named)], cand[~cand["name"].isin(named)]
                        .sort_values("deg", ascending=False).head(320)]).drop_duplicates("met_base")
    M = len(chosen)
    mu_m = rng.normal(16.0, 2.0, M)
    met = mu_m[:, None] + (rng.normal(scale=0.2, size=(M, 3)) @ latent) + rng.normal(scale=0.35, size=(M, n))
    eff_m = np.zeros(M)
    nm = chosen["name"].values
    for d in (PLANT["glycolysis_mets"], PLANT["bottleneck_mets"], PLANT["nonabundance_mets"]):
        for k, v in d.items():
            eff_m[nm == k] += v
    met += eff_m[:, None] * tumor[None, :]
    # batch effect in the metabolomics layer (location shift per feature)
    met += rng.normal(0.35, 0.15, M)[:, None] * (batch == "B")[None, :]

    def to_raw(X, miss=0.03):
        R = 2 ** X
        low = X < np.quantile(X, 0.08)
        R[low & (rng.random(X.shape) < 0.3)] = np.nan
        R[rng.random(X.shape) < miss * 0.2] = np.nan
        return R

    t_rna = pd.DataFrame(to_raw(rna_eff), columns=samples)
    t_rna.insert(0, "Gene", sym)
    pr = genes[keep_p]
    t_prot = pd.DataFrame(to_raw(prot[keep_p]), columns=samples)
    t_prot.insert(0, "Gene", pr["symbol"].values)
    t_prot.insert(0, "Protein", [u.split(";")[0] if u else s for u, s in zip(pr["uniprot"], pr["symbol"])])
    # metabolites: identify ~70% by HMDB, the rest by name only (tests name matching)
    by_id = rng.random(M) < 0.7
    t_met = pd.DataFrame(to_raw(met), columns=samples)
    t_met.insert(0, "HMDB", np.where(by_id, chosen["hmdb"].values, ""))
    t_met.insert(0, "Metabolite", nm)
    t_rna.to_csv(os.path.join(out_dir, "demo_transcriptomics.csv"), index=False)
    t_prot.to_csv(os.path.join(out_dir, "demo_proteomics.csv"), index=False)
    t_met.to_csv(os.path.join(out_dir, "demo_metabolomics.csv"), index=False)
    meta.to_csv(os.path.join(out_dir, "demo_metadata.csv"), index=False)
    truth = {k: (v if not isinstance(v, dict) else v) for k, v in PLANT.items()}
    truth["expected"] = {
        "throughput_up_pathway": "Glycolysis / Gluconeogenesis",
        "bottleneck_reaction_name": "L-Citrulline:L-aspartate ligase (AMP-forming)",
        "nonabundance_reaction_name": "L-glutamine amidohydrolase",
    }
    json.dump(truth, open(os.path.join(out_dir, "demo_ground_truth.json"), "w"), indent=2)
    return {"rna": t_rna.shape, "prot": t_prot.shape, "met": t_met.shape}


if __name__ == "__main__":
    from .kb import load_kb
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    print(make_demo(load_kb(), os.path.join(here, "demo_data")))


# ----------------------------------------------------------------------------- factorial demo (v1.7)
FACTORIAL = {
    "cancer": {
        "genes": {  # gene: (mRNA effect, protein effect) in log2 units, Cancer over Normal
            **{g: (0.9, 0.9) for g in PLANT["glycolysis_genes"]},
            "ASS1": (-1.3, -1.3),
            "IDO1": (1.2, 1.2),
            **{g: (0.0, 0.9) for g in PLANT["fao_genes"]},
            **{g: (0.9, 0.0) for g in PLANT["buffered_genes"]},
        },
        "mets": {**PLANT["glycolysis_mets"], **PLANT["bottleneck_mets"], "tryptophan": -0.9,
                 "L-formylkynurenine": 1.0, "kynurenine": 1.1},
    },
    "treatment": {  # glutaminase inhibitor: enzyme amount unchanged, its activity blocked
        "genes": {},
        "mets": {"glutamine": 1.3, "glutamate": -1.3},
    },
    "cancer_x_treatment": {  # the drug partly reverses the Warburg shift, in cancer samples only
        "genes": {g: (-0.6, -0.6) for g in PLANT["glycolysis_genes"]},
        "mets": {k: -0.5 * v for k, v in PLANT["glycolysis_mets"].items()},
    },
    "female": {  # steroid-hormone enzymes; aromatase product up, substrate down
        "genes": {"CYP19A1": (1.2, 1.2), "HSD17B1": (0.8, 0.8), "HSD3B1": (0.7, 0.7)},
        "mets": {"testosterone": -1.0, "19-hydroxytestosterone": 1.0, "estradiol-17beta": 1.0, "estrone": 0.8},
    },
}


def make_factorial_demo(kb, out_dir: str, seed: int = 13, n_per_cell: int = 5) -> dict:
    """Synthetic 40-sample study: Condition (Normal, Cancer) x Treatment (Control, Treatment) x Sex (Male,
    Female), 5 samples per combination, three batches, realistic missing values (8% at random in proteomics,
    below-detection censoring in every layer) and a batch shift in metabolomics batch C. Planted biology per
    factor is in FACTORIAL and written to ground_truth.json."""
    rng = np.random.default_rng(seed)
    os.makedirs(out_dir, exist_ok=True)
    rows = [(c, t, s) for c in ("Normal", "Cancer") for t in ("Control", "Treatment") for s in ("Male", "Female")
            for _ in range(n_per_cell)]
    n = len(rows)
    samples = [f"S{i + 1:02d}" for i in range(n)]
    cond = np.array([r[0] for r in rows])
    trt = np.array([r[1] for r in rows])
    sex = np.array([r[2] for r in rows])
    batch = np.array([["A", "B", "C"][i % 3] for i in rng.permutation(n)])
    meta = pd.DataFrame({"Sample": samples, "Group": [f"{c} {t}" for c, t in zip(cond, trt)], "Condition": cond,
                         "Treatment": trt, "Sex": sex, "Batch": batch, "Age": rng.integers(40, 82, n)})
    x = {"cancer": (cond == "Cancer").astype(float), "treatment": (trt == "Treatment").astype(float),
         "female": (sex == "Female").astype(float)}
    x["cancer_x_treatment"] = x["cancer"] * x["treatment"]

    genes = kb.genes[kb.genes["symbol"].ne("")].drop_duplicates("symbol").reset_index(drop=True)
    G = len(genes)
    sym = genes["symbol"].values
    gi = {s: i for i, s in enumerate(sym)}
    latent = rng.normal(size=(3, n))
    rna = rng.normal(8.0, 1.6, G)[:, None] + rng.normal(scale=0.25, size=(G, 3)) @ latent + \
        rng.normal(scale=0.35, size=(G, n))
    prot = rng.normal(20.0, 2.0, G)[:, None] + 0.55 * (rna - rna.mean(1, keepdims=True)) + \
        rng.normal(scale=0.3, size=(G, n))
    for f, spec in FACTORIAL.items():
        for g, (er, ep) in spec["genes"].items():
            if g in gi:
                rna[gi[g]] += er * x[f]
                prot[gi[g]] += ep * x[f]
    planted = {g for spec in FACTORIAL.values() for g in spec["genes"]} | {"GLS"}
    keep_p = (rng.random(G) < 0.65) | np.isin(sym, list(planted))

    m = kb.metabolites
    cur = set(m.loc[m["currency"], "met_base"])
    rx_ok = kb.reactions.loc[~(kb.reactions["is_transport"] | kb.reactions["is_exchange"]), "rxn"]
    deg = kb.stoich[kb.stoich["rxn"].isin(rx_ok) & ~kb.stoich["met_base"].isin(cur)].groupby("met_base").size()
    cand = m[m["met_base"].isin(deg.index) & ~m["currency"]].copy()
    cand["deg"] = cand["met_base"].map(deg)
    cand = cand[(cand["hmdb"] != "") | (cand["kegg"] != "")]
    named = {k for spec in FACTORIAL.values() for k in spec["mets"]}
    chosen = pd.concat([cand[cand["name"].isin(named)], cand[~cand["name"].isin(named)]
                        .sort_values("deg", ascending=False).head(330)]).drop_duplicates("met_base")
    M = len(chosen)
    nm = chosen["name"].values
    met = rng.normal(16.0, 2.0, M)[:, None] + rng.normal(scale=0.2, size=(M, 3)) @ latent + \
        rng.normal(scale=0.35, size=(M, n))
    for f, spec in FACTORIAL.items():
        for k, v in spec["mets"].items():
            met[nm == k] += v * x[f]
    met += rng.normal(0.4, 0.15, M)[:, None] * (batch == "C")[None, :]     # batch shift, not tied to a group

    def to_raw(X, mar=0.0):
        R = 2 ** X
        low = X < np.quantile(X, 0.08)
        R[low & (rng.random(X.shape) < 0.3)] = np.nan                         # below the detection limit
        R[rng.random(X.shape) < mar] = np.nan                                  # missing at random
        return R

    t_rna = pd.DataFrame(to_raw(rna, 0.005), columns=samples)
    t_rna.insert(0, "Gene", sym)
    pr = genes[keep_p]
    t_prot = pd.DataFrame(to_raw(prot[keep_p], 0.08), columns=samples)
    t_prot.insert(0, "Gene", pr["symbol"].values)
    t_prot.insert(0, "Protein", [u.split(";")[0] if u else s for u, s in zip(pr["uniprot"], pr["symbol"])])
    by_id = rng.random(M) < 0.7
    t_met = pd.DataFrame(to_raw(met, 0.01), columns=samples)
    t_met.insert(0, "HMDB", np.where(by_id, chosen["hmdb"].values, ""))
    t_met.insert(0, "Metabolite", nm)
    t_rna.to_csv(os.path.join(out_dir, "transcriptomics.csv"), index=False)
    t_prot.to_csv(os.path.join(out_dir, "proteomics.csv"), index=False)
    t_met.to_csv(os.path.join(out_dir, "metabolomics.csv"), index=False)
    meta.to_csv(os.path.join(out_dir, "metadata.csv"), index=False)
    truth = {f: {"genes": {g: list(v) for g, v in s["genes"].items()}, "mets": s["mets"]} for f, s in FACTORIAL.items()}
    truth["expected"] = {
        "Condition: Cancer over Normal": ["Glycolysis throughput up", "ASS1 candidate bottleneck",
                                        "IDO1 tryptophan -> kynurenine", "FAO protein-only", "cholesterol mRNA-only"],
        "Treatment: Treatment over Control": ["glutaminase non-abundance regulation (glutamine up, glutamate down)"],
        "Group: Cancer Treatment over Cancer Control": ["glutaminase step", "glycolysis partly reversed"],
        "Sex: Female over Male": ["aromatase CYP19A1 (testosterone down, estrone up)", "otherwise quiet"],
    }
    json.dump(truth, open(os.path.join(out_dir, "ground_truth.json"), "w"), indent=2)
    return {"rna": t_rna.shape, "prot": t_prot.shape, "met": t_met.shape, "samples": n}


def make_normalized_demo(src_dir: str, out_dir: str) -> dict:
    """The factorial demo as the Omnix modules would export it: log2 scale, normalized per sample, missing values
    left as gaps (no imputation), no batch correction, plus one Omnix sidecar per layer (omnix_io contract)."""
    from .omnix_io import sidecar
    os.makedirs(out_dir, exist_ok=True)
    spec = {"rna": ("transcriptomics", "BulkRNAAI Pro", ["Gene"], "log2, per-sample median scaling"),
            "prot": ("proteomics", "ProteoAI Pro", ["Protein", "Gene"], "log2, per-sample median centring"),
            "met": ("metabolomics", "MetaboAI Pro", ["Metabolite", "HMDB"],
                    "PQN (probabilistic quotient normalization), log2")}
    shapes = {}
    for k, (name, module, idc, how) in spec.items():
        df = pd.read_csv(os.path.join(src_dir, f"{name}.csv"))
        ids = df[idc]
        X = df.drop(columns=idc).astype(float)
        if k == "met":                                   # PQN on the raw scale, then log2
            ref = X.median(axis=1)
            q = X.div(ref, axis=0).median(axis=0)
            L = np.log2(X.div(q, axis=1))
        else:                                            # log2, then align every sample's median to the overall one
            L = np.log2(X)
            L = L - L.median(axis=0) + float(np.nanmedian(L.median(axis=0)))
        out = pd.concat([ids, L.round(5)], axis=1)
        out.to_csv(os.path.join(out_dir, f"{name}.csv"), index=False)
        info = sidecar(name, module, how, idc, f"{name}.csv", imputation="none", batch_corrected=False,
                       filtering="none (all features kept)", module_version="example",
                       created="2026-11-05T09:00:00")
        json.dump(info, open(os.path.join(out_dir, f"{name}_omnix.json"), "w"), indent=2)
        shapes[k] = out.shape
    import shutil
    shutil.copy(os.path.join(src_dir, "metadata.csv"), os.path.join(out_dir, "metadata.csv"))
    shutil.copy(os.path.join(src_dir, "ground_truth.json"), os.path.join(out_dir, "ground_truth.json"))
    return shapes
