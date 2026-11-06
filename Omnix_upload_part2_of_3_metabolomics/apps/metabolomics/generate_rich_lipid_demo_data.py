"""
generate_rich_lipid_demo_data.py - the Lipidomics version of the Rich Clinical Demo.

Produces two files (the sample metadata is shared with the metabolite Rich Clinical Demo):
  1. sample_peak_area_matrix_richdemo_lipidomics.csv - 127 lipids + 1 internal standard x 42 samples
  2. lipid_row_annotations_richdemo.csv              - Lipid_Class / Category / HMDB_ID / ISTD per lipid
  (sample_metadata_richdemo.csv - 36 biological samples in 4 diagnosis groups + 6 QC replicates, reused)

Design (mirrors generate_rich_demo_data.py):
  - Real lipid names and HMDB IDs from the lipids in metabolite_row_annotations_multimethod.csv.
  - Lognormal per-lipid baseline, small technical noise, a shared per-class "activity" factor per biological
    sample (lipids of one class co-vary), and diagnosis effects that follow known dysmetabolic lipid patterns:
    triacylglycerols, diacylglycerols, ceramides and acylcarnitines up, lysophosphatidylcholines and
    PUFA-containing phospholipids down in Prediabetic < Metabolic Syndrome / Type 2 Diabetes.
  - A per-sample injection/extraction factor applied to every row, including the internal standard
    ISTD_PC(15:0/18:1-d7), so targeted ISTD normalization has a real effect.
  - Tiered missingness (exact zeros) in biological samples, as in the other demos.
"""
import os
import re
import numpy as np
import pandas as pd

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
SEED = 1777
ISTD = "ISTD_PC(15:0/18:1-d7)"


def lipid_class(name: str) -> str:
    m = re.match(r"^(LysoPC|LysoPE|GlcCer|LacCer|Cer|PC|PE|PS|PA|PG|PI|SM|TG|DG|CE)\(", name)
    if m:
        return m.group(1)
    low = name.lower()
    if "carnitine" in low:
        return "Acylcarnitine"
    if "coa" in low:
        return "Acyl-CoA"
    if any(x in low for x in ("cholic", "cholate", "bile")):
        return "Bile acid"
    if low.endswith(" acid"):
        return "Fatty acid"
    if any(x in low for x in ("sterol", "sterone", "estradiol", "cortico", "pregnenolone", "progesterone",
                              "testosterone", "androst", "cortisol", "cortisone", "aldosterone")):
        return "Sterol / steroid"
    return "Other lipid"


CATEGORY = {
    "PC": "Glycerophospholipids", "PE": "Glycerophospholipids", "PS": "Glycerophospholipids",
    "PA": "Glycerophospholipids", "PG": "Glycerophospholipids", "PI": "Glycerophospholipids",
    "LysoPC": "Glycerophospholipids", "LysoPE": "Glycerophospholipids",
    "TG": "Glycerolipids", "DG": "Glycerolipids",
    "Cer": "Sphingolipids", "SM": "Sphingolipids", "GlcCer": "Sphingolipids", "LacCer": "Sphingolipids",
    "CE": "Sterol lipids", "Sterol / steroid": "Sterol lipids", "Bile acid": "Sterol lipids",
    "Fatty acid": "Fatty acyls", "Acylcarnitine": "Fatty acyls", "Acyl-CoA": "Fatty acyls",
    "Other lipid": "Other",
}

# Fold change (vs Healthy Control) per lipid class: (Prediabetic, Type 2 Diabetes, Metabolic Syndrome)
CLASS_EFFECT = {
    "TG": (1.35, 2.2, 2.0), "DG": (1.25, 1.8, 1.7), "Cer": (1.25, 1.7, 1.6), "Acylcarnitine": (1.2, 1.7, 1.5),
    "LysoPC": (0.8, 0.55, 0.6), "CE": (1.0, 0.85, 0.9),
}
PUFA = re.compile(r"20:5|22:6|20:4|22:5")


def generate():
    rng = np.random.default_rng(SEED)
    meta = pd.read_csv(os.path.join(BASE_DIR, "sample_metadata_richdemo.csv"))
    ann = pd.read_csv(os.path.join(BASE_DIR, "metabolite_row_annotations_multimethod.csv"))
    lip = ann[ann["Class"] == "Lipid"].drop_duplicates("Metabolite")
    names = lip["Metabolite"].tolist()
    hmdb = lip["HMDB_ID"].fillna("").tolist()
    classes = [lipid_class(n) for n in names]
    samples = meta["Sample"].tolist()
    is_qc = meta["IsQC"].astype(str).str.lower().eq("true").tolist()
    dx = meta["Diagnosis"].tolist()
    n, m = len(names), len(samples)
    bio = [i for i, q in enumerate(is_qc) if not q]

    data = rng.lognormal(mean=10.5, sigma=1.1, size=(n, 1)) * rng.lognormal(0, 0.12, size=(n, m))
    for c in sorted(set(classes)):            # co-regulated lipid classes
        rows = [i for i, k in enumerate(classes) if k == c]
        act = np.ones(m)
        act[bio] = rng.lognormal(0, 0.18, size=len(bio))
        data[rows, :] *= act
    groups = {g: [i for i, d in enumerate(dx) if d == g] for g in ("Prediabetic", "Type 2 Diabetes",
                                                                   "Metabolic Syndrome")}
    for i, (nm, c) in enumerate(zip(names, classes)):
        eff = CLASS_EFFECT.get(c)
        if eff is None and c in ("PC", "PE") and PUFA.search(nm):
            eff = (0.85, 0.6, 0.65)           # PUFA-containing phospholipids fall
        if eff is None and rng.random() < 0.2:
            f = rng.choice([rng.uniform(1.4, 2.2), rng.uniform(0.45, 0.7)])
            eff = (1 + (f - 1) * 0.5, f, 1 + (f - 1) * 0.8)
        if eff is None:
            continue
        jitter = rng.uniform(0.85, 1.15, size=3)
        for (g, idx), fc, j in zip(groups.items(), eff, jitter):
            data[i, idx] *= fc ** j
    # internal standard row + per-sample injection/extraction factor applied to everything
    istd = 25000 * rng.lognormal(0, 0.05, size=(1, m))
    inj = rng.lognormal(0, 0.12, size=(1, m))
    data = np.vstack([istd, data]) * inj
    # tiered missingness in biological samples (not the ISTD)
    for f in range(1, n + 1):
        t = rng.random()
        rate = rng.uniform(0, 0.15) if t < 0.72 else (rng.uniform(0.2, 0.5) if t < 0.92 else rng.uniform(0.55, 0.85))
        k = int(round(rate * len(bio)))
        if k:
            data[f, rng.choice(bio, size=k, replace=False)] = 0.0

    all_names = [ISTD] + names
    peak = pd.DataFrame(data, columns=samples)
    peak.insert(0, "Metabolite", all_names)
    peak.to_csv(os.path.join(BASE_DIR, "sample_peak_area_matrix_richdemo_lipidomics.csv"), index=False)
    rows = pd.DataFrame({
        "Metabolite": all_names,
        "Lipid_Class": ["Internal standard"] + classes,
        "Category": ["Internal standard"] + [CATEGORY[c] for c in classes],
        "HMDB_ID": [""] + hmdb,
        "ISTD": [ISTD] * (n + 1),
    })
    rows.to_csv(os.path.join(BASE_DIR, "lipid_row_annotations_richdemo.csv"), index=False)
    print(f"[richdemo-lipid] {n} lipids + 1 internal standard x {m} samples; classes: "
          f"{pd.Series(classes).value_counts().to_dict()}")


if __name__ == "__main__":
    generate()
