"""
Builds the SYNTHETIC example data set (generic format) shipped with Omnix Fluxomics.

Nothing here comes from a real experiment: metabolite names and formulas are taken from the public metabolite
library (data/metabolite_library.csv); every area is simulated from a simple labeling model with random,
seeded parameters. Design: two groups (shControl, GeneKD) × three labeling times (2, 8, 24 h) × 3 replicates,
one unlabeled reference and one blank, [U-13C6]-glucose tracer, data "already corrected" for natural abundance.

Run:  python make_synthetic_example.py      (writes synthetic_isotopologue_table.csv + synthetic_sample_metadata.csv)
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
LIB = HERE.parent.parent / "data" / "metabolite_library.csv"
rng = np.random.default_rng(20261007)

GROUPS = ["shControl", "GeneKD"]
TIMES = [2.0, 8.0, 24.0]
REPS = 3
# pathway-level labeling behaviour (fraction labeled at plateau, rate per h) and the GeneKD effect on it
PATHWAY_MODEL = {
    "GLY": (0.85, 0.9, 0.75), "PPP": (0.6, 0.5, 0.8), "TCA": (0.55, 0.2, 1.25), "HEX": (0.5, 0.15, 0.7),
    "NUC": (0.45, 0.06, 0.85), "PUR": (0.4, 0.06, 0.85), "PYR": (0.35, 0.05, 0.9), "AA": (0.25, 0.08, 1.1),
    "GLU": (0.4, 0.15, 1.2), "GLN": (0.15, 0.1, 1.0), "RED": (0.3, 0.07, 1.15), "OC": (0.2, 0.1, 0.9),
    "UREA": (0.15, 0.05, 1.0), "PL": (0.2, 0.04, 1.0), "FA": (0.15, 0.05, 1.0), "COF": (0.3, 0.04, 0.9),
}
WANT = ["glc", "g6p", "f6p", "fbp", "dhap", "ga3p", "3pg", "pep", "pyr", "lac", "gro3p", "g1p", "6pg", "r5p", "ru5p",
        "s7p", "e4p", "cit", "acon", "akg", "suc", "fum", "mal", "glu", "gln", "asp", "asn", "ala", "ser", "gly",
        "pro", "orn", "citr", "arg", "gsh", "gssg", "imp", "amp", "adp", "atp", "ump", "udp", "utp", "cmp", "oro",
        "udpglcnac", "udpgal", "glcn6p", "nad", "cho", "pcho", "c2car", "car", "met", "sam", "sah", "tau", "pglu",
        "gaba", "naa", "ita", "hyp", "xan", "urate", "ino"]


def carbons(formula: str) -> int:
    m = re.search(r"C(\d*)(?![a-z])", formula)
    return int(m.group(1) or 1) if m else 0


def main():
    L = pd.read_csv(LIB, dtype=str).fillna("")
    L = L[L.LibraryID.isin(WANT)].copy()
    L["n"] = L.Formula.map(carbons)
    L = L[(L.n >= 2) & (L.n <= 25)]
    samples, meta = [], []
    meta.append(dict(Sample="Unlabeled_Ref_1", Type="Unlabeled reference", Group="shControl", Treatment="",
                     Time="", TimeUnit="", Replicate="", Include=True))
    meta.append(dict(Sample="Blank_1", Type="Blank", Group="Blank", Treatment="", Time="", TimeUnit="", Replicate="",
                     Include=True))
    for g in GROUPS:
        for t in TIMES:
            for r in range(1, REPS + 1):
                s = f"{'Ctrl' if g == 'shControl' else g}_{t:g}h_R{r}"
                samples.append((s, g, t))
                meta.append(dict(Sample=s, Type="Labeled", Group=g, Treatment="", Time=t, TimeUnit="h", Replicate=r,
                                 Include=True))
    rows = []
    for m in L.itertuples():
        n = m.n
        code = (m.Pathways.split("|")[0] if m.Pathways else "AA")
        fmax, k, ko = PATHWAY_MODEL.get(code, (0.2, 0.08, 1.0))
        fmax = float(np.clip(fmax * rng.uniform(0.6, 1.25), 0.03, 0.95))
        k = k * rng.uniform(0.6, 1.5)
        ko_eff = ko * rng.uniform(0.9, 1.1)
        q = rng.uniform(0.8, 0.97) if code in ("GLY", "PPP", "HEX") else rng.uniform(0.35, 0.9)  # label completeness
        pool = 10 ** rng.uniform(6.3, 9.3)                 # total pool area (arbitrary units)
        pool_fc = 2 ** rng.normal(0, 0.6) if rng.random() < 0.4 else 1.0   # GeneKD changes some pool sizes
        method = "HILIC_Neg" if rng.random() < 0.55 else "RP_Pos"
        # isotopologue shape of the labeled molecules: binomial(n, q) restricted to ≥ 1 label
        i = np.arange(n + 1)
        from math import comb
        shape = np.array([comb(n, j) * q ** j * (1 - q) ** (n - j) for j in i])
        shape[0] = 0
        shape /= shape.sum()
        vals = {}
        for s, g, t in samples:
            kk = k * (ko_eff if g == "GeneKD" else 1.0)
            ff = min(fmax * (ko_eff ** 0.5 if g == "GeneKD" else 1.0), 0.97)
            lab = ff * (1 - np.exp(-kk * t)) * rng.normal(1, 0.06)
            lab = float(np.clip(lab, 0.002, 0.98))
            mid = (1 - lab) * (i == 0) + lab * shape * rng.lognormal(0, 0.08, n + 1)
            mid = mid / mid.sum()
            tot = pool * (pool_fc if g == "GeneKD" else 1.0) * rng.lognormal(0, 0.15)
            vals[s] = mid * tot
        ref = np.zeros(n + 1); ref[0] = pool * rng.lognormal(0, 0.1); ref[1:] = ref[0] * rng.uniform(1e-4, 8e-4, n)
        blank = np.full(n + 1, np.nan); blank[0] = pool * rng.uniform(1e-4, 2e-3)
        # keep isotopologues that reach ≥ 0.3 % in some labeled sample (others are not reported)
        V = pd.DataFrame(vals)
        frac = V / V.sum(axis=0)
        keep = [j for j in i if j == 0 or frac.loc[j].max() >= 0.003]
        for j in keep:
            row = dict(Metabolite=m.StandardName, Formula=m.Formula, Method=method, KEGG=m.KEGG,
                       Isotopologue=f"M+{j}", Unlabeled_Ref_1=round(float(ref[j]), 0),
                       Blank_1=(round(float(blank[j]), 0) if np.isfinite(blank[j]) else np.nan))
            for s, _, _ in samples:
                row[s] = round(float(max(V.loc[j, s], 1500.0)), 0)
            rows.append(row)
    T = pd.DataFrame(rows)
    T.to_csv(HERE / "synthetic_isotopologue_table.csv", index=False)
    pd.DataFrame(meta).to_csv(HERE / "synthetic_sample_metadata.csv", index=False)
    print(f"{T.Metabolite.nunique()} metabolites, {len(T)} isotopologue rows, {len(samples)} labeled samples")


if __name__ == "__main__":
    main()
