"""
make_istd_demo_data.py - the internal standard (ISTD) of each targeted demo dataset.

Every targeted dataset has ONE internal standard, used to normalize all of its features
(normalized = feature peak area / ISTD peak area, sample by sample):

  * the peak-area matrix has one ISTD row, spiked at a constant amount (it varies only by
    technical noise, ~5% CV);
  * the row-annotation file has an 'ISTD' column naming that internal standard for every
    metabolite/lipid of the dataset (the ISTD row names itself).

  Dataset                                     Internal standard
  Targeted metabolomics (standard demo)       ISTD_D4-Alanine
  Multi-dataset demo: Targeted Metabolomics   ISTD_D4-Alanine
  Multi-dataset demo: Targeted Lipidomics     ISTD_PC(15:0/18:1-d7)   (a lipid standard)
  Rich clinical demo                          ISTD_D4-Alanine

Analyte values are never changed, and an existing ISTD row with the right name is kept as
it is. The script is idempotent.

    python make_istd_demo_data.py
"""

import os
import zlib

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
METABOLOMICS_ISTD = "ISTD_D4-Alanine"
LIPIDOMICS_ISTD = "ISTD_PC(15:0/18:1-d7)"
FIX_NAMES = {"œâ-hydroxylaurate": "omega-Hydroxylaurate"}   # HMDB mojibake of "ω-hydroxylaurate"


def _is_istd(name) -> bool:
    return "istd" in str(name).lower()


def set_istd_row(matrix_csv, istd, tag):
    """Keep exactly one ISTD row (`istd`) at the top of the matrix. Returns the analyte names."""
    path = os.path.join(HERE, matrix_csv)
    df = pd.read_csv(path)
    df["Metabolite"] = df["Metabolite"].replace(FIX_NAMES)
    existing = df[df["Metabolite"] == istd]
    analytes = df[~df["Metabolite"].map(_is_istd)].reset_index(drop=True)
    if len(existing):
        row = existing.iloc[[0]]
    else:
        rng = np.random.default_rng(zlib.crc32(f"{tag}|{istd}".encode()))
        values = float(np.exp(rng.uniform(9.0, 10.5))) * rng.lognormal(0.0, 0.05, size=df.shape[1] - 1)
        row = pd.DataFrame([[istd, *values]], columns=df.columns)
    pd.concat([row, analytes], ignore_index=True).to_csv(path, index=False)
    return analytes["Metabolite"].tolist()


def annotate(ann, analytes, istd, method=None):
    """Set the ISTD column for `analytes` and add the ISTD's own row (before them)."""
    ann = ann.copy()
    ann["Metabolite"] = ann["Metabolite"].replace(FIX_NAMES)
    if method is not None:
        ann = ann[~((ann["Method"] == method) & ann["Metabolite"].map(_is_istd))]
    else:
        ann = ann[~ann["Metabolite"].map(_is_istd)]
    if "ISTD" not in ann.columns:
        ann["ISTD"] = ""
    ann["ISTD"] = ann["ISTD"].fillna("").astype(str)
    target = ann["Metabolite"].isin(analytes) if method is None else (ann["Method"] == method)
    missing = set(analytes) - set(ann.loc[target, "Metabolite"])
    assert not missing, f"analytes without annotation rows: {sorted(missing)[:5]}"
    ann.loc[target, "ISTD"] = istd
    row = {c: "" for c in ann.columns}
    row.update({"Metabolite": istd, "ISTD": istd})
    if "Class" in ann.columns:
        row["Class"] = "Internal standard"
    if "Method" in ann.columns:
        row["Method"] = method or ann.loc[target, "Method"].mode().iloc[0]
    # the ISTD row goes directly before the first analyte it serves
    first = ann.index[target][0]
    pos = ann.index.get_loc(first)
    return pd.concat([ann.iloc[:pos], pd.DataFrame([row]), ann.iloc[pos:]], ignore_index=True)


def main():
    # 1. Standard targeted metabolomics
    analytes = set_istd_row("sample_peak_area_matrix_targeted.csv", METABOLOMICS_ISTD, "targeted")
    p = os.path.join(HERE, "metabolite_row_annotations_targeted.csv")
    annotate(pd.read_csv(p), analytes, METABOLOMICS_ISTD).to_csv(p, index=False)

    # 2. Multi-dataset demo (one annotation file for all four datasets)
    p = os.path.join(HERE, "metabolite_row_annotations_multimethod.csv")
    ann = pd.read_csv(p)
    for method, istd in (("Targeted Metabolomics", METABOLOMICS_ISTD), ("Targeted Lipidomics", LIPIDOMICS_ISTD)):
        prefix = method.lower().replace(" ", "_")
        analytes = set_istd_row(f"sample_peak_area_matrix_multimethod_{prefix}.csv", istd, f"multimethod_{prefix}")
        ann = annotate(ann, analytes, istd, method=method)
    ann.to_csv(p, index=False)

    # 3. Rich clinical demo (one combined panel, one ISTD)
    analytes = set_istd_row("sample_peak_area_matrix_richdemo.csv", METABOLOMICS_ISTD, "richdemo")
    p = os.path.join(HERE, "metabolite_row_annotations_richdemo.csv")
    ann = annotate(pd.read_csv(p), analytes, METABOLOMICS_ISTD)
    ann.loc[ann["Metabolite"] == METABOLOMICS_ISTD, "Method"] = "Internal standard"
    ann.to_csv(p, index=False)
    print("ISTD rows and ISTD columns written.")


if __name__ == "__main__":
    main()
