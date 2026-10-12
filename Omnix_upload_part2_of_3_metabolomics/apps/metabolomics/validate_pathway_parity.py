"""
validate_metaboanalyst_parity.py - Reproduce a MetaboAnalyst compound-list ORA with MetaboAI Pro's
Option A and compare pathway by pathway.

Input: a MetaboAnalyst ORA output folder (datalist.csv, name_map.csv, msea_ora_result.csv; the run's
Rhistory.R should show SetMetabolomeFilter(mSet, F) and SetCurrentMsetLib(mSet, "kegg_pathway", 2)).
Runs Option A (run_msea_ora_on_list) with the KEGG library, the whole-library background
("All metabolites in selected library (MetaboAnalyst default)"), min set size 2 (MetaboAnalyst's
excludeNum = 2), hypergeometric test, Holm/FDR across all library sets -- twice:

  parity : only the queries MetaboAnalyst itself mapped (name_map Comment = 1) -- MetaboAnalyst does
           not resolve secondary HMDB accessions, so e.g. HMDB0000357 drops out there;
  app    : every query, with MetaboAI's own mapping (HMDB0000357 -> primary HMDB0000011).

Usage: python3 validate_metaboanalyst_parity.py /path/to/metaboanalyst_output [out.csv]
"""

import os
import sys

import numpy as np
import pandas as pd

from metabolomics_modules import pathway_analysis as pa
from metabolomics_modules import pathway_libraries as libs


def _sig(v, d=3):
    return float(f"{v:.{d}g}") if pd.notna(v) else np.nan


def compare(ma_dir: str):
    queries = [q for q in pd.read_csv(os.path.join(ma_dir, "datalist.csv"), header=None)[0].astype(str)
               if q.strip().lower() != "x"]
    nm = pd.read_csv(os.path.join(ma_dir, "name_map.csv"), dtype=str, keep_default_na=False)
    ma = pd.read_csv(os.path.join(ma_dir, "msea_ora_result.csv"), index_col=0)
    mapped = nm.loc[nm["Comment"] == "1", "Query"].tolist()
    lib = libs.get_library(libs.LIB_KEGG)
    runs = {}
    for label, ids in (("parity", mapped), ("app", queries)):
        res, info, _ = pa.run_msea_ora_on_list(ids, library=lib, reference=pa.REF_LIBRARY, min_size=2,
                                               adjust=pa.ADJUST_ALL)
        runs[label] = (res.set_index("Pathway"), info)
    rows = []
    for p, r in ma.iterrows():
        row = {"Pathway": p, "MA_total": int(r["total"]), "MA_expected": r["expected"], "MA_hits": int(r["hits"]),
               "MA_raw_p": r["Raw p"], "MA_Holm_p": r["Holm p"], "MA_FDR": r["FDR"]}
        for label, (res, _) in runs.items():
            a = res.loc[p] if p in res.index else None
            row.update({f"{label}_total": None if a is None else int(a["Pathway_Size"]),
                        f"{label}_expected": None if a is None else _sig(a["Expected"]),
                        f"{label}_hits": None if a is None else int(a["Overlap_Count"]),
                        f"{label}_raw_p": None if a is None else _sig(a["P-value"]),
                        f"{label}_Holm_p": None if a is None else _sig(a["Holm"]),
                        f"{label}_FDR": None if a is None else _sig(a["FDR"])})
        rows.append(row)
    df = pd.DataFrame(rows)
    extra = {label: sorted(set(res.index) - set(ma.index)) for label, (res, _) in runs.items()}
    return df, {k: v[1] for k, v in runs.items()}, extra


def _match(df, label):
    tot = df[f"{label}_total"] == df["MA_total"]
    hits = df[f"{label}_hits"] == df["MA_hits"]
    p = np.isclose(df[f"{label}_raw_p"].astype(float), df["MA_raw_p"], rtol=0.006)
    fdr = np.isclose(df[f"{label}_FDR"].astype(float), df["MA_FDR"], rtol=0.006, atol=0.0015)
    holm = np.isclose(df[f"{label}_Holm_p"].astype(float), df["MA_Holm_p"], rtol=0.006, atol=0.0015)
    return tot, hits, p, fdr, holm


def main(ma_dir, out_csv=None):
    df, infos, extra = compare(ma_dir)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_rows", 200)
    for label in ("parity", "app"):
        i = infos[label]
        tot, hits, p, fdr, holm = _match(df, label)
        print(f"\n=== {label}: N = {i['N']}, n = {i['n']}, sets tested = {i['n_sets_tested']}, "
              f"correction over m = {i['n_adjusted']}, calibrated = {i['calibrated']}")
        print(f"of {len(df)} MetaboAnalyst pathways: total match {tot.sum()}, hits match {hits.sum()}, "
              f"raw p match {p.sum()}, Holm match {holm.sum()}, FDR match {fdr.sum()}; "
              f"all five {(tot & hits & p & fdr & holm).sum()}")
        print(f"pathways with hits that MetaboAnalyst did not report: {extra[label] or 'none'}")
    cols = ["Pathway", "MA_total", "parity_total", "MA_hits", "parity_hits", "MA_raw_p", "parity_raw_p",
            "MA_FDR", "parity_FDR", "app_hits", "app_raw_p", "app_FDR"]
    print()
    print(df[cols].to_string(index=False))
    if out_csv:
        df.to_csv(out_csv, index=False)
    return df


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
