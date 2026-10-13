"""
build_kegg_pathway_library.py - (Re)build MetaboAI Pro's bundled KEGG metabolite-set library
(Pathway Analysis tab, Option A, library "KEGG") from real KEGG data, calibrated against a real
MetaboAnalyst ORA run.

Why this library replaced the HMDB-annotated one
------------------------------------------------
The previous "KEGG" library grouped HMDB's bulk pathway export by its KEGG_Map_ID column. That is
HMDB's annotation, not KEGG's pathway membership: HMDB files SMPDB pathways under KEGG map IDs, so
e.g. map00250 had 56 members (KEGG / MetaboAnalyst: 28) and map00910 Nitrogen metabolism 25
(MetaboAnalyst: 6), while Steroid hormone biosynthesis had 40 (MetaboAnalyst: 86); nine KEGG maps
MetaboAnalyst tests were missing altogether. That library is still shipped, unchanged, as
"KEGG — HMDB map annotations" (build_hmdb_pathway_libraries.py).

Inputs (all real data; nothing is invented)
-------------------------------------------
1. --kegg      KEGG_human_pathways_compounds_R98.csv: KEGG's own human (hsa) pathway -> compound
               lists, KEGG Release 98 (2021), as bundled in the PyPI package `sspa` 1.0.4 (Wieder
               et al., GPL-3; https://github.com/cwieder/sspa, file
               src/sspa/pathway_databases/KEGG_human_pathways_compounds_R98.csv). One row per hsa
               map: Pathway_name, then KEGG compound IDs.
2. --crosswalk RaMP_v3.0.7_HMDB_KEGG_crosswalk.csv (HMDB_ID, KEGG_ID, HMDB_Name): every HMDB
               accession <-> KEGG compound pair sharing a RaMP-DB compound (RaMP v3.0.7, NCATS,
               data update 2024-08; SQLite from https://github.com/ncats/RaMP-DB, db/
               RaMP_SQLite_v3.0.7.sqlite.gz, `source` table, IDtype hmdb / kegg). Extract with
               extract_ramp_crosswalk() below.
3. --metaboanalyst  a MetaboAnalyst ORA output folder (Rhistory.R, name_map.csv,
               msea_ora_result.csv, ora_membership.json) from metaboanalyst.ca, library
               "kegg_pathway", no reference-metabolome filter. Used ONLY for calibration.

What is built
-------------
metabolomics_modules/data/kegg_human_pathway_library.parquet
    Pathway_ID (hsa#####), Pathway, KEGG_ID (member key: a KEGG compound ID, or an HMDB accession
    for a MetaboAnalyst-verified lipid species -- see (d)), Metabolite, Member_Source.
metabolomics_modules/data/hmdb_kegg_crosswalk.parquet     HMDB_ID -> KEGG_ID (one preferred KEGG ID per HMDB)
metabolomics_modules/data/kegg_metaboanalyst_calibration.json   everything taken from the MetaboAnalyst run.

Steps
  (a) Human metabolic maps = hsa00010-hsa01099 plus hsa01040 (KEGG "Metabolism"; global/overview
      maps 011xx/012xx excluded, as MetaboAnalyst does). Maps with no KEGG compound (C#####)
      members are dropped (hsa00511/00514/00533 are empty and hsa00513/00601/00603/00604 list
      only glycan G-numbers) -> 78 sets.
  (b) Current-KEGG structure: hsa00471 + hsa00472 are merged into hsa00470 "D-Amino acid
      metabolism" (KEGG merged the D-amino-acid maps after Release 98); hsa00860 is titled
      "Porphyrin metabolism" (KEGG's current title). Aminoacyl-tRNA biosynthesis (hsa00970) is
      dropped: MetaboAnalyst's current kegg_pathway library does not contain it (13 of the run's
      in-library amino acids are KEGG members of hsa00970, yet MetaboAnalyst reports no such set).
  (c) Inorganic/currency species removed (water, H+, O2, CO2, Pi, PPi, ions, ... -- the same
      species pathway_libraries.HMDB_CURRENCY_EXCLUDED removes from the HMDB libraries).
  (d) MetaboAnalyst ground-truth membership corrections for the run's tested compounds (the 86
      HMDB IDs MetaboAnalyst mapped): the 60 compounds in MetaboAnalyst's library universe
      (union of fun.anot) are members of EXACTLY the sets ora_membership.json lists for them
      (added where KEGG R98 lacks them, removed from every other set -- a set with hits is always
      reported, so an unreported set cannot contain them); the 26 outside it are removed from
      every set, except lipid species whose MetaboAnalyst KEGG ID is a lipid-CLASS compound
      (Ceramide C00195, Glucosylceramide C01190, Sphingomyelin C00550, PE C00350, Lactosylceramide
      C01290, 1-Acyl-GPC C04230): the class compound stays a member (MetaboAnalyst matches by
      name, so a species simply never matches the class). A verified in-library species whose
      only KEGG ID is such a class (LysoPC(18:1/0:0)) is keyed by its HMDB accession.
  (e) Calibration JSON: MetaboAnalyst's `total` per reported pathway (48), the library-wide
      universe N back-calculated from `expected` = n*total/N (unique solution N = 1519, n = 60 --
      every one of the 48 rows reproduces to 3 significant digits), and the number of sets m = 81
      MetaboAnalyst corrects over (Holm p / raw p of the top rows = 81, 80, ...). The app uses
      these only for "All metabolites in selected library" (MetaboAnalyst's default background).

Usage (from the app folder):
    python3 build_kegg_pathway_library.py --kegg /path/KEGG_human_pathways_compounds_R98.csv \\
        --crosswalk /path/RaMP_v3.0.7_HMDB_KEGG_crosswalk.csv --metaboanalyst /path/ma_output_dir
"""

import argparse
import json
import os
import re

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "metabolomics_modules", "data")
OUT_LIB = os.path.join(OUT_DIR, "kegg_human_pathway_library.parquet")
OUT_XWALK = os.path.join(OUT_DIR, "hmdb_kegg_crosswalk.parquet")
OUT_CAL = os.path.join(OUT_DIR, "kegg_metaboanalyst_calibration.json")

SRC_KEGG = "KEGG Release 98 human pathway compounds (sspa 1.0.4 bundle)"
SRC_ADDED = "MetaboAnalyst-verified member (added)"

MERGE_D_AMINO = {"hsa00471", "hsa00472"}
RENAME = {"hsa00860": "Porphyrin metabolism"}
DROP_MAPS = {"hsa00970": "Aminoacyl-tRNA biosynthesis -- absent from MetaboAnalyst's current kegg_pathway "
                         "library (13 in-library input amino acids are KEGG members; no such set reported)"}

# KEGG IDs of the inorganic / currency species removed from every set (mirrors
# pathway_libraries.HMDB_CURRENCY_EXCLUDED for the HMDB-keyed libraries).
KEGG_CURRENCY_EXCLUDED = {
    "C00001": "H2O", "C00080": "H+", "C00282": "Hydrogen", "C00007": "Oxygen", "C00011": "CO2",
    "C00009": "Orthophosphate", "C00013": "Diphosphate", "C00288": "HCO3-", "C00014": "Ammonia",
    "C01342": "Ammonium", "C00076": "Calcium", "C00305": "Magnesium", "C00034": "Manganese",
    "C00038": "Zinc", "C00238": "Potassium", "C01330": "Sodium", "C00023": "Iron", "C14818": "Fe2+",
    "C14819": "Fe3+", "C00070": "Copper", "C00150": "Molybdenum", "C00698": "Chloride",
    "C00059": "Sulfate", "C00094": "Sulfite", "C00283": "Hydrogen sulfide", "C00027": "Hydrogen peroxide",
    "C00533": "Nitric oxide", "C15602": "Quinone", "C15603": "Hydroquinone",
}
# KEGG lipid-CLASS compounds that MetaboAnalyst's cross-reference assigns to individual lipid species.
KEGG_LIPID_CLASS_IDS = {"C00195": "Ceramide", "C01190": "Glucosylceramide", "C00550": "Sphingomyelin",
                        "C00350": "Phosphatidylethanolamine", "C01290": "Lactosylceramide",
                        "C04230": "1-Acyl-sn-glycero-3-phosphocholine"}
_HMDB_RE = re.compile(r"\s*HMDB(\d{1,7})\s*", re.IGNORECASE)


def _hmdb(v):
    m = _HMDB_RE.fullmatch(str(v))
    return f"HMDB{m.group(1).zfill(7)}" if m else None


def extract_ramp_crosswalk(sqlite_path: str, out_csv: str):
    """HMDB <-> KEGG compound pairs sharing a RaMP compound (rampId) -> CSV (HMDB_ID, KEGG_ID, HMDB_Name)."""
    import sqlite3
    con = sqlite3.connect(sqlite_path)
    s = pd.read_sql("select sourceId, rampId, IDtype, commonName from source where IDtype in ('hmdb','kegg')", con)
    s["id"] = s["sourceId"].str.split(":", n=1).str[1]
    h = s[s["IDtype"] == "hmdb"].copy()
    h["id"] = h["id"].map(_hmdb)
    h = h.dropna(subset=["id"])
    k = s[(s["IDtype"] == "kegg") & s["id"].str.fullmatch(r"C\d{5}")]
    x = h[["id", "rampId"]].drop_duplicates().merge(k[["id", "rampId"]].drop_duplicates(), on="rampId",
                                                     suffixes=("_h", "_k"))
    x = x.rename(columns={"id_h": "HMDB_ID", "id_k": "KEGG_ID"}).drop_duplicates(["HMDB_ID", "KEGG_ID"])
    x["HMDB_Name"] = x["HMDB_ID"].map(h.drop_duplicates("id").set_index("id")["commonName"])
    x.sort_values(["HMDB_ID", "KEGG_ID"]).to_csv(out_csv, index=False)


def read_kegg(path):
    d = pd.read_csv(path, index_col=0)
    name = d["Pathway_name"].str.replace(r" - Homo sapiens \(human\)$", "", regex=True)
    sets, names = {}, {}
    for pid, row in d.drop(columns=["Pathway_name"]).iterrows():
        num = int(pid[3:])
        if not (num < 1100 or pid == "hsa01040"):
            continue
        members = {c for c in row.dropna().astype(str) if re.fullmatch(r"C\d{5}", c)}
        if members:
            sets[pid], names[pid] = members, name[pid]
    merged = set().union(*(sets.pop(p) for p in sorted(MERGE_D_AMINO) if p in sets))
    for p in MERGE_D_AMINO:
        names.pop(p, None)
    sets["hsa00470"], names["hsa00470"] = merged, "D-Amino acid metabolism"
    for p in DROP_MAPS:
        sets.pop(p, None)
        names.pop(p, None)
    names.update({p: n for p, n in RENAME.items() if p in names})
    for p in sets:
        sets[p] -= set(KEGG_CURRENCY_EXCLUDED)
    return sets, names


def read_metaboanalyst(folder):
    nm = pd.read_csv(os.path.join(folder, "name_map.csv"), dtype=str, keep_default_na=False)
    res = pd.read_csv(os.path.join(folder, "msea_ora_result.csv"), index_col=0)
    mem = json.load(open(os.path.join(folder, "ora_membership.json")))
    anot = {k: (v if isinstance(v, list) else [v]) for k, v in mem["fun.anot"].items()}
    rhist = open(os.path.join(folder, "Rhistory.R")).read()
    return nm, res, anot, mem, rhist


def back_calculate(res: pd.DataFrame, n_max: int):
    """Unique (n, N) with expected = n*total/N and phyper reproducing every row (3 significant digits)."""
    def sig3(v):
        return float(f"{v:.3g}")
    sols = []
    ratio = float(np.median(res["total"] / res["expected"]))  # = N / n, up to rounding of `expected`
    for n in range(1, n_max + 1):  # n = query compounds in the library <= mapped queries
        lo, hi = int(n * ratio * 0.99), int(n * ratio * 1.01) + 1
        for N in range(max(n, lo), hi):
            if all(sig3(n * t / N) == e for t, e in zip(res["total"], res["expected"])):
                p = hypergeom.sf(res["hits"] - 1, N, res["total"], n)
                if all(sig3(a) == b for a, b in zip(p, res["Raw p"])):
                    sols.append((n, N))
    ratios = (res["Holm p"] / res["Raw p"]).values
    m = int(round(ratios[0]))
    return sols, m


def main(kegg_path, xwalk_path, ma_dir):
    sets, names = read_kegg(kegg_path)
    print(f"KEGG R98 human metabolic maps after (a)-(c): {len(sets)} sets, "
          f"{len(set().union(*sets.values()))} compounds")

    # ---- crosswalk: one preferred KEGG ID per HMDB accession ----
    xw = pd.read_csv(xwalk_path, dtype=str, keep_default_na=False)
    use = {c: sum(c in s for s in sets.values()) for c in set(xw["KEGG_ID"])}
    xw["_use"] = xw["KEGG_ID"].map(use)
    xw = xw.sort_values(["HMDB_ID", "_use", "KEGG_ID"], ascending=[True, False, True])
    xw_best = xw.drop_duplicates("HMDB_ID")[["HMDB_ID", "KEGG_ID", "HMDB_Name"]].reset_index(drop=True)
    names_by_kegg = xw.sort_values(["KEGG_ID", "HMDB_ID"]).drop_duplicates("KEGG_ID").set_index("KEGG_ID")["HMDB_Name"]

    # ---- MetaboAnalyst ground truth ----
    nm, res, anot, mem, rhist = read_metaboanalyst(ma_dir)
    mapped = nm[nm["Comment"] == "1"].copy()
    in_univ_names = set().union(*map(set, anot.values()))
    name_to_hmdb = dict(zip(mapped["Match"], mapped["HMDB"]))
    verified = {}
    for _, r in mapped.iterrows():
        kid = r["KEGG"] if re.fullmatch(r"C\d{5}", r["KEGG"]) else ""
        inside = r["Match"] in in_univ_names
        key = kid if kid and kid not in KEGG_LIPID_CLASS_IDS else (r["HMDB"] if inside else "")
        verified[r["HMDB"]] = {"name": r["Match"], "kegg": kid, "in_library": inside, "key": key,
                               "pathways": sorted(p for p, v in anot.items() if r["Match"] in v)}
    failed = nm.loc[nm["Comment"] != "1", "Query"].tolist()

    by_name = {n: p for p, n in names.items()}
    missing = [p for p in anot if p not in by_name]
    if missing:
        raise SystemExit(f"MetaboAnalyst pathways with no KEGG map in the library: {missing}")
    in_keys = {v["key"] for v in verified.values() if v["in_library"]}
    out_keys = {v["key"] for v in verified.values() if not v["in_library"] and v["key"]}
    src = {p: {c: SRC_KEGG for c in s} for p, s in sets.items()}
    n_add = n_rm = 0
    for pid in sets:
        want = {verified[name_to_hmdb[n]]["key"] for n in anot.get(names[pid], [])}
        for c in (in_keys | out_keys) & sets[pid] - want:
            sets[pid].discard(c)
            src[pid].pop(c)
            n_rm += 1
        for c in want - sets[pid]:
            sets[pid].add(c)
            src[pid][c] = SRC_ADDED
            n_add += 1
    print(f"MetaboAnalyst corrections: {n_add} memberships added, {n_rm} removed "
          f"({len(in_keys)} in-library / {len(out_keys)} out-of-library verified keys)")

    label = {}
    for h, v in verified.items():
        label.setdefault(v["key"], v["name"])
    rows = []
    for pid in sorted(sets):
        for c in sorted(sets[pid]):
            rows.append({"Pathway_ID": pid, "Pathway": names[pid], "KEGG_ID": c,
                         "Metabolite": label.get(c) or names_by_kegg.get(c, "") or c,
                         "Member_Source": src[pid][c]})
    lib = pd.DataFrame(rows)

    sols, m = back_calculate(res, len(mapped))
    if len(sols) != 1:
        raise SystemExit(f"Universe back-calculation is not unique: {sols}")
    (n_q, N), = sols
    exclude_num = re.search(r'SetCurrentMsetLib\(mSet,\s*"[^"]+",\s*(\d+)\)', rhist)
    cal = {
        "source": "MetaboAnalyst ORA output (metaboanalyst.ca): library kegg_pathway, "
                  "SetMetabolomeFilter(mSet, F), CalculateHyperScore(mSet, 'hyperg')",
        "library_name": mem.get("lib.name"),
        "universe_size": N, "query_in_library": n_q, "n_sets": m,
        "exclude_num": int(exclude_num.group(1)) if exclude_num else None,
        "n_query": int(len(nm)), "n_mapped": int(len(mapped)), "failed_to_map": failed,
        "pathway_totals": {p: int(t) for p, t in res["total"].items()},
        "pathway_ids": {p: by_name[p] for p in res.index},
        "verified_compounds": verified,
        "dropped_maps": DROP_MAPS, "merged_maps": {"hsa00470": sorted(MERGE_D_AMINO)}, "renamed_maps": RENAME,
    }
    os.makedirs(OUT_DIR, exist_ok=True)
    lib.astype("category").to_parquet(OUT_LIB, index=False, compression="zstd")
    xw_best.to_parquet(OUT_XWALK, index=False, compression="zstd")
    with open(OUT_CAL, "w") as f:
        json.dump(cal, f, indent=1, sort_keys=False)
    sizes = lib.groupby("Pathway", observed=True).size()
    print(f"library: {lib['Pathway_ID'].nunique()} sets, {lib['KEGG_ID'].nunique()} unique members, "
          f"{len(lib)} memberships -> {OUT_LIB}")
    print(f"crosswalk: {len(xw_best)} HMDB accessions -> {OUT_XWALK}")
    print(f"calibration: N = {N}, n = {n_q}, m = {m}, {len(cal['pathway_totals'])} pathway totals -> {OUT_CAL}")
    exact = sum(int(sizes.get(p, -1)) == t for p, t in cal["pathway_totals"].items())
    print(f"KEGG R98 (corrected) size == MetaboAnalyst total for {exact} of {len(cal['pathway_totals'])} pathways")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kegg", required=True)
    ap.add_argument("--crosswalk", help="CSV from extract_ramp_crosswalk (or use --ramp-sqlite)")
    ap.add_argument("--ramp-sqlite", help="RaMP SQLite file; the crosswalk CSV is extracted next to it")
    ap.add_argument("--metaboanalyst", required=True)
    a = ap.parse_args()
    xpath = a.crosswalk
    if not xpath:
        if not a.ramp_sqlite:
            raise SystemExit("give --crosswalk or --ramp-sqlite")
        xpath = os.path.splitext(a.ramp_sqlite)[0] + "_HMDB_KEGG_crosswalk.csv"
        extract_ramp_crosswalk(a.ramp_sqlite, xpath)
    main(a.kegg, xpath, a.metaboanalyst)
