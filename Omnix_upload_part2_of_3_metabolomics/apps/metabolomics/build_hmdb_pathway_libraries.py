"""
build_hmdb_pathway_libraries.py - (Re)build MetaboAI Pro's bundled KEGG and SMPDB metabolite-set
libraries (Pathway Analysis tab, Option A) from a real bulk HMDB pathway export.

Input: HMDB_Metabolite_Pathway_long.csv -- one row per HMDB metabolite-pathway membership, columns
    HMDB_ID, Metabolite, Pathway, SMPDB_ID, KEGG_Map_ID
(the <pathways> section of every HMDB metabocard; the export used for the shipped files has
815,749 rows, 54,282 metabolites, 49,628 pathway names).

Output (metabolomics_modules/data/, loaded once, lazily, by metabolomics_modules/pathway_libraries.py):
    kegg_pathway_library.parquet    Pathway_ID, Pathway, HMDB_Pathway_Names, ID_Source, HMDB_ID, Metabolite
                                    (served as the library "KEGG — HMDB map annotations (previous library)";
                                    the main "KEGG" library is now built from KEGG's own pathway compound
                                    lists by build_kegg_pathway_library.py, because HMDB's KEGG_Map_ID
                                    annotation does not reproduce KEGG / MetaboAnalyst set sizes)
    smpdb_pathway_library.parquet   Pathway_ID, Pathway, Category, KEGG_Map_ID, HMDB_ID, Metabolite
    hmdb_pathway_id_crosswalk.parquet  HMDB_ID, ID_Type ("SMPDB" | "KEGG_Map"), ID -- EVERY distinct
                                    HMDB_ID -> SMPDB_ID and HMDB_ID -> KEGG_Map_ID pair of the complete
                                    export (no category or species-variant filtering; SMPDB IDs in the
                                    7-digit form), used by the MSEA tab's ID Mapping table

How the two libraries are derived (all from the one file; nothing hand-curated except the
KEGG map titles and the SMPDB category rules below)
------------------------------------------------------------------------------------------
KEGG  = every row with a KEGG_Map_ID, grouped by KEGG_Map_ID (NOT by pathway name: HMDB lists
        several names under one map, e.g. map00250 = "Alanine, aspartate and glutamate
        metabolism" + SMPDB "Aspartate Metabolism" + SMPDB "Glutamate Metabolism"; no name maps
        to two map IDs). The set is the union of HMDB's memberships for that map. The display
        name is KEGG's official map title (KEGG_TITLES; three maps HMDB still cites --
        map00072, map00150, map00472 -- are retired in current KEGG and are labelled so); the
        HMDB names are kept in HMDB_Pathway_Names. The only rows with neither ID (29 rows,
        "pentose phosphate pathway") carry KEGG's exact map title for map00030 and are
        assigned to it (ID_Source = "name match to KEGG title").
SMPDB = every row with an SMPDB_ID, grouped by SMPDB_ID (IDs normalised to SMPDB's current
        7-digit form, SMP00057 -> SMP0000057; each ID has exactly one name). Each pathway gets
        an SMPDB-style Category from its name (classify() below; SMPDB's own type field is
        not in the export):
          Metabolic, Disease, Drug action, Drug metabolism, Signaling / physiological, and
          Metabolic (lipid/acylcarnitine species-specific) -- the 48,810 SMPDB pathways that are
          one-lipid-species variants of 4 parent pathways ("Cardiolipin Biosynthesis
          CL(...)", "De Novo Triacylglycerol Biosynthesis TG(...)", "Phosphatidylcholine /
          Phosphatidylethanolamine Biosynthesis PC/PE(...)") or single "Acylcarnitine X"
          pathways. Those species-specific variants are NOT written: each shares its parent's
          backbone metabolites and differs by one lipid species, so testing them would add
          ~49k near-duplicate sets (tens of thousands with identical members inside a
          dataset's detected background), bury real signal under the multiple-testing
          correction and inflate the table. The four generic parent pathways (SMP0014212,
          SMP0015896, SMP0020986, SMP0029731) ARE kept, as Metabolic.
        The four pathway names that occur under two SMPDB IDs get " [SMP…]" appended so set
        names stay unique.

Usage (from the app folder):
    python3 build_hmdb_pathway_libraries.py /path/to/HMDB_Metabolite_Pathway_long.csv
"""

import os
import re
import sys
import time

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "metabolomics_modules", "data")
OUT_KEGG = os.path.join(OUT_DIR, "kegg_pathway_library.parquet")
OUT_SMPDB = os.path.join(OUT_DIR, "smpdb_pathway_library.parquet")
OUT_XWALK = os.path.join(OUT_DIR, "hmdb_pathway_id_crosswalk.parquet")
COLUMNS = ["HMDB_ID", "Metabolite", "Pathway", "SMPDB_ID", "KEGG_Map_ID"]

# Official KEGG map titles (rest.kegg.jp/list/pathway) for every map ID the HMDB export cites.
KEGG_TITLES = {
    "map00010": "Glycolysis / Gluconeogenesis", "map00020": "Citrate cycle (TCA cycle)",
    "map00030": "Pentose phosphate pathway", "map00051": "Fructose and mannose metabolism",
    "map00052": "Galactose metabolism", "map00062": "Fatty acid elongation",
    "map00071": "Fatty acid degradation",
    "map00072": "Synthesis and degradation of ketone bodies (retired KEGG map)",
    "map00100": "Steroid biosynthesis", "map00120": "Primary bile acid biosynthesis",
    "map00130": "Ubiquinone and other terpenoid-quinone biosynthesis",
    "map00140": "Steroid hormone biosynthesis",
    "map00150": "Androgen and estrogen metabolism (retired KEGG map)",
    "map00190": "Oxidative phosphorylation", "map00230": "Purine metabolism",
    "map00232": "Caffeine metabolism", "map00240": "Pyrimidine metabolism",
    "map00250": "Alanine, aspartate and glutamate metabolism",
    "map00260": "Glycine, serine and threonine metabolism", "map00270": "Cysteine and methionine metabolism",
    "map00280": "Valine, leucine and isoleucine degradation", "map00310": "Lysine degradation",
    "map00330": "Arginine and proline metabolism", "map00340": "Histidine metabolism",
    "map00350": "Tyrosine metabolism", "map00360": "Phenylalanine metabolism",
    "map00380": "Tryptophan metabolism", "map00410": "beta-Alanine metabolism",
    "map00430": "Taurine and hypotaurine metabolism", "map00450": "Selenocompound metabolism",
    "map00472": "D-Arginine and D-ornithine metabolism (retired KEGG map)",
    "map00480": "Glutathione metabolism", "map00500": "Starch and sucrose metabolism",
    "map00520": "Amino sugar and nucleotide sugar metabolism", "map00561": "Glycerolipid metabolism",
    "map00562": "Inositol phosphate metabolism", "map00564": "Glycerophospholipid metabolism",
    "map00590": "Arachidonic acid metabolism", "map00592": "alpha-Linolenic acid metabolism",
    "map00620": "Pyruvate metabolism", "map00640": "Propanoate metabolism",
    "map00650": "Butanoate metabolism", "map00670": "One carbon pool by folate",
    "map00730": "Thiamine metabolism", "map00740": "Riboflavin metabolism",
    "map00750": "Vitamin B6 metabolism", "map00760": "Nicotinate and nicotinamide metabolism",
    "map00770": "Pantothenate and CoA biosynthesis", "map00780": "Biotin metabolism",
    "map00790": "Folate biosynthesis", "map00830": "Retinol metabolism",
    "map00860": "Porphyrin metabolism", "map00910": "Nitrogen metabolism",
    "map00920": "Sulfur metabolism", "map01040": "Biosynthesis of unsaturated fatty acids",
}

# ---- SMPDB category rules (name-based; documented in the module docstring) ----
CAT_METABOLIC = "Metabolic"
CAT_SPECIES = "Metabolic (lipid/acylcarnitine species-specific)"
CAT_DISEASE = "Disease"
CAT_DRUG_ACTION = "Drug action"
CAT_DRUG_METAB = "Drug metabolism"
CAT_SIGNAL = "Signaling / physiological"

SPECIES_FAMILIES = ("Cardiolipin Biosynthesis", "De Novo Triacylglycerol Biosynthesis",
                    "Phosphatidylcholine Biosynthesis", "Phosphatidylethanolamine Biosynthesis")
_SPECIES_RE = re.compile(r"^(%s)\s+\S+\(.*\)\s*$" % "|".join(re.escape(f) for f in SPECIES_FAMILIES))
_DRUG_ACTION_RE = re.compile(r"Action Pathway|Antihistamine Action|Inhibition of BCR-ABL", re.I)
_DRUG_METAB_RE = re.compile(r"Metabolism Pathway\s*$", re.I)
_DISEASE_RE = re.compile(
    r"Deficien|uria\b|urias?\b|acidura\b|emia\b|aemia\b|Disease|Syndrome|Disorder|Intolerance|Porphyri|"
    r"Leukodystrophy|Glycogenosis|Galactosemia|Xanthomatosis|Encephalopathy|Hyperplasia|acidosis|Oncogenic|"
    r"Cancer|Oncometabolite|Defect|Dystonia|malabsorption|Cystinosis|Gangliosidosis|Mucopolysaccharidosis|"
    r"dependency|Hypophosphatasia|Refsum|Desmosterolosis|Homocarnosinosis|Hyperinsulinism|Hypoacetylaspartia|"
    r"Chondrodysplasia|Barth|excess|^Triosephosphate isomerase$|Hyperlysinemia|Hyperglycinemia|Carnosinuria|"
    r"Tyrosinemia|Hyperphenylalan", re.I)
_SIGNAL_RE = re.compile(r"Signal|Receptor|Activation|Function$|Contraction|Coagulation|Production|Regulation|"
                        r"Transcription|Replication|Stimulation", re.I)
_HMDB_RE = re.compile(r"\s*HMDB(\d{1,7})\s*", re.IGNORECASE)
_SMP_RE = re.compile(r"\s*SMP(\d{1,7})\s*", re.IGNORECASE)


def classify_smpdb(name: str) -> str:
    """SMPDB-style pathway category from the pathway name (first matching rule wins)."""
    n = str(name).strip()
    if _SPECIES_RE.match(n) or n.startswith("Acylcarnitine "):
        return CAT_SPECIES
    if _DRUG_ACTION_RE.search(n):
        return CAT_DRUG_ACTION
    if _DRUG_METAB_RE.search(n):
        return CAT_DRUG_METAB
    if _DISEASE_RE.search(n):
        return CAT_DISEASE
    if _SIGNAL_RE.search(n):
        return CAT_SIGNAL
    return CAT_METABOLIC


def _norm_hmdb(v):
    m = _HMDB_RE.fullmatch(str(v))
    return f"HMDB{m.group(1).zfill(7)}" if m else None


def _norm_smp(v):
    m = _SMP_RE.fullmatch(str(v))
    return f"SMP{m.group(1).zfill(7)}" if m else ""


def read_export(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, dtype=str, keep_default_na=False)
    missing = [c for c in COLUMNS if c not in df.columns]
    if missing:
        raise SystemExit(f"{path}: missing column(s) {missing}; expected {COLUMNS}")
    df = df[COLUMNS].apply(lambda s: s.str.strip())
    df["HMDB_ID"] = df["HMDB_ID"].map(_norm_hmdb)
    return df[df["HMDB_ID"].notna() & (df["Pathway"] != "")]


def build_kegg(df: pd.DataFrame) -> pd.DataFrame:
    k = df[df["KEGG_Map_ID"] != ""].copy()
    k["ID_Source"] = "HMDB KEGG_Map_ID"
    title_to_map = {v.lower(): m for m, v in KEGG_TITLES.items()}
    noid = df[(df["KEGG_Map_ID"] == "") & (df["SMPDB_ID"] == "")].copy()
    noid["KEGG_Map_ID"] = noid["Pathway"].str.lower().map(title_to_map).fillna("")
    noid = noid[noid["KEGG_Map_ID"] != ""]
    noid["ID_Source"] = "name match to KEGG title"
    k = pd.concat([k, noid], ignore_index=True)
    unknown = sorted(set(k["KEGG_Map_ID"]) - set(KEGG_TITLES))
    if unknown:
        raise SystemExit(f"KEGG map IDs without a title in KEGG_TITLES: {unknown} -- add them.")
    names = k.groupby("KEGG_Map_ID")["Pathway"].agg(lambda s: "; ".join(sorted(set(s))))
    src = k.groupby("KEGG_Map_ID")["ID_Source"].agg(lambda s: "; ".join(sorted(set(s))))
    out = k.drop_duplicates(["KEGG_Map_ID", "HMDB_ID"])[["KEGG_Map_ID", "HMDB_ID", "Metabolite"]]
    out = out.rename(columns={"KEGG_Map_ID": "Pathway_ID"})
    out.insert(1, "Pathway", out["Pathway_ID"].map(KEGG_TITLES))
    out.insert(2, "HMDB_Pathway_Names", out["Pathway_ID"].map(names))
    out.insert(3, "ID_Source", out["Pathway_ID"].map(src))
    return out.sort_values(["Pathway_ID", "HMDB_ID"]).reset_index(drop=True)


def build_smpdb(df: pd.DataFrame):
    s = df[df["SMPDB_ID"] != ""].copy()
    s["Pathway_ID"] = s["SMPDB_ID"].map(_norm_smp)
    s = s[s["Pathway_ID"] != ""]
    per_id = s.drop_duplicates("Pathway_ID")[["Pathway_ID", "Pathway"]].copy()
    per_id["Category"] = per_id["Pathway"].map(classify_smpdb)
    cat_counts_all = per_id["Category"].value_counts().to_dict()
    dup = per_id["Pathway"].duplicated(keep=False)
    per_id.loc[dup, "Pathway"] = per_id.loc[dup, "Pathway"] + " [" + per_id.loc[dup, "Pathway_ID"] + "]"
    kegg_x = s[s["KEGG_Map_ID"] != ""].groupby("Pathway_ID")["KEGG_Map_ID"].agg(
        lambda x: "; ".join(sorted(set(x))))
    per_id["KEGG_Map_ID"] = per_id["Pathway_ID"].map(kegg_x).fillna("")
    keep = per_id[per_id["Category"] != CAT_SPECIES].set_index("Pathway_ID")
    out = s[s["Pathway_ID"].isin(keep.index)].drop_duplicates(["Pathway_ID", "HMDB_ID"])
    out = out[["Pathway_ID", "HMDB_ID", "Metabolite"]].copy()
    for c in ("Pathway", "Category", "KEGG_Map_ID"):
        out[c] = out["Pathway_ID"].map(keep[c])
    out = out[["Pathway_ID", "Pathway", "Category", "KEGG_Map_ID", "HMDB_ID", "Metabolite"]]
    return out.sort_values(["Pathway_ID", "HMDB_ID"]).reset_index(drop=True), cat_counts_all


def build_id_crosswalk(df: pd.DataFrame) -> pd.DataFrame:
    """Every distinct (HMDB_ID, SMPDB_ID) and (HMDB_ID, KEGG_Map_ID) pair of the complete export."""
    smp = pd.DataFrame({"HMDB_ID": df["HMDB_ID"], "ID_Type": "SMPDB", "ID": df["SMPDB_ID"].map(_norm_smp)})
    kg = pd.DataFrame({"HMDB_ID": df["HMDB_ID"], "ID_Type": "KEGG_Map", "ID": df["KEGG_Map_ID"]})
    out = pd.concat([smp, kg], ignore_index=True)
    out = out[out["ID"] != ""].drop_duplicates()
    return out.sort_values(["HMDB_ID", "ID_Type", "ID"]).reset_index(drop=True)


def _write(df: pd.DataFrame, path: str):
    df = df.copy()
    for c in df.columns:
        df[c] = df[c].astype("category")
    df.to_parquet(path, index=False, compression="zstd")


def main(path: str):
    t0 = time.time()
    df = read_export(path)
    print(f"read {len(df):,} rows / {df['HMDB_ID'].nunique():,} metabolites / "
          f"{df['Pathway'].nunique():,} pathway names in {time.time() - t0:.1f} s")
    os.makedirs(OUT_DIR, exist_ok=True)
    kegg = build_kegg(df)
    _write(kegg, OUT_KEGG)
    print(f"KEGG : {kegg['Pathway_ID'].nunique()} maps, {kegg['HMDB_ID'].nunique():,} metabolites, "
          f"{len(kegg):,} memberships -> {OUT_KEGG} ({os.path.getsize(OUT_KEGG) / 1024:.0f} KB)")
    smp, cats = build_smpdb(df)
    _write(smp, OUT_SMPDB)
    print("SMPDB categories in the export (pathways):", cats)
    kept = smp.drop_duplicates("Pathway_ID")["Category"].value_counts().to_dict()
    print(f"SMPDB: {smp['Pathway_ID'].nunique()} pathways written {kept}, {smp['HMDB_ID'].nunique():,} "
          f"metabolites, {len(smp):,} memberships -> {OUT_SMPDB} ({os.path.getsize(OUT_SMPDB) / 1024:.0f} KB)")
    xw = build_id_crosswalk(df)
    _write(xw, OUT_XWALK)
    print(f"ID crosswalk: {xw['HMDB_ID'].nunique():,} metabolites, "
          f"{int((xw['ID_Type'] == 'SMPDB').sum()):,} HMDB-SMPDB and {int((xw['ID_Type'] == 'KEGG_Map').sum()):,} "
          f"HMDB-KEGG map pairs -> {OUT_XWALK} ({os.path.getsize(OUT_XWALK) / 1024:.0f} KB)")
    met = smp[smp["Category"] == CAT_METABOLIC]
    print(f"SMPDB Metabolic: {met['Pathway_ID'].nunique()} pathways, {met['HMDB_ID'].nunique():,} metabolites")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    main(sys.argv[1])
