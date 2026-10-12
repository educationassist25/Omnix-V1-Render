"""
build_hmdb_gene_mapping.py - (Re)build MetaboAI Pro's bundled HMDB metabolite -> gene/protein table.

The Pathway Analysis tab's Option B (metabolite -> gene -> gene-set ORA) reads its metabolite-gene
associations from   metabolomics_modules/data/hmdb_metabolite_gene_mapping.parquet   (loaded once, lazily, by
metabolomics_modules/hmdb_gene_mapping.py). This script creates that file from a genuine HMDB export, in either
of two input formats:

  1. A tabular export (.xlsx / .csv / .tsv) with the columns
         HMDB_ID, Metabolite, Gene, Protein, UniProt
     one row per metabolite-protein association (e.g. HMDB_Metabolite_Gene_Mapping.xlsx).

  2. The official HMDB metabolite XML (hmdb_metabolites.xml, downloadable from
     https://hmdb.ca/downloads -> "All Metabolites" XML). It is streamed with iterparse (the
     file is several GB; memory stays small). For every <metabolite> the script reads
       <accession>, <name>, <secondary_accessions>/<accession>, and each
       <protein_associations>/<protein> -> <gene_name>, <name>, <uniprot_id>.
     The XML also provides SECONDARY (merged/retired) accessions, which the tabular export does
     not, so this route additionally writes
         metabolomics_modules/data/hmdb_secondary_accessions.parquet   (Secondary_ID -> HMDB_ID)
     which the app picks up automatically to resolve old accessions (e.g. HMDB0000357 ->
     HMDB0000011, 3-hydroxybutyric acid) in users' annotation files.

Usage (from the app folder):
    python3 build_hmdb_gene_mapping.py /path/to/HMDB_Metabolite_Gene_Mapping.xlsx
    python3 build_hmdb_gene_mapping.py /path/to/hmdb_metabolites.xml

Rows are kept verbatim (values stripped of surrounding whitespace; exact duplicate rows dropped;
a gene with no <gene_name> in the XML falls back to the protein name, as in the tabular export).
No filtering by organism or symbol validity is done here -- the app flags non-symbol labels at
runtime instead, so the bundled table stays a faithful copy of HMDB.
"""

import os
import re
import sys
import time

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "metabolomics_modules", "data")
OUT_MAPPING = os.path.join(OUT_DIR, "hmdb_metabolite_gene_mapping.parquet")
OUT_SECONDARY = os.path.join(OUT_DIR, "hmdb_secondary_accessions.parquet")
COLUMNS = ["HMDB_ID", "Metabolite", "Gene", "Protein", "UniProt"]
_HMDB_RE = re.compile(r"\s*HMDB(\d{1,7})\s*", re.IGNORECASE)


def _norm_hmdb(v):
    m = _HMDB_RE.fullmatch(str(v))
    return f"HMDB{m.group(1).zfill(7)}" if m else None


def read_tabular(path: str) -> pd.DataFrame:
    if path.lower().endswith((".xlsx", ".xls")):
        df = pd.read_excel(path, dtype=str, engine="openpyxl")
    else:
        df = pd.read_csv(path, dtype=str, sep=None, engine="python")
    missing = [c for c in COLUMNS if c not in df.columns]
    if missing:
        raise SystemExit(f"Input is missing required column(s): {missing}. Found: {list(df.columns)}")
    return df[COLUMNS]


def read_hmdb_xml(path: str):
    """Stream hmdb_metabolites.xml -> (association DataFrame, secondary-accession DataFrame)."""
    import xml.etree.ElementTree as ET

    rows, sec_rows = [], []
    ns = ""
    for event, elem in ET.iterparse(path, events=("start", "end")):
        if event == "start":
            if not ns and elem.tag.startswith("{"):
                ns = elem.tag.split("}")[0] + "}"
            continue
        if elem.tag != f"{ns}metabolite":
            continue
        acc = _norm_hmdb(elem.findtext(f"{ns}accession", ""))
        name = (elem.findtext(f"{ns}name", "") or "").strip()
        if acc:
            sec = elem.find(f"{ns}secondary_accessions")
            if sec is not None:
                for s in sec.findall(f"{ns}accession"):
                    sid = _norm_hmdb(s.text or "")
                    if sid and sid != acc:
                        sec_rows.append((sid, acc))
            pa = elem.find(f"{ns}protein_associations")
            if pa is not None:
                for p in pa.findall(f"{ns}protein"):
                    pname = (p.findtext(f"{ns}name", "") or "").strip()
                    gene = (p.findtext(f"{ns}gene_name", "") or "").strip() or pname
                    uni = (p.findtext(f"{ns}uniprot_id", "") or "").strip()
                    if gene:
                        rows.append((acc, name, gene, pname, uni))
        elem.clear()
    return (pd.DataFrame(rows, columns=COLUMNS),
            pd.DataFrame(sorted(set(sec_rows)), columns=["Secondary_ID", "HMDB_ID"]))


def clean(df: pd.DataFrame) -> pd.DataFrame:
    df = df.fillna("").astype(str)
    for c in COLUMNS:
        df[c] = df[c].str.strip()
    df["HMDB_ID"] = df["HMDB_ID"].map(_norm_hmdb)
    bad = df["HMDB_ID"].isna() | (df["Gene"] == "")
    if bad.any():
        print(f"  dropping {int(bad.sum())} row(s) with an invalid HMDB accession or empty gene")
    df = df[~bad].drop_duplicates()
    return df.sort_values(["HMDB_ID", "Gene", "UniProt"], kind="mergesort").reset_index(drop=True)


def write_parquet(df: pd.DataFrame, path: str):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    df.astype("category").to_parquet(path, index=False, compression="snappy")


def main(argv):
    if len(argv) != 2:
        raise SystemExit(__doc__)
    src = argv[1]
    t0 = time.time()
    print(f"Reading {src} ...")
    sec = None
    if src.lower().endswith(".xml"):
        df, sec = read_hmdb_xml(src)
    else:
        df = read_tabular(src)
    n_in = len(df)
    df = clean(df)
    write_parquet(df, OUT_MAPPING)
    print(f"  {n_in} input rows -> {len(df)} associations, {df['HMDB_ID'].nunique()} metabolites, "
          f"{df['Gene'].nunique()} distinct gene labels, {df['UniProt'].nunique()} UniProt accessions")
    print(f"  wrote {OUT_MAPPING} ({os.path.getsize(OUT_MAPPING) / 1e6:.2f} MB)")
    if sec is not None:
        sec.to_parquet(OUT_SECONDARY, index=False, compression="snappy")
        print(f"  wrote {OUT_SECONDARY} ({len(sec)} secondary accessions)")
    print(f"Done in {time.time() - t0:.1f} s.")


if __name__ == "__main__":
    main(sys.argv)
