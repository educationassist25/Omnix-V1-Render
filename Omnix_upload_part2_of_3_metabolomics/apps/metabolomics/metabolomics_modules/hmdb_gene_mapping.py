"""
hmdb_gene_mapping.py - Real HMDB metabolite -> gene / protein / UniProt associations for the
Pathway Analysis tab's Option B (metabolite -> gene -> gene-set ORA).

Source (authoritative)
----------------------
metabolomics_modules/data/hmdb_metabolite_gene_mapping.parquet -- a verbatim copy of a genuine HMDB export
(HMDB_Metabolite_Gene_Mapping.xlsx: columns HMDB_ID, Metabolite, Gene, Protein, UniProt; one row
per HMDB metabolite-protein association, i.e. the "Enzymes / Transporters" of each metabocard):
858,077 associations covering 22,849 HMDB metabolites, 6,124 gene labels and 7,281 UniProt
accessions. The FULL table is shipped (not a demo subset), so users' own datasets map too.
Example: HMDB0004952 Cer(d18:1/22:0) -> 70 associations (70 UniProt proteins, 69 distinct gene
symbols -- PLEKHA8 is linked through two UniProt entries), exactly as on hmdb.ca.
Rebuild it from a newer export -- or from the official hmdb_metabolites.xml, which additionally
yields secondary accessions -- with build_hmdb_gene_mapping.py in the app folder.

This REPLACES the small hand-curated metabolite->gene lists that pathway_reference_db.py used to
carry (typically 2-9 "representative" genes per metabolite, e.g. only CERS2/ASAH1/SGMS1 for
HMDB0004952); those lists have been removed.

What this module adds on top of the raw table (all transparent, none of it changes the data)
------------------------------------------------------------------------------------------
* Lazy, load-once access (module-level cache; the app has no st.cache_* convention): the table
  is read on first use (~0.1 s, 1.2 MB file) and indexed by HMDB accession.
* Secondary-accession resolution: HMDB merges/retires accessions. A dataset can therefore carry
  an old ID that has no row in the export (the Rich Clinical Demo's 3-hydroxybutyric acid is
  HMDB0000357, a secondary accession of HMDB0000011 -- verified on hmdb.ca). VERIFIED_SECONDARY
  lists the ones checked by hand; if metabolomics_modules/data/hmdb_secondary_accessions.parquet exists (built
  from hmdb_metabolites.xml) every HMDB secondary accession is resolved.
* Name index over all 22,849 HMDB common names (case/punctuation-insensitive), used as a last
  resort when a metabolite has no HMDB ID from the annotation table or the curated reference.
* Gene-label quality flag: a few HMDB "gene" values are not gene symbols (e.g. microbial
  '(REFSEQ) ... ' descriptions, 'UNC13B variant protein', 'ppar gamma2'); is_gene_symbol()
  identifies them so Option B can keep them in the mapping table but out of the ORA input.
* Symbol harmonization for set matching: HMDB still uses some retired HGNC symbols (GBA, MUT,
  PPAP2A, COL4A3BP, ATP5A1, ...) whereas gene-set collections may use the current ones (GBA1,
  MMUT, PLPP1, CERT1, ATP5F1A, ...) -- or, like MSigDB v7.0 Hallmark, a mix. RETIRED_TO_CURRENT is
  a UniProt-anchored alias table (each retired symbol occurs in the HMDB export; the UniProt
  accession and protein name given there identify the gene). symbol_for_collection() picks,
  per collection, whichever form the collection actually contains; the HMDB symbol is always
  shown unchanged next to it.
"""

import os
import re

import numpy as np
import pandas as pd

_DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
MAPPING_PATH = os.path.join(_DATA_DIR, "hmdb_metabolite_gene_mapping.parquet")
SECONDARY_PATH = os.path.join(_DATA_DIR, "hmdb_secondary_accessions.parquet")

SOURCE_LABEL = ("HMDB metabolite-protein associations (bundled export "
                "HMDB_Metabolite_Gene_Mapping.xlsx -> metabolomics_modules/data/hmdb_metabolite_gene_mapping.parquet)")
EXPECTED_COUNTS = {"rows": 858077, "metabolites": 22849}  # of the bundled export (sanity check)

# Secondary (merged) accession -> primary accession, each checked on hmdb.ca.
VERIFIED_SECONDARY = {
    "HMDB0000357": "HMDB0000011",  # 3-Hydroxybutyric acid (hmdb.ca/metabolites/HMDB0000011 lists HMDB0000357)
}

# Retired HGNC symbol used by HMDB -> current HGNC symbol. Each retired symbol occurs in the bundled
# HMDB export, and the UniProt accession listed is the one HMDB gives for it (verified against the
# bundled table; 6 are unreviewed TrEMBL entries named after the same gene, e.g. Q6IBN4 'PECI protein').
RETIRED_TO_CURRENT = {
    # symbol: (current symbol, UniProt accession in the HMDB export)
    "GBA": ("GBA1", "P04062"), "MUT": ("MMUT", "P22033"), "G6PC": ("G6PC1", "P35575"),
    "ADSS": ("ADSS2", "P30520"), "ADSSL1": ("ADSS1", "Q8N142"), "AGPAT6": ("GPAT4", "Q86UL3"),
    "AGPAT9": ("GPAT3", "Q53EU6"), "NT5C3": ("NT5C3A", "Q9H0P0"), "EPT1": ("SELENOI", "Q9C0D9"),
    "PPAP2A": ("PLPP1", "O14494"), "PPAP2C": ("PLPP2", "O43688"), "PPAP2B": ("PLPP3", "O14495"),
    "PPAPDC1A": ("PLPP4", "Q5VZY2"), "PPAPDC1B": ("PLPP5", "Q8NEB5"), "PPAPDC2": ("PLPP6", "Q8IY26"),
    "PPAPDC3": ("PLPP7", "Q8NBV4"), "LPPR4": ("PLPPR4", "Q7Z2D5"), "LPPR5": ("PLPPR5", "Q32ZL2"),
    "COL4A3BP": ("CERT1", "Q9Y5P4"), "GLTPD1": ("CPTP", "Q5TA50"), "ABP1": ("AOC1", "P19801"),
    "ATP5A1": ("ATP5F1A", "P25705"), "ATP5B": ("ATP5F1B", "P06576"), "ATP5C1": ("ATP5F1C", "P36542"),
    "ATP5D": ("ATP5F1D", "P30049"), "ATP5E": ("ATP5F1E", "P56381"), "ATP5O": ("ATP5PO", "P48047"),
    "ATP5F1": ("ATP5PB", "P24539"), "ATP5H": ("ATP5PD", "O75947"), "ATP5J": ("ATP5PF", "P18859"),
    "ATP5G1": ("ATP5MC1", "P05496"), "ATP5G3": ("ATP5MC3", "P48201"),
    "ATP5I": ("ATP5ME", "P56385"), "ATP5J2": ("ATP5MF", "P56134"),
    "PECI": ("ECI2", "Q6IBN4"), "ADFP": ("PLIN2", "Q6FHZ7"),
    "ERO1L": ("ERO1A", "Q96HE7"), "ERO1LB": ("ERO1B", "Q86YB8"), "CTGF": ("CCN2", "P29279"),
    "DAK": ("TKFC", "Q3LXA3"), "FUK": ("FCSK", "Q8N0W3"), "CCBL1": ("KYAT1", "Q16773"),
    "CCBL2": ("KYAT3", "Q6YP21"), "GUCY1A3": ("GUCY1A1", "Q02108"), "GUCY1B3": ("GUCY1B1", "Q02153"),
    "SQRDL": ("SQOR", "Q9Y6N5"), "LEPRE1": ("P3H1", "Q32P28"), "LEPREL1": ("P3H2", "Q8IVL5"),
    "LEPREL2": ("P3H3", "Q8IVL6"), "ADRBK1": ("GRK2", "P25098"), "ADRBK2": ("GRK3", "P35626"),
    "MGEA5": ("OGA", "O60502"), "NAPRT1": ("NAPRT", "Q6XQN6"), "CECR1": ("ADA2", "Q9NZK5"),
    "PROSC": ("PLPBP", "O94903"), "AGXT2L1": ("ETNPPL", "Q8TBG4"), "AGXT2L2": ("PHYKPL", "Q8IUZ5"),
    "CARKD": ("NAXD", "Q8IW45"), "SC5DL": ("SC5D", "O75845"),
    "PLA2G16": ("PLAAT3", "P53816"), "HRASLS2": ("PLAAT2", "Q9NWW9"), "GLT25D1": ("COLGALT1", "Q8NBJ5"),
    "GLT25D2": ("COLGALT2", "Q8IYK4"), "EPB49": ("DMTN", "Q08495"), "SKIV2L2": ("MTREX", "P42285"),
    "MLTK": ("MAP3K20", "Q9NYL2"), "MLL": ("KMT2A", "Q03164"), "MLL2": ("KMT2D", "O14686"),
    "MLL3": ("KMT2C", "Q8NEZ4"), "WBP7": ("KMT2B", "Q9UMN6"), "WHSC1": ("NSD2", "O96028"),
    "WHSC1L1": ("NSD3", "Q9BZ95"), "MYST3": ("KAT6A", "A5PKX7"), "PFTK1": ("CDK14", "O94921"),
    "MST4": ("STK26", "Q9P289"), "PIK4CB": ("PI4KB", "Q9UBF8"), "IHPK3": ("IP6K3", "Q5TAQ4"),
    "PAK7": ("PAK5", "Q9P286"), "MB21D1": ("CGAS", "Q8N884"), "RFWD2": ("COP1", "Q05CT6"),
    "TENC1": ("TNS2", "Q63HR2"), "ACPP": ("ACP3", "P15309"), "ALPPL2": ("ALPG", "P10696"),
    "AGPHD1": ("HYKK", "A2RU49"), "PEO1": ("TWNK", "Q96RR1"), "ATPBD4": ("DPH6", "Q7L8W6"),
    "C22orf28": ("RTCB", "Q9Y3I0"), "NADKD1": ("NADK2", "Q4G0N4"), "NO66": ("RIOX1", "Q9H6W3"),
    "MINA": ("RIOX2", "Q8IUF8"), "PYCRL": ("PYCR3", "Q53H96"), "MOSC2": ("MTARC2", "Q969Z3"),
    "PET112": ("GATB", "O75879"), "QTRTD1": ("QTRT2", "Q9H974"), "WBSCR17": ("GALNT17", "Q6IS24"),
    "WBSCR22": ("BUD23", "O43709"), "C5orf4": ("FAXDC2", "Q96IV6"), "PAPD4": ("TENT2", "Q6PIY7"),
    "TUBB2C": ("TUBB4B", "P68371"), "TUBB4": ("TUBB4A", "P04350"), "COPG": ("COPG1", "Q9Y678"),
    "SEPT10": ("SEPTIN10", "Q9P0V9"), "HMHA1": ("ARHGAP45", "Q92619"), "FARSLA": ("FARSA", "Q6IBR2"),
    "CARS": ("CARS1", "P49589"), "LARS": ("LARS1", "Q9P2J5"), "MARS": ("MARS1", "P56192"),
    "NARS": ("NARS1", "O43776"), "QARS": ("QARS1", "P47897"), "SARS": ("SARS1", "P49591"),
    "VARS": ("VARS1", "P26640"), "YARS": ("YARS1", "P54577"), "HARS": ("HARS1", "P12081"),
    "GARS": ("GARS1", "P41250"), "TARS": ("TARS1", "P26639"), "IARS": ("IARS1", "P41252"),
    "AARS": ("AARS1", "P49588"), "DARS": ("DARS1", "P14868"), "EPRS": ("EPRS1", "P07814"),
    "KARS": ("KARS1", "Q15046"), "RARS": ("RARS1", "P54136"), "WARS": ("WARS1", "P23381"),
    "TARSL2": ("TARS3", "A2RTX5"),
}
CURRENT_TO_RETIRED = {}
for _old, (_new, _u) in RETIRED_TO_CURRENT.items():
    CURRENT_TO_RETIRED.setdefault(_new, []).append(_old)

_SYMBOL_RE = re.compile(r"[A-Za-z0-9][A-Za-z0-9\-\.@_]*")
_HMDB_RE = re.compile(r"\s*HMDB(\d{1,7})\s*", re.IGNORECASE)

_STATE = {}


def normalize_name(name) -> str:
    return re.sub(r"[^a-z0-9]", "", str(name).lower())


def normalize_hmdb_id(value):
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return None
    m = _HMDB_RE.fullmatch(str(value))
    return f"HMDB{m.group(1).zfill(7)}" if m else None


def is_gene_symbol(label) -> bool:
    """True for an HGNC-style symbol (e.g. 'CERS2', 'MT-ND1', 'DKFZp434L0435'); False for free-text labels."""
    return bool(_SYMBOL_RE.fullmatch(str(label)))


def _load():
    """Read the bundled table once and build the indexes (module-level cache)."""
    if "table" in _STATE:
        return _STATE
    if not os.path.exists(MAPPING_PATH):
        raise FileNotFoundError(
            f"Bundled HMDB gene-mapping table not found at {MAPPING_PATH}. Rebuild it with "
            "build_hmdb_gene_mapping.py (see that script's docstring).")
    # Columns stay categorical (dictionary-encoded, as stored): converting all 858k rows to Python
    # strings would cost ~2 s; only the small per-metabolite slices are converted, on access.
    df = pd.read_parquet(MAPPING_PATH)
    rows_by_id = df.groupby("HMDB_ID", sort=False, observed=True).indices
    names = df.drop_duplicates("HMDB_ID").astype(str).set_index("HMDB_ID")["Metabolite"]
    name_index = {}
    for hid, nm in names.items():
        name_index.setdefault(normalize_name(nm), []).append(hid)
    secondary = dict(VERIFIED_SECONDARY)
    if os.path.exists(SECONDARY_PATH):
        sec = pd.read_parquet(SECONDARY_PATH)
        secondary.update(dict(zip(sec["Secondary_ID"].astype(str), sec["HMDB_ID"].astype(str))))
    _STATE.update({"table": df, "rows_by_id": rows_by_id, "names": names.to_dict(),
                   "name_index": name_index, "secondary": secondary})
    return _STATE


def table() -> pd.DataFrame:
    """The full association table (HMDB_ID, Metabolite, Gene, Protein, UniProt; categorical dtype)."""
    return _load()["table"]


def summary() -> dict:
    s = _load()
    t = s["table"]
    return {"associations": len(t), "metabolites": len(s["rows_by_id"]), "genes": t["Gene"].nunique(),
            "uniprot": t["UniProt"].nunique(), "secondary_accessions": len(s["secondary"]),
            "source": SOURCE_LABEL}


def resolve_accession(hmdb_id):
    """(accession to use, note). Primary accessions pass through; known secondary ones are resolved."""
    hid = normalize_hmdb_id(hmdb_id)
    if hid is None:
        return None, ""
    s = _load()
    if hid in s["rows_by_id"]:
        return hid, ""
    prim = s["secondary"].get(hid)
    if prim:
        return prim, f"{hid} is an HMDB secondary accession of {prim}; using {prim}"
    return hid, ""


def has_associations(hmdb_id) -> bool:
    return hmdb_id in _load()["rows_by_id"]


def hmdb_name(hmdb_id) -> str:
    return _load()["names"].get(hmdb_id, "")


def associations(hmdb_id) -> pd.DataFrame:
    """All HMDB protein associations of one metabolite: columns Gene, Protein, UniProt (HMDB order)."""
    s = _load()
    cache = s.setdefault("assoc_cache", {})
    if hmdb_id not in cache:
        idx = s["rows_by_id"].get(hmdb_id)
        cache[hmdb_id] = (pd.DataFrame(columns=["Gene", "Protein", "UniProt"]) if idx is None else
                          s["table"].iloc[idx][["Gene", "Protein", "UniProt"]].astype(str).reset_index(drop=True))
    return cache[hmdb_id].copy()


def genes(hmdb_id, symbols_only: bool = False) -> list:
    """Distinct gene labels of one metabolite (optionally only HGNC-style symbols)."""
    g = list(dict.fromkeys(associations(hmdb_id)["Gene"]))
    return [x for x in g if is_gene_symbol(x)] if symbols_only else g


def lookup_by_name(name):
    """HMDB accession for an exact (normalized) HMDB common name; None if absent or ambiguous."""
    if name is None:
        return None
    hits = _load()["name_index"].get(normalize_name(name), [])
    return hits[0] if len(hits) == 1 else None


def symbol_for_collection(symbol: str, collection_genes) -> str:
    """
    The form of `symbol` to use when matching against a gene-set collection: the HMDB symbol
    itself if the collection contains it, otherwise its current (or retired) HGNC equivalent if
    the collection contains that, otherwise the HMDB symbol unchanged. `collection_genes` is a
    set of UPPER-CASE symbols (matching is case-insensitive, as in the ORA engine).
    """
    s = str(symbol)
    u = s.upper()
    if collection_genes is None or u in collection_genes:
        return s
    cur = RETIRED_TO_CURRENT.get(s)
    if cur and cur[0].upper() in collection_genes:
        return cur[0]
    for old in CURRENT_TO_RETIRED.get(s, ()):
        if old.upper() in collection_genes:
            return old
    return s
