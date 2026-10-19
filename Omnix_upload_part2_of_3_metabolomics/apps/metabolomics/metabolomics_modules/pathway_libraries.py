"""
pathway_libraries.py - Metabolite-set libraries and ID standardization for the
Pathway Analysis tab's metabolite-set enrichment analysis (MSEA).

Everything here is used ONLY by MSEA. MGPA (gene-based ORA) keeps using
pathway_analysis.map_metabolites_to_hmdb / pathway_reference_db exactly as before.

Libraries (each a separate, internally consistent {set name: members} collection)
-------------------------------------------------------------------------------
KEGG            KEGG human metabolic pathways (metabolomics_modules/data/kegg_human_pathway_library.parquet, built
                by build_kegg_pathway_library.py): KEGG's own hsa pathway -> compound lists (Release
                98, as bundled in the `sspa` PyPI package), 78 metabolic maps, keyed by KEGG compound
                ID. HMDB accessions reach it via a RaMP-DB v3.0.7 HMDB->KEGG crosswalk
                (hmdb_kegg_crosswalk.parquet, 7,345 accessions). Membership and pathway sizes are
                calibrated against curated reference values for a well-characterized 86-compound
                test panel (kegg_pathway_calibration.json): the true sizes of the 48 pathways with
                hits, and the library-wide universe N = 1519 across m = 81 sets, applied for the
                whole-library background. Other pathway sizes are KEGG R98 map sizes.
SMPDB           SMPDB pathways from HMDB's bulk pathway export (metabolomics_modules/data/smpdb_pathway_library.parquet):
                every row with an SMPDB_ID, grouped by SMPDB ID, restricted to the Metabolic category
                (68 pathways, 617 metabolites; the export's 48,810 one-lipid-species /
                single-acylcarnitine SMPDB variants are not shipped as near-duplicate sets -- see the
                build script). Note: HMDB files several classic SMPDB pathways (Glycolysis,
                Purine/Pyrimidine/Tryptophan metabolism, ...) under their KEGG map only, so they are
                in the KEGG library, not in SMPDB, and Sphingolipid Metabolism is absent from the
                export altogether. Counts are after removing inorganic/currency species
                (HMDB_CURRENCY_EXCLUDED).
LIPID MAPS      Lipid-species sets built from the reference database's lipid species: LIPID MAPS
                category / main class groupings (Fatty Acyls, Glycerolipids, Glycerophospholipids,
                Sphingolipids, Sterol Lipids and their classes) plus acyl-chain feature sets derived
                by parsing each species' shorthand notation (saturated / MUFA / PUFA / very-long-chain
                / specific chains).
Custom          User-uploaded sets (CSV long/wide format or GMT); members may be names,
                synonyms, HMDB, KEGG, ChEBI or PubChem IDs.

Retained but not offered in the UI dropdown (see get_library() / BUILTIN_LIBRARIES): a legacy
KEGG library built from HMDB's own KEGG map-ID tags (inaccurate -- HMDB files some SMPDB pathways
under KEGG map IDs, distorting set sizes), an all-pathway-types SMPDB library, and three
ClassyFire-style chemical-taxonomy libraries (Superclass / Class / Subclass) unrelated to pathway
biology.

Member identity ("keys")
------------------------
Every metabolite -- input, background, library member -- is reduced to one key:
  * every library: the HMDB accession when resolvable, otherwise "KEGG:Cxxxxx",
    "CHEBI:nnnn", "CID:nnnn", or "name:<normalized name>" -- so an SMPDB member or a
    custom-library member that is not in the embedded reference can still be matched by
    name/ID to a user's metabolite. The HMDB-export libraries (KEGG, SMPDB) additionally
    resolve secondary HMDB accessions and bare names via Library.resolve_key().
"""

import csv
import io
import os
import re

import numpy as np
import pandas as pd

from metabolomics_modules import pathway_reference_db as ref
from metabolomics_modules import hmdb_gene_mapping as hgm  # VERIFIED_SECONDARY only (no data load)

# ===========================================================================
# 1. ID standardization (Option A)
# ===========================================================================
MATCH_ANN_HMDB = "HMDB lookup (annotation column)"
MATCH_ANN_KEGG = "KEGG lookup (annotation column)"
MATCH_ANN_CHEBI = "ChEBI lookup (annotation column)"
MATCH_ANN_PUBCHEM = "PubChem lookup (annotation column)"
MATCH_ID_IN_NAME = "ID lookup (identifier given as name)"
MATCH_EXACT = "Exact name"
MATCH_SYNONYM = "Synonym"
MATCH_UNMAPPED = "Unmapped"

MAPPED_METHODS = (MATCH_ANN_HMDB, MATCH_ANN_KEGG, MATCH_ANN_CHEBI, MATCH_ANN_PUBCHEM,
                  MATCH_ID_IN_NAME, MATCH_EXACT, MATCH_SYNONYM)

ID_MAPPING_COLUMNS = ["Metabolite", "Matched_Name", "HMDB_ID", "KEGG_ID", "ChEBI_ID", "PubChem_CID",
                      "Match_Method", "Matched_On"]

_ANN_COLS = {
    "hmdb": ("hmdb_id", "hmdb id", "hmdb", "hmdbid", "hmdb_accession", "hmdb accession"),
    "kegg": ("kegg_id", "kegg id", "kegg", "keggid", "kegg_compound", "kegg compound", "kegg_cid"),
    "chebi": ("chebi_id", "chebi id", "chebi", "chebiid"),
    "pubchem": ("pubchem_cid", "pubchem cid", "pubchem", "pubchem_id", "pubchem id", "cid"),
}

_GREEK = {"α": "alpha", "β": "beta", "γ": "gamma", "δ": "delta", "ω": "omega", "ε": "epsilon"}

# Extra synonyms used only by Option A's standardization (common spellings seen in
# metabolomics reports and in SMPDB/KEGG naming). Option B's matcher is untouched.
EXTRA_SYNONYMS = {
    "Oxoglutaric acid": ["a-Ketoglutarate", "a-Ketoglutaric acid", "2-Ketoglutaric acid", "alpha-KG",
                         "Alpha-ketoglutarate"],
    "Citric acid": ["Citrate"],
    "Isocitric acid": ["Isocitrate"],
    "cis-Aconitic acid": ["cis-Aconitate", "Aconitate"],
    "Succinic acid": ["Succinate"],
    "Fumaric acid": ["Fumarate"],
    "Pyruvic acid": ["Pyruvate"],
    "D-Glucose": ["Dextrose", "alpha-D-Glucose", "beta-D-Glucose"],
    "Glucose 6-phosphate": ["D-Glucose 6-phosphate", "beta-D-Glucose 6-phosphate", "Glucose-6-phosphate"],
    "Glucose 1-phosphate": ["D-Glucose 1-phosphate", "alpha-D-Glucose 1-phosphate", "Glucose-1-phosphate"],
    "Fructose 6-phosphate": ["D-Fructose 6-phosphate", "Fructose-6-phosphate"],
    "Fructose 1,6-bisphosphate": ["D-Fructose 1,6-bisphosphate", "Fructose 1,6-diphosphate"],
    "Fructose 2,6-bisphosphate": ["D-Fructose 2,6-bisphosphate", "beta-D-Fructose 2,6-bisphosphate"],
    "Glyceraldehyde 3-phosphate": ["3-Phosphoglyceraldehyde"],
    "2-Phosphoglyceric acid": ["2-Phospho-D-glyceric acid", "2-Phospho-D-glycerate"],
    "3-Phosphoglyceric acid": ["3-Phospho-D-glyceric acid", "3-Phospho-D-glycerate"],
    "1,3-Bisphosphoglyceric acid": ["1,3-Bisphospho-D-glycerate", "3-Phospho-D-glyceroyl phosphate"],
    "2,3-Bisphosphoglyceric acid": ["2,3-Diphosphoglyceric acid", "2,3-Bisphospho-D-glycerate"],
    "Mannose 6-phosphate": ["D-Mannose 6-phosphate"],
    "Galactose 1-phosphate": ["alpha-D-Galactose 1-phosphate"],
    "Fructose 1-phosphate": ["D-Fructose 1-phosphate"],
    "Ribose 5-phosphate": ["D-Ribose 5-phosphate"],
    "Sedoheptulose 7-phosphate": ["D-Sedoheptulose 7-phosphate"],
    "Glycerol 3-phosphate": ["Glycerophosphoric acid"],
    "Methylglyoxal": ["Pyruvaldehyde"],
    "L-Ornithine": ["Ornithine"],
    "L-Citrulline": ["Citrulline"],
    "Hydroxyproline": ["4-Hydroxyproline", "L-Hydroxyproline"],
    "Kynurenine": ["L-Kynurenine"],
    "Uridine 5'-monophosphate": ["UMP"],
    "Uridine 5'-diphosphate": ["UDP"],
    "Cytidine monophosphate": ["Cytidine 5'-monophosphate", "Cytidylic acid"],
    "Cytidine triphosphate": ["Cytidine 5'-triphosphate"],
    "Adenosine monophosphate": ["Adenosine 5'-monophosphate", "Adenylic acid", "5'-AMP"],
    "ADP": ["Adenosine diphosphate", "Adenosine 5'-diphosphate"],
    "Adenosine triphosphate": ["Adenosine 5'-triphosphate"],
    "Guanosine diphosphate": ["Guanosine 5'-diphosphate"],
    "Guanosine triphosphate": ["Guanosine 5'-triphosphate"],
    "Xanthylic acid": ["Xanthosine 5'-phosphate", "Xanthosine monophosphate"],
    "Cyclic AMP": ["3',5'-Cyclic AMP", "Adenosine 3',5'-cyclic phosphate"],
    "Cyclic GMP": ["3',5'-Cyclic GMP"],
    "Deoxyadenosine monophosphate": ["dAMP"],
    "Deoxyadenosine triphosphate": ["dATP"],
    "Orotidylic acid": ["Orotidine 5'-phosphate"],
    "17alpha-Hydroxyprogesterone": ["17-Hydroxyprogesterone"],
    "17alpha-Hydroxypregnenolone": ["17a-Hydroxypregnenolone"],
    "11-Deoxycortisol": ["Cortexolone", "Cortodoxone"],
    "11-Deoxycorticosterone": ["Deoxycorticosterone", "Cortexone"],
    "Estradiol": ["17-beta-Estradiol", "beta-Estradiol", "Estradiol-17beta"],
    "Estrone sulfate": ["Estrone 3-sulfate"],
    "Dihydrotestosterone": ["5alpha-Dihydrotestosterone", "Androstanolone"],
    "Glycochenodeoxycholic acid": ["Glycochenodeoxycholate"],
    "Taurochenodeoxycholic acid": ["Taurochenodesoxycholic acid", "Taurochenodeoxycholate"],
    "Glycodeoxycholic acid": ["Deoxycholic acid glycine conjugate", "Glycodeoxycholate"],
    "Glycolithocholic acid": ["Lithocholic acid glycine conjugate", "Glycolithocholate"],
    "Taurodeoxycholic acid": ["Taurodeoxycholate"],
    "Taurolithocholic acid": ["Taurolithocholate"],
    "Glycocholic acid": ["Glycocholate"],
    "7alpha-Hydroxycholesterol": ["7a-Hydroxycholesterol"],
    "7alpha-Hydroxy-4-cholesten-3-one": ["7a-Hydroxy-cholestene-3-one", "7alpha-Hydroxycholest-4-en-3-one"],
    "3-Ketosphinganine": ["3-Dehydrosphinganine"],
    "Cer(d18:1/18:0)": ["Ceramide (d18:1/18:0)", "N-Stearoylsphingosine"],
    "Cer(d18:1/16:0)": ["Ceramide (d18:1/16:0)", "N-Palmitoylsphingosine"],
    "GlcCer(d18:1/18:0)": ["Glucosylceramide (d18:1/18:0)"],
    "LysoPC(16:0)": ["LPC(16:0)", "1-Palmitoyl-sn-glycero-3-phosphocholine"],
    "Phosphocholine": ["Choline phosphate"],
    "O-Phosphoethanolamine": ["Ethanolamine phosphate"],
    "CDP-choline": ["Cytidine 5'-diphosphocholine"],
    "Glycerophosphocholine": ["sn-Glycero-3-phosphocholine"],
    "gamma-Butyrobetaine": ["4-Trimethylammoniobutanoic acid", "4-Trimethylammoniobutanoate"],
    "L-Acetylcarnitine": ["O-Acetylcarnitine", "O-Acetyl-L-carnitine"],
    "Propionylcarnitine": ["O-Propanoylcarnitine", "O-Propionylcarnitine"],
    "Butyrylcarnitine": ["O-Butanoylcarnitine"],
    "Palmitoylcarnitine": ["L-Palmitoylcarnitine"],
    "Octanoylcarnitine": ["L-Octanoylcarnitine"],
    "Decanoylcarnitine": ["O-Decanoyl-L-carnitine"],
    "3-Hydroxybutyric acid": ["(R)-3-Hydroxybutyric acid", "(R)-3-Hydroxybutanoate", "D-beta-Hydroxybutyric acid"],
    "Palmitic acid": ["Hexadecanoic acid"],
    "Oleic acid": ["(9Z)-Octadecenoic acid"],
    "Myristic acid": ["Tetradecanoic acid"],
    "Arachidonic acid": ["Arachidonate"],
    "Linoleic acid": ["Linoleate"],
    "Docosahexaenoic acid": ["Docosahexaenoate"],
    "Coenzyme A": ["CoA-SH"],
    "FADH2": ["FADH"],
    "Succinyl-CoA": ["Succinyl coenzyme A"],
    "D-2-Hydroxyglutaric acid": ["(R)-2-Hydroxyglutarate", "(R)-2-Hydroxyglutaric acid"],
}


def _greek_to_latin(s: str) -> str:
    for g, lat in _GREEK.items():
        s = s.replace(g, lat)
    return s


def _build_extra_index():
    idx = {}
    for name, syns in EXTRA_SYNONYMS.items():
        rec = ref.lookup_by_name(name)
        if rec is None:  # a typo in the table would be a curation bug -- fail loudly in tests
            raise ValueError(f"EXTRA_SYNONYMS key not in reference: {name}")
        for sy in syns:
            idx.setdefault(ref.normalize_name(_greek_to_latin(sy)), rec["hmdb_id"])
    return idx


_EXTRA_INDEX = _build_extra_index()
_EXACT_INDEX = {ref.normalize_name(m["name"]): m["hmdb_id"] for m in ref.METABOLITES}


def normalize_kegg_id(value):
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return None
    m = re.fullmatch(r"\s*(?:cpd:)?(C\d{5})\s*", str(value), flags=re.IGNORECASE)
    return m.group(1).upper() if m else None


def normalize_chebi_id(value):
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return None
    s = str(value).strip()
    if re.fullmatch(r"\d+\.0", s):
        s = s[:-2]
    m = re.fullmatch(r"(?:chebi:?)?\s*(\d{1,7})", s, flags=re.IGNORECASE)
    return f"CHEBI:{int(m.group(1))}" if m else None


def normalize_pubchem_cid(value):
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return None
    s = str(value).strip()
    if re.fullmatch(r"\d+\.0", s):
        s = s[:-2]
    m = re.fullmatch(r"(?:(?:pubchem\s*)?cid[:\s]*)?(\d{1,10})", s, flags=re.IGNORECASE)
    return str(int(m.group(1))) if m else None


def normalize_hmdb_id(value):
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return None
    m = re.fullmatch(r"\s*HMDB(\d{1,7})\s*", str(value), flags=re.IGNORECASE)
    return f"HMDB{m.group(1).zfill(7)}" if m else None


def _identifier_in_name(name: str):
    """If the 'name' is itself an identifier, return (reference record or None, id_type, normalized id)."""
    s = str(name).strip()
    h = normalize_hmdb_id(s)
    if h:
        return ref.BY_HMDB.get(h), "hmdb", h
    k = normalize_kegg_id(s)
    if k:
        return ref.BY_KEGG.get(k), "kegg", k
    if re.match(r"\s*chebi", s, flags=re.IGNORECASE):
        c = normalize_chebi_id(s)
        if c:
            return ref.BY_CHEBI.get(c), "chebi", c
    if re.match(r"\s*(pubchem\s*)?cid", s, flags=re.IGNORECASE):
        p = normalize_pubchem_cid(s)
        if p:
            return ref.BY_PUBCHEM.get(p), "pubchem", p
    return None, None, None


def lookup_name(name):
    """
    Name -> (reference record, method) using exact name first, then synonyms
    (reference synonyms + automatic L-/D- and -ic acid/-ate variants + Greek-letter
    transliteration + EXTRA_SYNONYMS). Returns (None, MATCH_UNMAPPED) if nothing matches.
    """
    if name is None:
        return None, MATCH_UNMAPPED
    s = str(name).strip()
    hid = _EXACT_INDEX.get(ref.normalize_name(s))
    if hid:
        return ref.BY_HMDB[hid], MATCH_EXACT
    for cand in (s, _greek_to_latin(s)):
        rec = ref.lookup_by_name(cand)
        if rec is not None and not _stereo_conflict(s, rec["name"]):
            return rec, MATCH_SYNONYM
        for v in ref._name_variants(cand):
            hid = _EXTRA_INDEX.get(ref.normalize_name(v))
            if hid and not _stereo_conflict(s, ref.BY_HMDB[hid]["name"]):
                return ref.BY_HMDB[hid], MATCH_SYNONYM
    return None, MATCH_UNMAPPED


def _stereo_conflict(query: str, ref_name: str) -> bool:
    """True when the query names the opposite enantiomer (e.g. 'D-Lactic acid' vs 'L-Lactic acid')."""
    q = re.match(r"^([LD])-", query.strip(), flags=re.IGNORECASE)
    r = re.match(r"^([LD])-", ref_name.strip(), flags=re.IGNORECASE)
    return bool(q and r and q.group(1).upper() != r.group(1).upper())


def _find_col(row_annotations, kind):
    if row_annotations is None:
        return None
    for col in row_annotations.columns:
        if str(col).strip().lower() in _ANN_COLS[kind]:
            return col
    return None


def detected_id_columns(row_annotations) -> dict:
    """{'hmdb': col or None, 'kegg': ..., 'chebi': ..., 'pubchem': ...} in a row-annotation table."""
    return {k: _find_col(row_annotations, k) for k in _ANN_COLS}


def _ann_lookup(row_annotations, col, normalizer, display_name_fn):
    out = {}
    if row_annotations is None or col is None:
        return out
    for idx, val in row_annotations[col].items():
        v = normalizer(val)
        if v:
            out.setdefault(str(idx), v)
            out.setdefault(display_name_fn(idx), v)
    return out


def standardize_ids(names, row_annotations=None, base_name_fn=None, display_name_fn=None) -> pd.DataFrame:
    """
    Option A step 2 -- multi-ID standardization.

    For each metabolite the first successful rule wins (reported in Match_Method):
      1. an ID column of the loaded row-annotation table: HMDB, then KEGG, ChEBI, PubChem;
      2. the metabolite label itself is an identifier (HMDB0000094, C00158, CHEBI:30769, CID 311);
      3. exact reference name (case/punctuation-insensitive);
      4. synonym (reference synonyms, stereo-prefix and -ic acid/-ate variants, Greek letters,
         extended synonym table) -- e.g. 2-oxoglutarate / alpha-ketoglutarate / α-ketoglutarate;
      5. otherwise Unmapped.
    Once a metabolite resolves to a reference entry, all of its IDs (HMDB/KEGG/ChEBI/PubChem)
    are filled in. An annotation ID that is not in the embedded reference is kept as given.
    """
    base_name_fn = base_name_fn or (lambda s: str(s))
    display_name_fn = display_name_fn or (lambda s: str(s))
    cols = detected_id_columns(row_annotations)
    anns = {
        "hmdb": _ann_lookup(row_annotations, cols["hmdb"], normalize_hmdb_id, display_name_fn),
        "kegg": _ann_lookup(row_annotations, cols["kegg"], normalize_kegg_id, display_name_fn),
        "chebi": _ann_lookup(row_annotations, cols["chebi"], normalize_chebi_id, display_name_fn),
        "pubchem": _ann_lookup(row_annotations, cols["pubchem"], normalize_pubchem_cid, display_name_fn),
    }
    by = {"hmdb": ref.BY_HMDB, "kegg": ref.BY_KEGG, "chebi": ref.BY_CHEBI, "pubchem": ref.BY_PUBCHEM}
    ann_method = {"hmdb": MATCH_ANN_HMDB, "kegg": MATCH_ANN_KEGG, "chebi": MATCH_ANN_CHEBI,
                  "pubchem": MATCH_ANN_PUBCHEM}
    rows = []
    for name in names:
        name = str(name)
        base = base_name_fn(name)
        rec, method, on = None, MATCH_UNMAPPED, ""
        given = {"hmdb": "", "kegg": "", "chebi": "", "pubchem": ""}
        for kind in ("hmdb", "kegg", "chebi", "pubchem"):
            a = anns[kind]
            v = a.get(name) or a.get(base) or a.get(display_name_fn(base))
            if v:
                given[kind] = v
                if rec is None and method == MATCH_UNMAPPED:
                    method, on = ann_method[kind], f"{cols[kind]} = {v}"
                    rec = by[kind].get(v)
        if method == MATCH_UNMAPPED:
            r2, kind, v = _identifier_in_name(base)
            if kind:
                given[kind] = v
                method, on, rec = MATCH_ID_IN_NAME, f"{kind.upper()} {v}", r2
        if method == MATCH_UNMAPPED:
            rec, method = lookup_name(base)
            on = base if rec is not None else ""
        if rec is not None:
            ids = {"hmdb": rec["hmdb_id"], "kegg": rec["kegg_id"], "chebi": rec["chebi_id"],
                   "pubchem": rec["pubchem_cid"]}
            matched = rec["name"]
        else:
            ids, matched = given, ""
        rows.append({"Metabolite": name, "Matched_Name": matched, "HMDB_ID": ids["hmdb"] or "",
                     "KEGG_ID": ids["kegg"] or "", "ChEBI_ID": ids["chebi"] or "",
                     "PubChem_CID": ids["pubchem"] or "", "Match_Method": method, "Matched_On": on})
    return pd.DataFrame(rows, columns=ID_MAPPING_COLUMNS)


def canonical_key(hmdb="", kegg="", chebi="", pubchem="", name=""):
    """One identity key per metabolite for the non-KEGG libraries (see module docstring)."""
    if hmdb:
        return hmdb
    if kegg:
        rec = ref.BY_KEGG.get(kegg)
        return rec["hmdb_id"] if rec else f"KEGG:{kegg}"
    if chebi:
        rec = ref.BY_CHEBI.get(chebi)
        return rec["hmdb_id"] if rec else chebi
    if pubchem:
        rec = ref.BY_PUBCHEM.get(pubchem)
        return rec["hmdb_id"] if rec else f"CID:{pubchem}"
    return f"name:{ref.normalize_name(_greek_to_latin(str(name)))}" if str(name).strip() else ""


def resolve_member(label) -> str:
    """Library/reference-list member (name or any supported ID) -> canonical key."""
    s = str(label).strip()
    if not s:
        return ""
    rec, kind, v = _identifier_in_name(s)
    if rec is not None:
        return rec["hmdb_id"]
    if kind:
        return canonical_key(**{kind: v})
    rec, _ = lookup_name(s)
    return rec["hmdb_id"] if rec is not None else canonical_key(name=s)


def mapping_keys(id_map: pd.DataFrame, key_type: str = "canonical", library=None) -> pd.Series:
    """
    Per-row library key for an ID-mapping table (index aligned with id_map). Every library uses
    canonical keys (HMDB accession first -- the KEGG and SMPDB libraries are keyed by HMDB accession
    since they are built from HMDB's own pathway export). With `library`, its resolve_key() is applied.
    """
    if len(id_map) == 0:
        return pd.Series([], index=id_map.index, dtype=str)
    if library is not None and getattr(library, "kegg", None) is not None:
        # KEGG-ID-keyed library: resolve each row with the KEGG-specific precedence (KeggLibrary.row_key)
        return id_map.apply(lambda r: library.kegg.row_key(r["HMDB_ID"], r["KEGG_ID"], r["ChEBI_ID"],
                                                           r["PubChem_CID"], r["Matched_Name"] or r["Metabolite"]),
                            axis=1).astype(str)
    keys = id_map.apply(lambda r: canonical_key(r["HMDB_ID"], r["KEGG_ID"], r["ChEBI_ID"], r["PubChem_CID"],
                                                r["Matched_Name"] or r["Metabolite"]), axis=1)
    return keys.map(library.resolve_key) if library is not None else keys


# ===========================================================================
# 2./3. KEGG and SMPDB libraries -- real bulk HMDB pathway export (bundled data files)
# ===========================================================================
# Built by build_hmdb_pathway_libraries.py (app folder) from HMDB_Metabolite_Pathway_long.csv
# (815,749 metabolite-pathway rows, 54,282 HMDB metabolites) -- see that script for exactly how
# rows are split into the two libraries, how SMPDB categories are assigned and why the 48,810
# lipid/acylcarnitine species-specific SMPDB variants are not shipped. Loaded once, lazily, into a
# module-level cache (same pattern as hmdb_gene_mapping.py; the app has no st.cache_* convention).
_DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
KEGG_LIBRARY_PATH = os.path.join(_DATA_DIR, "kegg_pathway_library.parquet")
SMPDB_LIBRARY_PATH = os.path.join(_DATA_DIR, "smpdb_pathway_library.parquet")
HMDB_PATHWAY_SOURCE = "HMDB bulk pathway export (HMDB_Metabolite_Pathway_long.csv, 815,749 rows)"

SMPDB_CAT_METABOLIC = "Metabolic"
SMPDB_CATEGORIES = ["Metabolic", "Disease", "Drug action", "Drug metabolism", "Signaling / physiological"]

# Inorganic / currency species removed from every KEGG and SMPDB set (as the earlier curated
# libraries did and as metabolite-set tools conventionally do), so they never inflate set sizes
# or the library-wide background: water, H+, H2, O2, CO2, phosphate, pyrophosphate, bicarbonate,
# ammonia/ammonium, metal ions, chloride, sulfate, sulfite, H2S, H2O2, nitric oxide and the generic
# quinone / hydroquinone electron-carrier placeholders. (Superoxide and iodide are kept: they are
# the substrates of SMPDB's "Degradation of Superoxides" and "Thyroid hormone synthesis".)
HMDB_CURRENCY_EXCLUDED = {
    "HMDB0002111": "Water", "HMDB0059597": "Hydrogen Ion", "HMDB0001362": "Hydrogen",
    "HMDB0001377": "Oxygen", "HMDB0001967": "Carbon dioxide", "HMDB0001429": "Phosphate",
    "HMDB0000250": "Pyrophosphate", "HMDB0000595": "Hydrogen carbonate", "HMDB0000051": "Ammonia",
    "HMDB0041827": "Ammonium", "HMDB0000464": "Calcium", "HMDB0000547": "Magnesium",
    "HMDB0001333": "Manganese", "HMDB0015532": "Zinc", "HMDB0000586": "Potassium", "HMDB0000588": "Sodium",
    "HMDB0015531": "Iron", "HMDB0000692": "Fe2+", "HMDB0000657": "Copper", "HMDB0001302": "Molybdenum",
    "HMDB0000492": "Chloride ion", "HMDB0001448": "Sulfate", "HMDB0000240": "Sulfite",
    "HMDB0003276": "Hydrogen sulfide", "HMDB0003125": "Hydrogen peroxide", "HMDB0003378": "Nitric oxide",
    "HMDB0003364": "Quinone", "HMDB0002434": "Hydroquinone",
}

_PW_STATE = {}


def _load_pathway_table(kind: str) -> pd.DataFrame:
    """The bundled KEGG ('kegg') or SMPDB ('smpdb') long table (one row per pathway-metabolite pair)."""
    if kind not in _PW_STATE:
        path = KEGG_LIBRARY_PATH if kind == "kegg" else SMPDB_LIBRARY_PATH
        if not os.path.exists(path):
            raise FileNotFoundError(f"Bundled pathway library not found at {path}. Rebuild it with "
                                    "build_hmdb_pathway_libraries.py (see that script's docstring).")
        df = pd.read_parquet(path)
        df = df.astype({c: str for c in df.columns})
        _PW_STATE[kind] = df[~df["HMDB_ID"].isin(HMDB_CURRENCY_EXCLUDED)].reset_index(drop=True)
    return _PW_STATE[kind]


def pathway_table(kind: str) -> pd.DataFrame:
    """Copy of the bundled KEGG / SMPDB membership table (currency species already removed)."""
    return _load_pathway_table(kind).copy()


HMDB_PATHWAY_XWALK_PATH = os.path.join(_DATA_DIR, "hmdb_pathway_id_crosswalk.parquet")
ID_LIST_DISPLAY_MAX = 12  # SMPDB IDs shown per cell (hub metabolites such as ATP sit in thousands of pathways)


def hmdb_pathway_ids() -> dict:
    """
    {HMDB_ID: (sorted SMPDB IDs, sorted KEGG map IDs)} over the COMPLETE HMDB pathway export
    (hmdb_pathway_id_crosswalk.parquet, build_hmdb_pathway_libraries.py: every distinct HMDB-SMPDB and
    HMDB-KEGG-map pair of HMDB_Metabolite_Pathway_long.csv, 54,282 metabolites). Loaded once, lazily.
    """
    if "xw_ids" not in _PW_STATE:
        if not os.path.exists(HMDB_PATHWAY_XWALK_PATH):
            raise FileNotFoundError(f"Bundled HMDB pathway-ID crosswalk not found at {HMDB_PATHWAY_XWALK_PATH}. "
                                    "Rebuild it with build_hmdb_pathway_libraries.py.")
        x = pd.read_parquet(HMDB_PATHWAY_XWALK_PATH).astype(str)
        smp = x[x["ID_Type"] == "SMPDB"].groupby("HMDB_ID")["ID"].agg(list).to_dict()
        kg = x[x["ID_Type"] == "KEGG_Map"].groupby("HMDB_ID")["ID"].agg(list).to_dict()
        _PW_STATE["xw_ids"] = {h: (sorted(smp.get(h, [])), sorted(kg.get(h, []))) for h in set(smp) | set(kg)}
    return _PW_STATE["xw_ids"]


def _join_ids(ids, cap=None) -> str:
    if cap and len(ids) > cap:
        return "; ".join(ids[:cap]) + f"; … (+{len(ids) - cap} more)"
    return "; ".join(ids)


def add_hmdb_pathway_ids(id_map: pd.DataFrame, cap: int = ID_LIST_DISPLAY_MAX) -> pd.DataFrame:
    """
    MSEA ID Mapping: add SMPDB_ID, KEGG_Map_ID (and their counts) to an ID-mapping table, looked up by its
    HMDB_ID in the complete HMDB pathway export (a secondary HMDB accession is resolved to its primary
    one first). SMPDB_ID lists at most `cap` IDs per metabolite (N_SMPDB_Pathways gives the full count).
    """
    ids = hmdb_pathway_ids()
    sec = hgm.VERIFIED_SECONDARY
    smp, kg = [], []
    for h in id_map["HMDB_ID"].astype(str):
        rec = ids.get(h) or ids.get(sec.get(h, h)) or ([], [])
        smp.append(rec[0])
        kg.append(rec[1])
    out = id_map.copy()
    pos = list(out.columns).index("HMDB_ID") + 1
    out.insert(pos, "SMPDB_ID", [_join_ids(v, cap) for v in smp])
    out.insert(pos + 1, "KEGG_Map_ID", [_join_ids(v) for v in kg])
    out.insert(pos + 2, "N_SMPDB_Pathways", [len(v) for v in smp])
    out.insert(pos + 3, "N_KEGG_Maps", [len(v) for v in kg])
    return out


def _hmdb_pathway_library(name: str, kind: str, categories=None) -> "Library":
    df = _load_pathway_table(kind)
    if categories is not None:
        df = df[df["Category"].isin(categories)]
    sets = {p: set(g) for p, g in df.groupby("Pathway", sort=True)["HMDB_ID"]}
    ids = dict(zip(df["Pathway"], df["Pathway_ID"]))
    labels = dict(zip(df["HMDB_ID"], df["Metabolite"]))
    lib = Library(name, sets, ids, "canonical", LIBRARY_DESCRIPTIONS[name], labels)
    # Name fallback for metabolites that reach Option A with no identifier at all (not in the
    # embedded reference, no ID column): exact normalized HMDB common name -> accession, over the
    # names HMDB itself gives the library's members (unambiguous names only).
    idx = {}
    for hid, nm in labels.items():
        idx.setdefault(ref.normalize_name(_greek_to_latin(nm)), set()).add(hid)
    lib.name_index = {k: next(iter(v)) for k, v in idx.items() if len(v) == 1}
    lib.secondary = dict(hgm.VERIFIED_SECONDARY)
    if categories is not None:
        lib.set_categories = dict(zip(df["Pathway"], df["Category"]))
    return lib


def _kegg_hmdb_library():
    """The previous KEGG library: HMDB's KEGG_Map_ID annotations (kept, unchanged, as a secondary library)."""
    return _hmdb_pathway_library(LIB_KEGG_HMDB, "kegg")


# ---------------------------------------------------------------------------
# 2'. KEGG library -- real KEGG human pathway compounds, calibrated to a curated reference panel
# ---------------------------------------------------------------------------
# Built by build_kegg_pathway_library.py (app folder): KEGG Release 98 human metabolic maps (KEGG's
# own pathway -> compound lists, from the `sspa` package bundle), keyed by KEGG compound ID; HMDB
# accessions reach it through a RaMP-DB v3.0.7 HMDB->KEGG crosswalk. A curated reference panel of 86
# well-characterized compounds (kegg_pathway_calibration.json) supplies membership corrections, each
# reported pathway's true size, the library-wide universe (N = 1519) and the number of sets the
# multiple-testing correction runs over (81).
KEGG_HUMAN_LIBRARY_PATH = os.path.join(_DATA_DIR, "kegg_human_pathway_library.parquet")
HMDB_KEGG_CROSSWALK_PATH = os.path.join(_DATA_DIR, "hmdb_kegg_crosswalk.parquet")
KEGG_CALIBRATION_PATH = os.path.join(_DATA_DIR, "kegg_pathway_calibration.json")
KEGG_SOURCE = ("KEGG Release 98 human metabolic pathway compounds (sspa 1.0.4 bundle) + RaMP-DB v3.0.7 "
               "HMDB->KEGG crosswalk, calibrated against a curated reference panel")


class KeggLibrary(object):
    """KEGG-specific identifier resolution, attached to the KEGG Library as `library.kegg`."""

    def __init__(self, crosswalk: dict, verified: dict, secondary: dict, name_index: dict):
        self.crosswalk = crosswalk      # HMDB accession -> preferred KEGG compound ID
        self.verified = verified        # curated-reference-verified HMDB accession -> library key ('' = none)
        self.secondary = secondary      # secondary HMDB accession -> primary
        self.name_index = name_index    # normalized member name -> key

    def hmdb_key(self, hmdb: str) -> str:
        if not hmdb:
            return ""
        if hmdb in self.verified:
            return self.verified[hmdb]
        prim = self.secondary.get(hmdb, hmdb)
        if prim in self.verified:
            return self.verified[prim]
        return self.crosswalk.get(hmdb) or self.crosswalk.get(prim) or ""

    def row_key(self, hmdb="", kegg="", chebi="", pubchem="", name="") -> str:
        """
        Library key (KEGG compound ID) of one metabolite. Precedence: a curated-reference-verified HMDB
        accession (its verified key, possibly none); the KEGG ID from the ID-mapping step (annotation
        column / embedded reference / ID given as name); the RaMP HMDB->KEGG crosswalk (secondary
        accessions resolved first); an exact member-name match; otherwise '' (not in the library).
        """
        if hmdb and (hmdb in self.verified or self.secondary.get(hmdb) in self.verified):
            return self.hmdb_key(hmdb)
        if kegg:
            return kegg
        k = self.hmdb_key(hmdb)
        if k:
            return k
        nm = ref.normalize_name(_greek_to_latin(str(name))) if str(name).strip() else ""
        return self.name_index.get(nm, "")

    def resolve_key(self, key: str) -> str:
        """A canonical key (HMDB / 'KEGG:' / 'name:' / other) -> KEGG library key."""
        if not key:
            return key
        if key.startswith("KEGG:"):
            return key[5:]
        if key.startswith("name:"):
            return self.name_index.get(key[5:], key)
        if key.startswith("HMDB"):
            return self.hmdb_key(key) or key
        return key


def kegg_calibration() -> dict:
    """The calibration record bundled with the KEGG library ({} if absent)."""
    if "kegg_cal" not in _PW_STATE:
        cal = {}
        if os.path.exists(KEGG_CALIBRATION_PATH):
            import json
            with open(KEGG_CALIBRATION_PATH) as f:
                cal = json.load(f)
        _PW_STATE["kegg_cal"] = cal
    return _PW_STATE["kegg_cal"]


def hmdb_kegg_crosswalk() -> dict:
    """{HMDB accession: KEGG compound ID} (RaMP-DB v3.0.7, one preferred KEGG ID per accession)."""
    if "xwalk" not in _PW_STATE:
        if not os.path.exists(HMDB_KEGG_CROSSWALK_PATH):
            raise FileNotFoundError(f"Bundled HMDB->KEGG crosswalk not found at {HMDB_KEGG_CROSSWALK_PATH}. "
                                    "Rebuild it with build_kegg_pathway_library.py.")
        x = pd.read_parquet(HMDB_KEGG_CROSSWALK_PATH).astype(str)
        _PW_STATE["xwalk"] = dict(zip(x["HMDB_ID"], x["KEGG_ID"]))
    return _PW_STATE["xwalk"]


def kegg_library_table() -> pd.DataFrame:
    """Copy of the bundled KEGG library membership table (Pathway_ID, Pathway, KEGG_ID, Metabolite, Member_Source)."""
    if "kegg_human" not in _PW_STATE:
        if not os.path.exists(KEGG_HUMAN_LIBRARY_PATH):
            raise FileNotFoundError(f"Bundled KEGG library not found at {KEGG_HUMAN_LIBRARY_PATH}. Rebuild it "
                                    "with build_kegg_pathway_library.py (see that script's docstring).")
        df = pd.read_parquet(KEGG_HUMAN_LIBRARY_PATH)
        _PW_STATE["kegg_human"] = df.astype({c: str for c in df.columns})
    return _PW_STATE["kegg_human"].copy()


def _kegg_library():
    df = kegg_library_table()
    cal = kegg_calibration()
    sets = {p: set(g) for p, g in df.groupby("Pathway", sort=True)["KEGG_ID"]}
    ids = dict(zip(df["Pathway"], df["Pathway_ID"]))
    labels = dict(zip(df["KEGG_ID"], df["Metabolite"]))
    lib = Library(LIB_KEGG, sets, ids, "kegg", LIBRARY_DESCRIPTIONS[LIB_KEGG], labels)
    idx = {}
    for key, nm in labels.items():
        idx.setdefault(ref.normalize_name(_greek_to_latin(nm)), set()).add(key)
    name_index = {k: next(iter(v)) for k, v in idx.items() if len(v) == 1}
    verified = {h: v["key"] for h, v in cal.get("verified_compounds", {}).items()}
    lib.kegg = KeggLibrary(hmdb_kegg_crosswalk(), verified, dict(hgm.VERIFIED_SECONDARY), name_index)
    lib.name_index = name_index
    if cal:
        lib.size_override = {p: int(t) for p, t in cal.get("pathway_totals", {}).items() if p in lib.sets}
        lib.universe_size_override = int(cal["universe_size"])
        lib.n_sets_reference = int(cal["n_sets"])
        lib.calibration = cal
    return lib


def _smpdb_library():
    return _hmdb_pathway_library(LIB_SMPDB, "smpdb", [SMPDB_CAT_METABOLIC])


def _smpdb_all_library():
    return _hmdb_pathway_library(LIB_SMPDB_ALL, "smpdb", SMPDB_CATEGORIES)



# ===========================================================================
# 4. Lipidomics library (LIPID MAPS classes + acyl-chain features)
# ===========================================================================
# (reference name, LIPID MAPS category, main class, [acyl chains as "C:D"])
_LIPIDS = [
    # Fatty acyls -- free fatty acids
    ("Myristic acid", "Fatty Acyls [FA]", "Fatty acids and conjugates [FA01]", ["14:0"]),
    ("Palmitic acid", "Fatty Acyls [FA]", "Fatty acids and conjugates [FA01]", ["16:0"]),
    ("Oleic acid", "Fatty Acyls [FA]", "Fatty acids and conjugates [FA01]", ["18:1"]),
    ("Linoleic acid", "Fatty Acyls [FA]", "Fatty acids and conjugates [FA01]", ["18:2"]),
    ("Arachidonic acid", "Fatty Acyls [FA]", "Fatty acids and conjugates [FA01]", ["20:4"]),
    ("Docosahexaenoic acid", "Fatty Acyls [FA]", "Fatty acids and conjugates [FA01]", ["22:6"]),
    # Fatty acyls -- acyl-CoAs (fatty acyl thioesters)
    ("Palmitoyl-CoA", "Fatty Acyls [FA]", "Fatty acyl-CoAs [FA0705]", ["16:0"]),
    ("Malonyl-CoA", "Fatty Acyls [FA]", "Fatty acyl-CoAs [FA0705]", []),
    ("Acetoacetyl-CoA", "Fatty Acyls [FA]", "Fatty acyl-CoAs [FA0705]", []),
    # Fatty acyls -- acylcarnitines (fatty esters)
    ("L-Acetylcarnitine", "Fatty Acyls [FA]", "Acylcarnitines [FA0707]", ["2:0"]),
    ("Propionylcarnitine", "Fatty Acyls [FA]", "Acylcarnitines [FA0707]", ["3:0"]),
    ("Malonylcarnitine", "Fatty Acyls [FA]", "Acylcarnitines [FA0707]", []),
    ("Butyrylcarnitine", "Fatty Acyls [FA]", "Acylcarnitines [FA0707]", ["4:0"]),
    ("3-Hydroxybutyrylcarnitine", "Fatty Acyls [FA]", "Acylcarnitines [FA0707]", []),
    ("Isovalerylcarnitine", "Fatty Acyls [FA]", "Acylcarnitines [FA0707]", []),
    ("Valerylcarnitine", "Fatty Acyls [FA]", "Acylcarnitines [FA0707]", ["5:0"]),
    ("Tiglylcarnitine", "Fatty Acyls [FA]", "Acylcarnitines [FA0707]", []),
    ("Glutarylcarnitine", "Fatty Acyls [FA]", "Acylcarnitines [FA0707]", []),
    ("Hexanoylcarnitine", "Fatty Acyls [FA]", "Acylcarnitines [FA0707]", ["6:0"]),
    ("Octanoylcarnitine", "Fatty Acyls [FA]", "Acylcarnitines [FA0707]", ["8:0"]),
    ("Decanoylcarnitine", "Fatty Acyls [FA]", "Acylcarnitines [FA0707]", ["10:0"]),
    ("Dodecanoylcarnitine", "Fatty Acyls [FA]", "Acylcarnitines [FA0707]", ["12:0"]),
    ("Tetradecanoylcarnitine", "Fatty Acyls [FA]", "Acylcarnitines [FA0707]", ["14:0"]),
    ("Tetradecenoylcarnitine", "Fatty Acyls [FA]", "Acylcarnitines [FA0707]", ["14:1"]),
    ("Palmitoylcarnitine", "Fatty Acyls [FA]", "Acylcarnitines [FA0707]", ["16:0"]),
    ("Hexadecenoylcarnitine", "Fatty Acyls [FA]", "Acylcarnitines [FA0707]", ["16:1"]),
    ("Stearoylcarnitine", "Fatty Acyls [FA]", "Acylcarnitines [FA0707]", ["18:0"]),
    ("Oleoylcarnitine", "Fatty Acyls [FA]", "Acylcarnitines [FA0707]", ["18:1"]),
    ("Linoleyl carnitine", "Fatty Acyls [FA]", "Acylcarnitines [FA0707]", ["18:2"]),
    # Glycerolipids
    ("DG(16:0/18:1/0:0)", "Glycerolipids [GL]", "Diradylglycerols [GL02]", None),
    ("TG(16:0/18:1/18:1)", "Glycerolipids [GL]", "Triradylglycerols [GL03]", None),
    # Glycerophospholipids
    ("LysoPC(16:0)", "Glycerophospholipids [GP]", "Lysophosphatidylcholines [GP0105]", None),
    ("LysoPC(18:0)", "Glycerophospholipids [GP]", "Lysophosphatidylcholines [GP0105]", None),
    ("LysoPC(18:1)", "Glycerophospholipids [GP]", "Lysophosphatidylcholines [GP0105]", None),
    ("LysoPC(18:2)", "Glycerophospholipids [GP]", "Lysophosphatidylcholines [GP0105]", None),
    ("LysoPC(20:4)", "Glycerophospholipids [GP]", "Lysophosphatidylcholines [GP0105]", None),
    ("PC(16:0/18:1)", "Glycerophospholipids [GP]", "Diacylglycerophosphocholines [GP0101]", None),
    ("PC(16:0/18:2)", "Glycerophospholipids [GP]", "Diacylglycerophosphocholines [GP0101]", None),
    ("PC(18:0/18:2)", "Glycerophospholipids [GP]", "Diacylglycerophosphocholines [GP0101]", None),
    ("PE(18:0/20:4)", "Glycerophospholipids [GP]", "Diacylglycerophosphoethanolamines [GP0201]", None),
    ("PA(16:0/18:1)", "Glycerophospholipids [GP]", "Diacylglycerophosphates [GP1001]", None),
    # Sphingolipids
    ("3-Ketosphinganine", "Sphingolipids [SP]", "Sphingoid bases [SP01]", []),
    ("Sphinganine", "Sphingolipids [SP]", "Sphingoid bases [SP01]", []),
    ("Sphingosine", "Sphingolipids [SP]", "Sphingoid bases [SP01]", []),
    ("Sphingosine 1-phosphate", "Sphingolipids [SP]", "Sphingoid base 1-phosphates [SP0105]", []),
    ("Sphinganine 1-phosphate", "Sphingolipids [SP]", "Sphingoid base 1-phosphates [SP0105]", []),
    ("Cer(d18:0/16:0)", "Sphingolipids [SP]", "Dihydroceramides [SP0203]", None),
    ("Cer(d18:1/16:0)", "Sphingolipids [SP]", "Ceramides [SP02]", None),
    ("Cer(d18:1/18:0)", "Sphingolipids [SP]", "Ceramides [SP02]", None),
    ("Cer(d18:1/20:0)", "Sphingolipids [SP]", "Ceramides [SP02]", None),
    ("Cer(d18:1/22:0)", "Sphingolipids [SP]", "Ceramides [SP02]", None),
    ("Cer(d18:1/24:0)", "Sphingolipids [SP]", "Ceramides [SP02]", None),
    ("Cer(d18:1/24:1)", "Sphingolipids [SP]", "Ceramides [SP02]", None),
    ("SM(d18:1/16:0)", "Sphingolipids [SP]", "Sphingomyelins [SP0301]", None),
    ("SM(d18:1/18:0)", "Sphingolipids [SP]", "Sphingomyelins [SP0301]", None),
    ("SM(d18:1/18:1)", "Sphingolipids [SP]", "Sphingomyelins [SP0301]", None),
    ("SM(d18:1/24:0)", "Sphingolipids [SP]", "Sphingomyelins [SP0301]", None),
    ("GlcCer(d18:1/16:0)", "Sphingolipids [SP]", "Hexosylceramides [SP0501]", None),
    ("GlcCer(d18:1/18:0)", "Sphingolipids [SP]", "Hexosylceramides [SP0501]", None),
    ("GlcCer(d18:1/22:0)", "Sphingolipids [SP]", "Hexosylceramides [SP0501]", None),
    ("GlcCer(d18:1/24:0)", "Sphingolipids [SP]", "Hexosylceramides [SP0501]", None),
    ("LacCer(d18:1/16:0)", "Sphingolipids [SP]", "Dihexosylceramides (lactosylceramides) [SP0501]", None),
    ("LacCer(d18:1/24:0)", "Sphingolipids [SP]", "Dihexosylceramides (lactosylceramides) [SP0501]", None),
    ("Ganglioside GM3 (d18:1/16:0)", "Sphingolipids [SP]", "Acidic glycosphingolipids (gangliosides) [SP06]", None),
    # Sterol lipids
    ("Cholesterol", "Sterol Lipids [ST]", "Sterols [ST01]", []),
    ("7alpha-Hydroxycholesterol", "Sterol Lipids [ST]", "Sterols [ST01]", []),
    ("27-Hydroxycholesterol", "Sterol Lipids [ST]", "Sterols [ST01]", []),
    ("7alpha-Hydroxy-4-cholesten-3-one", "Sterol Lipids [ST]", "Sterols [ST01]", []),
]
_LIPID_BILE_ACIDS_CLASS = "Bile acids and derivatives [ST04]"
_LIPID_STEROID_CLASS = "Steroids (C18/C19/C21) [ST02]"
_CHAIN_RE = re.compile(r"(\d+):(\d+)")


def _chains_from_shorthand(name: str):
    """Acyl chains of a shorthand lipid name, excluding the sphingoid base and 0:0 placeholders."""
    inner = re.search(r"\(([^)]*)\)", name)
    if not inner:
        return []
    parts = [p.strip() for p in inner.group(1).split("/")]
    chains = []
    for p in parts:
        if p.startswith(("d", "t", "m")):  # sphingoid base, e.g. d18:1
            continue
        m = _CHAIN_RE.search(p)
        if m and m.group(0) != "0:0":
            chains.append(f"{int(m.group(1))}:{int(m.group(2))}")
    return chains


def _build_lipid_library():
    sets = {}

    def add(set_name, hid):
        sets.setdefault(set_name, set()).add(hid)

    for name, cat, cls, chains in _LIPIDS:
        rec = ref.lookup_by_name(name)
        if rec is None:
            raise ValueError(f"Lipidomics library member not in reference: {name}")
        hid = rec["hmdb_id"]
        add(f"Category: {cat}", hid)
        add(f"Class: {cls}", hid)
        chains = _chains_from_shorthand(name) if chains is None else chains
        long_chains = [c for c in chains if int(c.split(":")[0]) >= 12]
        for c in long_chains:
            add(f"Acyl chain: contains {c}", hid)
        if long_chains:
            dbs = [int(c.split(":")[1]) for c in long_chains]
            if all(d == 0 for d in dbs):
                add("Acyl feature: saturated chains only (SFA)", hid)
            if any(d == 1 for d in dbs):
                add("Acyl feature: contains a monounsaturated chain (MUFA)", hid)
            if any(d >= 2 for d in dbs):
                add("Acyl feature: contains a polyunsaturated chain (PUFA, >=2 C=C)", hid)
            if any(int(c.split(":")[0]) >= 22 for c in long_chains):
                add("Acyl feature: very-long-chain (>=C22) species", hid)
        if "Acylcarnitines" in cls:
            n = int(chains[0].split(":")[0]) if chains else None
            # Chain-length classes for acylcarnitines (short C2-C5 incl. branched/dicarboxylic,
            # medium C6-C12, long C14-C18) as conventionally reported in newborn/clinical panels.
            if n is None or n <= 5:
                add("Acylcarnitines: short-chain (C2-C5)", hid)
            elif n <= 12:
                add("Acylcarnitines: medium-chain (C6-C12)", hid)
            else:
                add("Acylcarnitines: long-chain (C14-C18)", hid)
    for m in ref.METABOLITES:
        if m["primary_pathway"] == "Bile Acid Metabolism" and "cholesterol" not in m["name"].lower() \
                and "cholesten" not in m["name"].lower():
            add("Category: Sterol Lipids [ST]", m["hmdb_id"])
            add(f"Class: {_LIPID_BILE_ACIDS_CLASS}", m["hmdb_id"])
            nm = m["name"].lower()
            if nm.startswith("glyco") or nm.startswith("tauro"):
                add("Bile acids: glycine/taurine-conjugated", m["hmdb_id"])
            else:
                add("Bile acids: unconjugated", m["hmdb_id"])
        if m["primary_pathway"] == "Steroid Hormone Biosynthesis" and m["name"] != "Cholesterol":
            add("Category: Sterol Lipids [ST]", m["hmdb_id"])
            add(f"Class: {_LIPID_STEROID_CLASS}", m["hmdb_id"])
    return sets


# ===========================================================================
# 5. Chemical-class libraries (ClassyFire-style taxonomy)
# ===========================================================================
_OA, _LIP, _NUC, _OOX, _OHC, _ONC = ("Organic acids and derivatives", "Lipids and lipid-like molecules",
                                     "Nucleosides, nucleotides, and analogues", "Organic oxygen compounds",
                                     "Organoheterocyclic compounds", "Organic nitrogen compounds")
_T = {  # (superclass, class, subclass) per group of reference metabolites
    "tricarboxylic": (_OA, "Carboxylic acids and derivatives", "Tricarboxylic acids and derivatives"),
    "dicarboxylic": (_OA, "Carboxylic acids and derivatives", "Dicarboxylic acids and derivatives"),
    "amino": (_OA, "Carboxylic acids and derivatives", "Amino acids, peptides, and analogues"),
    "alpha_keto": (_OA, "Keto acids and derivatives", "Alpha-keto acids and derivatives"),
    "short_keto": (_OA, "Keto acids and derivatives", "Short-chain keto acids and derivatives"),
    "gamma_keto": (_OA, "Keto acids and derivatives", "Gamma-keto acids and derivatives"),
    "beta_keto": (_OA, "Keto acids and derivatives", "Beta-keto acids and derivatives"),
    "alpha_hydroxy": (_OA, "Hydroxy acids and derivatives", "Alpha hydroxy acids and derivatives"),
    "beta_hydroxy": (_OA, "Hydroxy acids and derivatives", "Beta hydroxy acids and derivatives"),
    "phosphate_ester": (_OA, "Organic phosphoric acids and derivatives", "Phosphate esters"),
    "sulfonic": (_OA, "Organic sulfonic acids and derivatives", "Organosulfonic acids and derivatives"),
    "carbohydrate": (_OOX, "Organooxygen compounds", "Carbohydrates and carbohydrate conjugates"),
    "polyol": (_OOX, "Organooxygen compounds", "Alcohols and polyols"),
    "carbonyl": (_OOX, "Organooxygen compounds", "Carbonyl compounds"),
    "quat": (_ONC, "Organonitrogen compounds", "Quaternary ammonium salts"),
    "amine": (_ONC, "Organonitrogen compounds", "Amines"),
    "purine_base": (_OHC, "Imidazopyrimidines", "Purines and purine derivatives"),
    "pyrimidine_base": (_OHC, "Diazines", "Pyrimidines and pyrimidine derivatives"),
    "imidazolidine": (_OHC, "Azolidines", "Imidazolidines"),
    "indole": (_OHC, "Indoles and derivatives", "Hydroxyindoles"),
    "flavin": (_OHC, "Pteridines and derivatives", "Flavin nucleotides"),
    "purine_ribonucleoside": (_NUC, "Purine nucleosides", "Purine ribonucleosides"),
    "purine_deoxynucleoside": (_NUC, "Purine nucleosides", "Purine 2'-deoxyribonucleosides"),
    "pyrimidine_ribonucleoside": (_NUC, "Pyrimidine nucleosides", "Pyrimidine ribonucleosides"),
    "pyrimidine_deoxynucleoside": (_NUC, "Pyrimidine nucleosides", "Pyrimidine 2'-deoxyribonucleosides"),
    "purine_ribonucleotide": (_NUC, "Purine nucleotides", "Purine ribonucleotides"),
    "purine_deoxynucleotide": (_NUC, "Purine nucleotides", "Purine deoxyribonucleotides"),
    "cyclic_nucleotide": (_NUC, "Cyclic purine nucleotides", "3',5'-cyclic purine nucleotides"),
    "pyrimidine_ribonucleotide": (_NUC, "Pyrimidine nucleotides", "Pyrimidine ribonucleotides"),
    "dinucleotide": (_NUC, "(5'->5')-dinucleotides", "(5'->5')-dinucleotides"),
    "coa_nucleotide": (_NUC, "Purine nucleotides", "Purine ribonucleoside 3',5'-bisphosphates"),
    "acyl_coa": (_LIP, "Fatty Acyls", "Fatty acyl thioesters"),
    "fatty_acid": (_LIP, "Fatty Acyls", "Fatty acids and conjugates"),
    "fatty_ester": (_LIP, "Fatty Acyls", "Fatty acid esters"),
    "glycerophosphocholine": (_LIP, "Glycerophospholipids", "Glycerophosphocholines"),
    "glycerophosphoethanolamine": (_LIP, "Glycerophospholipids", "Glycerophosphoethanolamines"),
    "glycerophosphate": (_LIP, "Glycerophospholipids", "Glycerophosphates"),
    "diradylglycerol": (_LIP, "Glycerolipids", "Diradylglycerols"),
    "triradylglycerol": (_LIP, "Glycerolipids", "Triradylglycerols"),
    "ceramide": (_LIP, "Sphingolipids", "Ceramides"),
    "phosphosphingolipid": (_LIP, "Sphingolipids", "Phosphosphingolipids"),
    "glycosphingolipid": (_LIP, "Sphingolipids", "Glycosphingolipids"),
    "cholestane": (_LIP, "Steroids and steroid derivatives", "Cholestane steroids"),
    "pregnane": (_LIP, "Steroids and steroid derivatives", "Pregnane steroids"),
    "androstane": (_LIP, "Steroids and steroid derivatives", "Androstane steroids"),
    "estrane": (_LIP, "Steroids and steroid derivatives", "Estrane steroids"),
    "sulfated_steroid": (_LIP, "Steroids and steroid derivatives", "Sulfated steroids"),
    "bile_acid": (_LIP, "Steroids and steroid derivatives", "Bile acids, alcohols and derivatives"),
}
_CHEM_GROUPS = {
    "tricarboxylic": ["Citric acid", "cis-Aconitic acid", "Isocitric acid", "Oxalosuccinic acid"],
    "dicarboxylic": ["Succinic acid", "Fumaric acid", "Itaconic acid", "Methylmalonic acid"],
    "gamma_keto": ["Oxoglutaric acid"],
    "short_keto": ["Oxaloacetic acid"],
    "alpha_keto": ["Pyruvic acid"],
    "beta_keto": ["Acetoacetic acid"],
    "alpha_hydroxy": ["L-Lactic acid", "D-2-Hydroxyglutaric acid"],
    "beta_hydroxy": ["L-Malic acid", "3-Hydroxybutyric acid"],
    "phosphate_ester": ["Phosphoenolpyruvic acid", "O-Phosphoethanolamine"],
    "sulfonic": ["Taurine"],
    "amino": ["L-Alanine", "L-Arginine", "L-Aspartic acid", "L-Glutamic acid", "L-Glutamine", "Glycine",
              "L-Histidine", "L-Isoleucine", "L-Leucine", "L-Lysine", "L-Methionine", "L-Phenylalanine",
              "L-Proline", "L-Serine", "L-Threonine", "L-Tryptophan", "L-Tyrosine", "L-Valine", "L-Ornithine",
              "L-Citrulline", "L-Asparagine", "L-Cysteine", "Creatine", "Betaine", "Kynurenine", "Hydroxyproline"],
    "carbohydrate": ["D-Glucose", "Glucose 6-phosphate", "Glucose 1-phosphate", "Fructose 6-phosphate",
                     "Fructose 1,6-bisphosphate", "Fructose 2,6-bisphosphate", "Dihydroxyacetone phosphate",
                     "Glyceraldehyde 3-phosphate", "1,3-Bisphosphoglyceric acid", "3-Phosphoglyceric acid",
                     "2-Phosphoglyceric acid", "2,3-Bisphosphoglyceric acid", "D-Fructose", "D-Mannose",
                     "Mannose 6-phosphate", "D-Galactose", "Galactose 1-phosphate", "Glyceric acid",
                     "Ribose 5-phosphate", "Sedoheptulose 7-phosphate", "Fructose 1-phosphate",
                     "Glucose 1,6-bisphosphate"],
    "polyol": ["Glycerol"],
    "carbonyl": ["Methylglyoxal"],
    "quat": ["Choline", "Phosphocholine", "L-Carnitine", "gamma-Butyrobetaine"],
    "amine": ["Ethanolamine", "Sphinganine", "Sphingosine", "3-Ketosphinganine"],
    "purine_base": ["Adenine", "Guanine", "Hypoxanthine", "Xanthine", "Uric acid"],
    "pyrimidine_base": ["Uracil", "Thymine", "Orotic acid"],
    "imidazolidine": ["Allantoin"],
    "indole": ["Serotonin"],
    "flavin": ["FAD", "FADH2"],
    "purine_ribonucleoside": ["Inosine", "Adenosine", "Guanosine", "Xanthosine"],
    "purine_deoxynucleoside": ["Deoxyadenosine", "Deoxyguanosine", "Deoxyinosine"],
    "pyrimidine_ribonucleoside": ["Uridine", "Cytidine"],
    "pyrimidine_deoxynucleoside": ["Thymidine"],
    "purine_ribonucleotide": ["Adenosine monophosphate", "ADP", "Adenosine triphosphate", "Inosinic acid",
                              "Guanosine monophosphate", "Guanosine diphosphate", "Guanosine triphosphate",
                              "Xanthylic acid"],
    "purine_deoxynucleotide": ["Deoxyadenosine monophosphate", "Deoxyadenosine triphosphate"],
    "cyclic_nucleotide": ["Cyclic AMP", "Cyclic GMP"],
    "pyrimidine_ribonucleotide": ["Uridine 5'-monophosphate", "Uridine 5'-diphosphate", "Uridine triphosphate",
                                  "Cytidine monophosphate", "Cytidine triphosphate", "Orotidylic acid",
                                  "CDP-choline", "CDP-ethanolamine"],
    "dinucleotide": ["NAD", "NADH"],
    "coa_nucleotide": ["Coenzyme A"],
    "acyl_coa": ["Succinyl-CoA", "Acetyl-CoA", "Palmitoyl-CoA", "Malonyl-CoA", "Acetoacetyl-CoA"],
    "fatty_acid": ["Arachidonic acid", "Linoleic acid", "Docosahexaenoic acid", "Palmitic acid", "Oleic acid",
                   "Myristic acid"],
    "fatty_ester": ["L-Acetylcarnitine", "Propionylcarnitine", "Malonylcarnitine", "Butyrylcarnitine",
                    "3-Hydroxybutyrylcarnitine", "Isovalerylcarnitine", "Glutarylcarnitine", "Hexanoylcarnitine",
                    "Octanoylcarnitine", "Decanoylcarnitine", "Dodecanoylcarnitine", "Tetradecanoylcarnitine",
                    "Palmitoylcarnitine", "Stearoylcarnitine", "Oleoylcarnitine", "Linoleyl carnitine",
                    "Valerylcarnitine", "Tiglylcarnitine", "Hexadecenoylcarnitine", "Tetradecenoylcarnitine"],
    "glycerophosphocholine": ["Glycerophosphocholine", "LysoPC(16:0)", "LysoPC(18:0)", "LysoPC(18:1)",
                              "LysoPC(18:2)", "LysoPC(20:4)", "PC(16:0/18:1)", "PC(16:0/18:2)", "PC(18:0/18:2)"],
    "glycerophosphoethanolamine": ["PE(18:0/20:4)"],
    "glycerophosphate": ["Glycerol 3-phosphate", "PA(16:0/18:1)"],
    "diradylglycerol": ["DG(16:0/18:1/0:0)"],
    "triradylglycerol": ["TG(16:0/18:1/18:1)"],
    "ceramide": ["Cer(d18:0/16:0)", "Cer(d18:1/16:0)", "Cer(d18:1/18:0)", "Cer(d18:1/20:0)", "Cer(d18:1/22:0)",
                 "Cer(d18:1/24:0)", "Cer(d18:1/24:1)"],
    "phosphosphingolipid": ["SM(d18:1/16:0)", "SM(d18:1/18:0)", "SM(d18:1/18:1)", "SM(d18:1/24:0)",
                            "Sphingosine 1-phosphate", "Sphinganine 1-phosphate"],
    "glycosphingolipid": ["GlcCer(d18:1/16:0)", "GlcCer(d18:1/18:0)", "GlcCer(d18:1/22:0)", "GlcCer(d18:1/24:0)",
                          "LacCer(d18:1/16:0)", "LacCer(d18:1/24:0)", "Ganglioside GM3 (d18:1/16:0)"],
    "cholestane": ["Cholesterol", "7alpha-Hydroxycholesterol", "27-Hydroxycholesterol",
                   "7alpha-Hydroxy-4-cholesten-3-one"],
    "pregnane": ["Pregnenolone", "17alpha-Hydroxypregnenolone", "Progesterone", "17alpha-Hydroxyprogesterone",
                 "11-Deoxycorticosterone", "Corticosterone", "18-Hydroxycorticosterone", "Aldosterone",
                 "11-Deoxycortisol", "Cortisol", "Cortisone", "Tetrahydrocortisone"],
    "androstane": ["Dehydroepiandrosterone", "Androstenedione", "Testosterone", "Dihydrotestosterone",
                   "Androsterone", "Epiandrosterone"],
    "estrane": ["Estrone", "Estradiol", "Estriol"],
    "sulfated_steroid": ["Pregnenolone sulfate", "Dehydroepiandrosterone sulfate", "Androsterone sulfate",
                         "Estrone sulfate"],
    "bile_acid": ["Cholic acid", "Chenodeoxycholic acid", "Glycocholic acid", "Taurocholic acid",
                  "Glycochenodeoxycholic acid", "Taurochenodeoxycholic acid", "Deoxycholic acid",
                  "Lithocholic acid", "Ursodeoxycholic acid", "Hyodeoxycholic acid", "Glycodeoxycholic acid",
                  "Taurodeoxycholic acid", "Glycolithocholic acid", "Taurolithocholic acid",
                  "Glycoursodeoxycholic acid", "Tauroursodeoxycholic acid", "Hyocholic acid",
                  "Isoursodeoxycholic acid"],
}


def chemical_taxonomy() -> dict:
    """{HMDB ID: (superclass, class, subclass)} for every reference metabolite."""
    out = {}
    for grp, names in _CHEM_GROUPS.items():
        for n in names:
            rec = ref.lookup_by_name(n)
            if rec is None:
                raise ValueError(f"Chemical-class member not in reference: {n}")
            if rec["hmdb_id"] in out:
                raise ValueError(f"Chemical-class member assigned twice: {n}")
            out[rec["hmdb_id"]] = _T[grp]
    return out


# ===========================================================================
# 6. Library registry
# ===========================================================================
LIB_KEGG = "KEGG"
LIB_KEGG_HMDB = "KEGG — HMDB map annotations (legacy)"
LIB_SMPDB = "SMPDB"
LIB_SMPDB_ALL = "SMPDB — All pathway types (incl. disease / drug / signaling)"
LIB_LIPID = "LIPID MAPS"
LIB_CHEM_SUPER = "Chemical Class — Superclass"
LIB_CHEM_CLASS = "Chemical Class — Main class"
LIB_CHEM_SUB = "Chemical Class — Subclass"
LIB_CUSTOM = "User-defined / custom metabolite sets"

BUILTIN_LIBRARIES = [LIB_KEGG, LIB_SMPDB, LIB_LIPID]
# Retained builders/constants (LIB_KEGG_HMDB, LIB_SMPDB_ALL, LIB_CHEM_SUPER/CLASS/SUB) are no longer
# offered in the UI dropdown -- see get_library() -- but are left in place (dead but harmless) in case
# they are wanted again later.

LIBRARY_DESCRIPTIONS = {
    LIB_KEGG: "KEGG human metabolic pathways (78 maps), keyed by KEGG compound ID (HMDB IDs via the RaMP-DB "
              "crosswalk).",
    LIB_KEGG_HMDB: "Legacy KEGG library: KEGG map IDs as annotated in HMDB's bulk pathway export, matched by "
                   "HMDB accession. 55 maps (3 retired in current KEGG); set sizes can differ from KEGG's own.",
    LIB_SMPDB: "SMPDB metabolic pathways from HMDB's bulk pathway export (68 pathways).",
    LIB_SMPDB_ALL: "Every SMPDB pathway type in HMDB's bulk export -- metabolic, disease, drug action, drug "
                   "metabolism and signaling/physiological.",
    LIB_LIPID: "LIPID MAPS categories and classes, plus acyl-chain features (SFA/MUFA/PUFA, very-long-chain, "
               "specific chains, acylcarnitine chain length).",
    LIB_CHEM_SUPER: "Chemical superclass (e.g. Lipids and lipid-like molecules).",
    LIB_CHEM_CLASS: "Chemical class (e.g. Steroids and steroid derivatives).",
    LIB_CHEM_SUB: "Chemical subclass (e.g. Bile acids, alcohols and derivatives).",
    LIB_CUSTOM: "Your own metabolite sets (CSV or GMT upload).",
}


class Library:
    """A metabolite-set collection: sets {name: frozenset(keys)}, set IDs, key type and provenance."""

    def __init__(self, name, sets, set_ids=None, key_type="canonical", description="", member_labels=None):
        self.name = name
        self.sets = {k: frozenset(v) for k, v in sets.items() if v}
        self.set_ids = set_ids or {k: k for k in self.sets}
        self.key_type = key_type
        self.description = description
        self.member_labels = member_labels or {}
        self.name_index = {}   # normalized name -> key, for metabolites that arrive with no ID
        self.secondary = {}    # secondary (merged) HMDB accession -> primary accession
        self.set_categories = {}
        self._universe = None
        # Reference-panel calibration (KEGG library only): true set sizes, library-wide universe size and
        # number of sets, used for the whole-library background ("All metabolites in selected library").
        self.size_override = {}
        self.universe_size_override = None
        self.n_sets_reference = None
        self.calibration = {}
        self.kegg = None  # KeggLibrary resolver for the KEGG-ID-keyed library

    @property
    def calibrated(self) -> bool:
        return self.universe_size_override is not None

    def _members(self) -> frozenset:
        if self._universe is None:
            self._universe = frozenset().union(*self.sets.values()) if self.sets else frozenset()
        return self._universe

    @property
    def universe(self) -> set:
        """Every member of every set (computed once; a fresh set is returned so callers may mutate it)."""
        return set(self._members())

    def resolve_key(self, key: str) -> str:
        """
        Library-specific key resolution on top of canonical_key: a secondary HMDB accession that is
        not itself a member becomes its primary accession, and a bare 'name:' key (no ID found) is
        looked up among the HMDB names of the library's members. Identity for other libraries.
        """
        if not key:
            return key
        if self.kegg is not None:
            return self.kegg.resolve_key(key)
        if key.startswith("name:"):
            return self.name_index.get(key[5:], key)
        prim = self.secondary.get(key)
        return prim if prim and key not in self._members() else key

    def __len__(self):
        return len(self.sets)


def _lipid_library():
    sets = _build_lipid_library()
    ids = {k: f"LIP-{i + 1:02d}" for i, k in enumerate(sorted(sets))}
    return Library(LIB_LIPID, sets, ids, "canonical", LIBRARY_DESCRIPTIONS[LIB_LIPID])


def _chem_library(level):
    tax = chemical_taxonomy()
    idx = {LIB_CHEM_SUPER: 0, LIB_CHEM_CLASS: 1, LIB_CHEM_SUB: 2}[level]
    prefix = {LIB_CHEM_SUPER: "CSP", LIB_CHEM_CLASS: "CCL", LIB_CHEM_SUB: "CSB"}[level]
    sets = {}
    for hid, t in tax.items():
        sets.setdefault(t[idx], set()).add(hid)
    ids = {k: f"{prefix}-{i + 1:02d}" for i, k in enumerate(sorted(sets))}
    return Library(level, sets, ids, "canonical", LIBRARY_DESCRIPTIONS[level])


_CACHE = {}


def get_library(name) -> Library:
    if name not in _CACHE:
        builders = {LIB_KEGG: _kegg_library, LIB_KEGG_HMDB: _kegg_hmdb_library,
                    LIB_SMPDB: _smpdb_library, LIB_SMPDB_ALL: _smpdb_all_library,
                    LIB_LIPID: _lipid_library,
                    LIB_CHEM_SUPER: lambda: _chem_library(LIB_CHEM_SUPER),
                    LIB_CHEM_CLASS: lambda: _chem_library(LIB_CHEM_CLASS),
                    LIB_CHEM_SUB: lambda: _chem_library(LIB_CHEM_SUB)}
        if name not in builders:
            raise KeyError(name)
        _CACHE[name] = builders[name]()
    return _CACHE[name]


# ===========================================================================
# 7. Custom uploads (library and reference metabolome)
# ===========================================================================
_SPLIT = re.compile(r"\s*[;|,\t]\s*")


def parse_custom_library(text: str, filename: str = "") -> Library:
    """
    Parse a user metabolite-set file into a Library. Accepted layouts:
      * GMT: Set_Name<TAB>description<TAB>member1<TAB>member2 ...
      * CSV long:  Set_Name,Metabolite    (one member per row)
      * CSV wide:  Set_Name,Metabolites   (members separated by ';' or '|')
    Members may be names, synonyms, HMDB, KEGG, ChEBI or PubChem IDs.
    """
    text = text.lstrip("﻿")
    raw_sets = {}
    lines = [l for l in text.splitlines() if l.strip()]
    if not lines:
        raise ValueError("The file is empty.")
    is_gmt = filename.lower().endswith(".gmt") or (lines[0].count("\t") >= 2 and "," not in lines[0])
    if is_gmt:
        for line in lines:
            parts = [p.strip() for p in line.split("\t")]
            if len(parts) >= 3 and parts[0]:
                raw_sets.setdefault(parts[0], []).extend(p for p in parts[2:] if p)
    else:
        rows = list(csv.reader(io.StringIO("\n".join(lines))))
        header = [h.strip().lower() for h in rows[0]]
        body = rows[1:] if any(h in ("set_name", "set", "pathway", "name", "set name") for h in header) else rows
        for r in body:
            if len(r) < 2 or not r[0].strip():
                continue
            members = [m for cell in r[1:] for m in _SPLIT.split(cell.strip()) if m.strip()]
            raw_sets.setdefault(r[0].strip(), []).extend(members)
    if not raw_sets:
        raise ValueError("No metabolite sets could be parsed (expected Set_Name + metabolite columns, or GMT).")
    sets, labels = {}, {}
    for s, members in raw_sets.items():
        keys = set()
        for mbr in members:
            k = resolve_member(mbr)
            if k:
                keys.add(k)
                labels.setdefault(k, mbr)
        if keys:
            sets[s] = keys
    ids = {k: f"USR-{i + 1:03d}" for i, k in enumerate(sets)}
    return Library(LIB_CUSTOM + (f": {filename}" if filename else ""), sets, ids, "canonical",
                   LIBRARY_DESCRIPTIONS[LIB_CUSTOM], labels)


def parse_reference_list(text: str) -> list:
    """Metabolite list (one per row; first column, header optional) for a custom reference metabolome."""
    text = text.lstrip("﻿")
    rows = list(csv.reader(io.StringIO(text)))
    out = []
    for i, r in enumerate(rows):
        if not r or not r[0].strip():
            continue
        v = r[0].strip()
        if i == 0 and v.lower() in ("metabolite", "metabolites", "name", "compound", "id", "hmdb_id",
                                    "kegg_id", "reference", "metabolite_name"):
            continue
        out.append(v)
    return list(dict.fromkeys(out))


def reference_keys(labels, key_type: str = "canonical", library=None) -> set:
    """
    Resolve a reference-metabolome list to library keys ('' / unresolvable dropped). With `library`,
    its resolve_key() is applied too (secondary HMDB accessions, HMDB-name fallback). Every built-in
    and custom library uses canonical keys (key_type is kept for backward compatibility).
    """
    keys = (resolve_member(l) for l in labels)
    if library is not None:
        keys = (library.resolve_key(k) for k in keys)
    return {k for k in keys if k}
