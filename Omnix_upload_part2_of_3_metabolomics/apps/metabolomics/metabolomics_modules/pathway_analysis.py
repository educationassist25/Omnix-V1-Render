"""
pathway_analysis.py - the two enrichment tabs of MetaboAI Pro.

Tab "Metabolite-Set Enrichment Analysis (MSEA)" (render_msea_tab) -- metabolite-set over-representation
analysis:
  1. Input Metabolite Selection -- from the existing Statistics-tab results only (significant by P-value /
     by FDR / all tested; optional Log2FC direction)
  2. ID Mapping -- HMDB_ID (annotation HMDB/KEGG/ChEBI/PubChem columns, IDs given as names, exact names,
     synonyms), plus SMPDB_ID and KEGG_Map_ID for that HMDB_ID from the complete HMDB compound-pathway
     export (pathway_libraries.add_hmdb_pathway_ids; metabolomics_modules/data/hmdb_pathway_id_crosswalk.parquet)
  3. Reference Metabolome -- the selected metabolite-set library (KEGG, SMPDB, LIPID MAPS or user-uploaded
     sets; see pathway_libraries.py) applied to ONE background by default: every metabolite in the selected
     library (REF_LIBRARY), the conventional default background for library-wide over-representation
     analysis when no custom reference metabolome is uploaded. A single off-by-default checkbox restricts
     this to only the metabolites detected/measured in the dataset (REF_DETECTED), for a more conservative,
     platform-specific analysis.
  4. Enrichment Analysis -- hypergeometric / one-sided Fisher; Holm and BH-FDR across every set of the
     library within the size limits (ADJUST_ALL: 0-hit sets count as P = 1 in the correction)
  5. Results -- "Pathway Results -- Library: <library> · Reference: <the reference actually used>",
     enrichment dot plot (x = -log10 P, colour = P-value, size = enrichment ratio; SVG download),
     CSV/Excel downloads.
  The engine still supports a custom reference (run_msea_ora, run_msea_ora_on_list), used by
  validate_pathway_parity.py to reproduce a curated reference-panel run.

Tab "Metabolite-Gene Pathway Analysis (MGPA)" (render_mgpa_tab):
  1. Input Metabolite Selection -> 2. HMDB Mapping -> 3. Metabolite-Gene Mapping & Unique Gene List (REAL HMDB
  protein associations, metabolomics_modules/hmdb_gene_mapping.py: 858,077 associations / 22,849 metabolites) ->
  4. Gene-Set ORA Settings, mirroring the Proteomics app's (ProteoAI Pro) "Over-Representation Analysis
  (ORA)": Organism, Gene-set database (GO BP/MF/CC, KEGG, Reactome, WikiPathways via Enrichr; MSigDB
  Hallmark bundled offline), Background / universe (all genes of the detected metabolites | entire
  gene-set library), Min / Max gene set size, Run button -> 5. Results (ProteoAI's columns, CSV, analysis
  metadata) -> 6. Pathway Dot Plot (ProteoAI's ora_dotplot) -> 7. Downloads. Engine:
  gene_set_enrichment.run_ora_hypergeometric (ProteoAI's, verbatim: hypergeometric, BH-FDR).

Statistics
----------
Over-representation is tested one-sided (enrichment only) with the
hypergeometric distribution, for each pathway:

    N = |background (universe)|           K = |pathway ∩ background|
    n = |input list ∩ background|         k = |input ∩ pathway ∩ background|
    P = P(X >= k) = hypergeom.sf(k - 1, N, K, n)

which is identical to a one-sided Fisher's exact test on the 2x2 table
[[k, n-k], [K-k, N-K-n+k]] (both are offered; they give the same p-value).
Enrichment ratio (fold enrichment) = (k / n) / (K / N) = k / Expected with Expected = n K / N; for
gene-based ORA the clusterProfiler-style Gene ratio = k / n is reported as well. MSEA P-values are
Holm- and Benjamini-Hochberg-adjusted (the app's stats_analysis._fdr), always corrected across every set
of the library within the size limits (0-hit sets = P 1) -- the engine also supports the tested-sets
family (ADJUST_TESTED). Equivalent to R's phyper(k - 1, K, N - K, n, lower.tail = FALSE), p.adjust
"holm" and "fdr"; with the whole-library background N = |library|, n = input compounds in the library.
MGPA uses ProteoAI's BH-FDR across the tested gene sets (gene_set_enrichment.py).

This module only READS the existing statistics results; it never recomputes
or modifies them, and it does not touch any other tab's session state.
"""

import math
import re
import textwrap

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.figure
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy.stats import hypergeom, fisher_exact

from metabolomics_modules import pathway_reference_db as ref
from metabolomics_modules import pathway_libraries as libs
from metabolomics_modules import hmdb_gene_mapping as hgm
from metabolomics_modules import stats_analysis, utils, dataset_manager
from metabolomics_modules import gene_set_enrichment as gse

matplotlib.use("Agg")

MODE_P = "Significant by P-value"
MODE_FDR = "Significant by FDR"
MODE_ALL = "All metabolites (from t-test/ANOVA)"

METHOD_A = "Metabolite-Set Enrichment Analysis (MSEA)"
METHOD_B = "Metabolite–Gene Pathway Analysis (MGPA)"

TEST_HYPERGEOM = "Hypergeometric test (one-sided)"
TEST_FISHER = "Fisher's exact test (one-sided)"

BG_DATASET = "dataset"
BG_REFERENCE = "reference"

LIST_SEP = "; "
STATUS_MAPPED = "Mapped"
STATUS_NO_GENES = "HMDB ID found — no protein/gene associations in HMDB"
STATUS_NO_PATHWAY = STATUS_NO_GENES  # backward-compatible alias (pre-HMDB-table name)
STATUS_UNMAPPED = "Unmapped (no HMDB ID)"
GENE_LABEL_SYMBOL = "Gene symbol"
GENE_LABEL_OTHER = "Non-symbol HMDB label (excluded from ORA)"

RESULT_COLUMNS_METABOLITE = ["Pathway_ID", "Pathway", "P-value", "Holm", "FDR", "Enrichment_Ratio", "Expected",
                             "Overlap_Count", "Pathway_Size", "Input_Size", "Background_Size",
                             "Overlapping_Metabolites"]
RESULT_COLUMNS_GENE = ["Pathway_ID", "Pathway", "P-value", "FDR", "Gene_Ratio", "Enrichment_Ratio",
                       "Overlap_Count", "Pathway_Size", "Input_Size", "Background_Size", "Overlapping_Genes"]

_HMDB_COL_CANDIDATES = ("hmdb_id", "hmdb id", "hmdb", "hmdbid", "hmdb_accession", "hmdb accession")
_DUP_SUFFIX = re.compile(r"^(.*) \((\d+)\)$")


class PathwayAnalysisError(Exception):
    pass


# ===========================================================================
# 1. Input metabolite selection (reads existing Statistics-tab results)
# ===========================================================================
def available_stats_sources(state) -> dict:
    """
    Every statistics result already computed in the Statistics tab, keyed by a
    display label. Nothing is recomputed. Each entry records which columns hold
    the raw and BH-adjusted p-values, so selection works identically for the
    two-group table ('p-value', FDR/BH label, Log2FC) and the ANOVA table
    ('ANOVA p-value', FDR/BH label).
    """
    sources = {}
    two = state.get("stats_result")
    groups = state.get("stats_result_groups")
    if two is not None and len(two) and groups:
        g_a, g_b = groups
        fdr_col = "BH P Value" if "BH P Value" in two.columns else "FDR"
        sources[f"Two-group test: {g_a} vs {g_b}"] = {
            "df": two, "p_col": "p-value", "fdr_col": fdr_col, "kind": "two_group",
            "comparison": f"{g_a}_vs_{g_b}", "fc_col": "Log2FC" if "Log2FC" in two.columns else None,
            "groups": (g_a, g_b),
        }
    anova = state.get("anova_result")
    if anova is not None and len(anova):
        used = state.get("anova_groups_used") or []
        fdr_col = "BH P Value" if "BH P Value" in anova.columns else "FDR"
        label = "ANOVA: " + (", ".join(map(str, used)) if used else "multi-group")
        sources[label] = {
            "df": anova, "p_col": "ANOVA p-value", "fdr_col": fdr_col, "kind": "anova",
            "comparison": "ANOVA_" + "_vs_".join(map(str, used)) if used else "ANOVA",
            "fc_col": None, "groups": tuple(used),
        }
    return sources


def select_metabolites(stats_df: pd.DataFrame, p_col: str, fdr_col: str, mode: str,
                       threshold: float = 0.05, direction: str = "both", fc_col: str = None) -> pd.DataFrame:
    """
    Return the subset of the statistics table used as ORA input.
      mode MODE_P   : p-value <= threshold
      mode MODE_FDR : BH-adjusted p-value <= threshold
      mode MODE_ALL : every tested metabolite (no filter)
    Thresholds are inclusive (<=), matching the Heatmap tab's cutoff convention.
    direction ('both' | 'up' | 'down') optionally restricts a two-group result
    by the sign of Log2FC. NaN p-values never pass a significance filter.
    """
    df = stats_df.copy()
    if mode == MODE_P:
        mask = df[p_col] <= threshold
    elif mode == MODE_FDR:
        mask = df[fdr_col] <= threshold
    else:
        mask = pd.Series(True, index=df.index)
    if fc_col and fc_col in df.columns and direction in ("up", "down"):
        mask &= (df[fc_col] > 0) if direction == "up" else (df[fc_col] < 0)
    return df[mask.fillna(False)]


# ===========================================================================
# 2. HMDB mapping
# ===========================================================================
def normalize_hmdb_id(value):
    """'HMDB00094' / 'hmdb0000094' / ' HMDB0000094 ' -> 'HMDB0000094'; anything else -> None."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return None
    m = re.fullmatch(r"\s*HMDB(\d{1,7})\s*", str(value), flags=re.IGNORECASE)
    return f"HMDB{m.group(1).zfill(7)}" if m else None


def find_hmdb_column(row_annotations):
    """Name of the HMDB accession column in a row-annotation table (case/space-insensitive), or None."""
    if row_annotations is None:
        return None
    for col in row_annotations.columns:
        if str(col).strip().lower() in _HMDB_COL_CANDIDATES:
            return col
    return None


def _annotation_hmdb_lookup(row_annotations) -> dict:
    col = find_hmdb_column(row_annotations)
    if col is None:
        return {}
    lookup = {}
    for idx, val in row_annotations[col].items():
        hid = normalize_hmdb_id(val)
        if hid:
            lookup.setdefault(str(idx), hid)
            lookup.setdefault(utils.display_feature_name(idx), hid)
    return lookup


def _base_name(name: str) -> str:
    """Undo utils.with_display_feature_names' ' (2)' de-duplication suffix."""
    m = _DUP_SUFFIX.match(str(name))
    return m.group(1) if m else str(name)


HMDB_MAPPING_COLUMNS = ["Metabolite", "HMDB_ID", "HMDB_Name", "N_HMDB_Proteins", "N_HMDB_Genes",
                        "Reference_Name", "Mapped_Pathways", "N_Pathways", "Mapping_Source", "Mapping_Status",
                        "Note"]


def map_metabolites_to_hmdb(names, row_annotations=None) -> pd.DataFrame:
    """
    Option B step 2: map metabolite names to HMDB accessions and to their REAL HMDB
    protein/gene associations (metabolomics_modules/hmdb_gene_mapping.py -- the bundled HMDB export,
    858,077 associations for 22,849 metabolites).

    Accession precedence (reported per row in Mapping_Source):
      1. the loaded row-annotation table's HMDB ID column (authoritative when present);
      2. otherwise the curated reference by name/synonym (case/punctuation-insensitive,
         L-/D- prefix and '-ic acid'/'-ate' tolerant);
      3. otherwise an exact (normalized) match to an HMDB common name in the HMDB export.
    A secondary (merged) HMDB accession is resolved to its primary accession (Note column).
    Status: Mapped = the accession has >= 1 protein association in HMDB; an accession with none
    is kept but flagged (STATUS_NO_GENES); a name with no accession at all is Unmapped.
    Mapped_Pathways are the curated reference's (Option A-style) pathway categories, shown for
    information only -- Option B's genes come exclusively from the HMDB table.
    """
    ann = _annotation_hmdb_lookup(row_annotations)
    rows = []
    for name in names:
        name = str(name)
        base = _base_name(name)
        ann_id = ann.get(name) or ann.get(base) or ann.get(utils.display_feature_name(base))
        ref_rec = ref.lookup_by_name(base)
        notes = []
        if ann_id:
            given, source = ann_id, "Row annotation (HMDB_ID column)"
            if ref_rec is not None and ref_rec["hmdb_id"] != ann_id:
                notes.append(f"Annotation HMDB ID used; curated name match would give "
                             f"{ref_rec['hmdb_id']} ({ref_rec['name']})")
        elif ref_rec is not None:
            given, source = ref_rec["hmdb_id"], "Curated reference (name/synonym match)"
        else:
            given = hgm.lookup_by_name(base)
            source = "HMDB common-name match" if given else "—"
        hid, sec_note = hgm.resolve_accession(given) if given else (None, "")
        if sec_note:
            notes.append(sec_note)
        rec = ref.BY_HMDB.get(hid) or (ref.BY_HMDB.get(given) if given else None)
        pathways = rec["pathways"] if rec is not None else []
        if hid and hgm.has_associations(hid):
            status = STATUS_MAPPED
            assoc = hgm.associations(hid)
            n_prot, n_gene = len(assoc), assoc["Gene"].nunique()
        else:
            status = STATUS_NO_GENES if hid else STATUS_UNMAPPED
            n_prot = n_gene = 0
        rows.append({
            "Metabolite": name, "HMDB_ID": hid or "", "HMDB_Name": hgm.hmdb_name(hid) if hid else "",
            "N_HMDB_Proteins": n_prot, "N_HMDB_Genes": n_gene,
            "Reference_Name": rec["name"] if rec is not None else "",
            "Mapped_Pathways": LIST_SEP.join(pathways), "N_Pathways": len(pathways),
            "Mapping_Source": source, "Mapping_Status": status, "Note": "; ".join(notes),
        })
    return pd.DataFrame(rows, columns=HMDB_MAPPING_COLUMNS)


def _hmdb_to_input_names(mapping_df: pd.DataFrame) -> dict:
    out = {}
    for _, r in mapping_df[mapping_df["HMDB_ID"].isin(ref.BY_HMDB)].iterrows():
        out.setdefault(r["HMDB_ID"], []).append(r["Metabolite"])
    return out


# ===========================================================================
# 3. ORA engine (shared by Option A and Option B)
# ===========================================================================
def hypergeom_pvalue(k: int, N: int, K: int, n: int) -> float:
    """One-sided P(X >= k), X ~ Hypergeometric(N, K, n)."""
    if k <= 0:
        return 1.0
    return float(min(1.0, max(0.0, hypergeom.sf(k - 1, N, K, n))))


def fisher_pvalue(k: int, N: int, K: int, n: int) -> float:
    """One-sided (greater) Fisher's exact test on [[k, n-k], [K-k, N-K-n+k]] -- equals hypergeom_pvalue."""
    if k <= 0:
        return 1.0
    table = [[k, n - k], [K - k, N - K - n + k]]
    return float(fisher_exact(table, alternative="greater")[1])


def hypergeom_pvalue_exact(k: int, N: int, K: int, n: int) -> float:
    """Exact combinatorial form, sum_{i>=k} C(K,i) C(N-K,n-i) / C(N,n) -- used for validation."""
    total = math.comb(N, n)
    return sum(math.comb(K, i) * math.comb(N - K, n - i) for i in range(k, min(n, K) + 1)) / total


def holm_adjust(pvals, m: int = None) -> np.ndarray:
    """
    Holm step-down adjusted p-values (R p.adjust(method="holm")): sorted ascending,
    p_(i) * (m - i + 1), running maximum, capped at 1. m defaults to len(pvals); a larger m
    treats the extra (untested) hypotheses as p = 1, exactly as R does when they are included.
    """
    p = np.asarray(pvals, dtype=float)
    k = len(p)
    if k == 0:
        return p
    m = max(int(m or k), k)
    order = np.argsort(p, kind="mergesort")
    adj = np.minimum(1.0, np.maximum.accumulate((m - np.arange(k)) * p[order]))
    out = np.empty(k)
    out[order] = adj
    return out


def bh_adjust(pvals, m: int = None) -> np.ndarray:
    """Benjamini-Hochberg via stats_analysis._fdr, over m hypotheses (extra ones = p 1, as in R)."""
    p = np.asarray(pvals, dtype=float)
    k = len(p)
    m = max(int(m or k), k)
    if k == 0:
        return p
    return np.asarray(stats_analysis._fdr(np.concatenate([p, np.ones(m - k)])))[:k]


ADJUST_TESTED = "tested"   # correct across sets with >= 1 overlapping input
ADJUST_ALL = "all"         # correct across every set within the size limits, incl. 0-hit sets


def run_ora(input_items, universe, sets: dict, set_ids: dict = None, min_size: int = 2,
            max_size: int = None, test: str = TEST_HYPERGEOM, item_labels: dict = None,
            gene_mode: bool = False, set_size_override: dict = None, universe_size: int = None,
            adjust: str = ADJUST_TESTED, n_sets_family: int = None):
    """
    Generic over-representation analysis.

    input_items : iterable of identifiers (HMDB IDs or gene symbols)
    universe    : iterable -- the background; inputs and sets are intersected with it
    sets        : {pathway name: iterable of identifiers}
    set_ids     : {pathway name: pathway ID} (defaults to the name itself)
    item_labels : {identifier: display label} for the overlap column
    set_size_override : {pathway name: K} -- known true set sizes used instead of |set ∩ universe|
                  (the KEGG library's reference-panel calibration; whole-library background only)
    universe_size : N to use instead of |universe| (same calibration)
    adjust      : ADJUST_TESTED -- FDR/Holm across the sets tested (>= 1 overlap);
                  ADJUST_ALL -- across every set within the size limits (0-overlap sets count as
                  P = 1); n_sets_family, if larger, sets that count (m)
    Returns (results DataFrame sorted by P-value, stats dict).
    """
    U = set(universe)
    S = set(input_items) & U
    N, n = (int(universe_size) if universe_size else len(U)), len(S)
    set_size_override = set_size_override or {}
    set_ids = set_ids or {}
    item_labels = item_labels or {}
    pfun = fisher_pvalue if test == TEST_FISHER else hypergeom_pvalue
    rows = []
    n_size_filtered = n_no_overlap = 0
    for name, members in sets.items():
        G = set(members) & U
        K = int(set_size_override.get(name, len(G)))
        if K == 0:
            continue
        if K < min_size or (max_size and K > max_size):
            n_size_filtered += 1
            continue
        hits = S & G
        k = len(hits)
        if k == 0:
            n_no_overlap += 1
            continue
        fold = (k / n) / (K / N) if n and N else np.nan
        row = {"Pathway_ID": set_ids.get(name, name), "Pathway": name, "P-value": pfun(k, N, K, n)}
        if gene_mode:
            row["Gene_Ratio"] = k / n
        else:
            row["Expected"] = n * K / N if N else np.nan
        row.update({
            "Enrichment_Ratio": fold, "Overlap_Count": k, "Pathway_Size": K, "Input_Size": n,
            "Background_Size": N,
            ("Overlapping_Genes" if gene_mode else "Overlapping_Metabolites"):
                LIST_SEP.join(sorted(item_labels.get(h, h) for h in hits)),
        })
        rows.append(row)
    cols = RESULT_COLUMNS_GENE if gene_mode else RESULT_COLUMNS_METABOLITE
    res = pd.DataFrame(rows)
    m = len(rows) if adjust != ADJUST_ALL else max(len(rows) + n_no_overlap, int(n_sets_family or 0))
    if len(res):
        if gene_mode:
            res["FDR"] = stats_analysis._fdr(res["P-value"].values)
        else:
            res["FDR"] = bh_adjust(res["P-value"].values, m)
            res["Holm"] = holm_adjust(res["P-value"].values, m)
        res = res[cols].sort_values(["P-value", "Enrichment_Ratio"], ascending=[True, False],
                                    kind="mergesort").reset_index(drop=True)
    else:
        res = pd.DataFrame(columns=cols)
    info = {"N": N, "n": n, "n_input_total": len(set(input_items)), "n_sets_total": len(sets),
            "n_sets_tested": len(res), "n_sets_size_filtered": n_size_filtered,
            "n_sets_no_overlap": n_no_overlap, "n_adjusted": m if len(res) else 0, "adjust": adjust,
            "universe_members": len(U)}
    return res, info


def metabolite_pathway_sets(background: str, background_hmdb=None) -> dict:
    """{pathway: set(HMDB IDs)} from the embedded reference (optionally restricted to a background)."""
    sets = {}
    for m in ref.METABOLITES:
        if background == BG_DATASET and background_hmdb is not None and m["hmdb_id"] not in background_hmdb:
            continue
        for pw in m["pathways"]:
            sets.setdefault(pw, set()).add(m["hmdb_id"])
    return sets


def run_metabolite_ora(selected_map: pd.DataFrame, background_map: pd.DataFrame, background: str = BG_DATASET,
                       min_size: int = 2, max_size: int = None, test: str = TEST_HYPERGEOM):
    """Legacy metabolite-level ORA on the curated reference's pathway categories (not used by the UI)."""
    sel_ids = set(selected_map.loc[selected_map["HMDB_ID"].isin(ref.BY_HMDB), "HMDB_ID"])
    if background == BG_DATASET:
        universe = set(background_map.loc[background_map["HMDB_ID"].isin(ref.BY_HMDB), "HMDB_ID"])
    else:
        universe = {m["hmdb_id"] for m in ref.METABOLITES if m["pathways"]}
    sets = metabolite_pathway_sets(BG_REFERENCE)
    set_ids = {pw: ref.PATHWAYS.get(pw, {}).get("id", pw) for pw in sets}
    labels = {h: " / ".join(v) for h, v in _hmdb_to_input_names(selected_map).items()}
    return run_ora(sel_ids, universe, sets, set_ids=set_ids, min_size=min_size, max_size=max_size,
                   test=test, item_labels=labels, gene_mode=False)


# ===========================================================================
# 3a'. MSEA -- metabolite-set over-representation analysis (library + reference metabolome)
# ===========================================================================
REF_DETECTED = "All detected/measured metabolites"
REF_LIBRARY = "All metabolites in selected library (default)"
REF_CUSTOM = "Upload custom reference metabolome"


def standardize_for_option_a(names, row_annotations=None) -> pd.DataFrame:
    """Option A step 2: multi-ID standardization (names/synonyms/HMDB/KEGG/ChEBI/PubChem)."""
    return libs.standardize_ids(names, row_annotations, base_name_fn=_base_name,
                                display_name_fn=utils.display_feature_name)


def library_membership(id_map: pd.DataFrame, library) -> pd.Series:
    """Per-row library key ('' if the metabolite has no key or is not in any set of the library)."""
    keys = libs.mapping_keys(id_map, library.key_type, library)
    univ = library.universe
    return keys.where(keys.isin(univ), "")


def run_msea_ora(sel_idmap: pd.DataFrame, bg_idmap: pd.DataFrame, library, reference: str = REF_DETECTED,
                 custom_reference_keys=None, min_size: int = 2, max_size: int = None,
                 test: str = TEST_HYPERGEOM, adjust: str = ADJUST_TESTED):
    """
    MSEA: over-representation of the selected metabolites in each set of `library`.

    Universe (reference metabolome), per `reference`:
      REF_DETECTED : library keys of every tested metabolite in the dataset (bg_idmap)
      REF_LIBRARY  : every member of every set of the library (default)
      REF_CUSTOM   : custom_reference_keys  ∩  library members
    Input = selected metabolites' library keys ∩ universe. Sets are intersected with the universe
    (Pathway_Size = K within the background), then run_ora applies the hypergeometric / one-sided
    Fisher test, Holm and BH-FDR (across the sets tested, or with adjust=ADJUST_ALL across every
    set of the library within the size limits).

    Reference-panel calibration (library.calibrated, i.e. the KEGG library) applies to REF_LIBRARY
    only: N = the calibrated library universe (1519 compounds), K = the true size of every calibrated
    pathway (other sets: their KEGG size), and with ADJUST_ALL the correction runs over the calibrated
    81 sets: phyper(k - 1, K, N - K, n, lower.tail = FALSE) with p.adjust(..., "holm") / p.adjust(...,
    "fdr").
    """
    lib_univ = library.universe
    sel_keys = library_membership(sel_idmap, library)
    bg_keys = library_membership(bg_idmap, library)
    detected = {k for k in bg_keys if k}
    if reference == REF_LIBRARY:
        universe = set(lib_univ)
    elif reference == REF_CUSTOM:
        universe = set(custom_reference_keys or ()) & lib_univ
    else:
        universe = detected
    input_keys = {k for k in sel_keys if k}
    labels = {}
    member_names = getattr(library, "member_labels", {}) or {}
    for k, nm, matched in zip(sel_keys, sel_idmap["Metabolite"], sel_idmap["Matched_Name"]):
        if k:
            nm = str(nm)
            if libs._identifier_in_name(nm)[1]:  # an ID given as the label -> show the compound name too
                cname = matched or member_names.get(k, "")
                nm = f"{nm} ({cname})" if cname and cname != nm else nm
            labels.setdefault(k, []).append(nm)
    labels = {k: " / ".join(v) for k, v in labels.items()}
    calibrated = reference == REF_LIBRARY and getattr(library, "calibrated", False)
    res, info = run_ora(input_keys, universe, library.sets, set_ids=library.set_ids, min_size=min_size,
                        max_size=max_size, test=test, item_labels=labels, gene_mode=False,
                        set_size_override=library.size_override if calibrated else None,
                        universe_size=library.universe_size_override if calibrated else None,
                        adjust=adjust, n_sets_family=library.n_sets_reference if calibrated else None)
    info.update({
        "calibrated": calibrated,
        "library": library.name, "reference": reference, "n_library_members": len(lib_univ),
        "n_detected_in_library": len(detected), "n_selected_in_library": len(input_keys),
        "n_selected_not_in_universe": len(input_keys - universe),
        "selected_not_in_universe": sorted(labels.get(k, k) for k in input_keys - universe),
        "n_custom_reference_in_library": len(universe) if reference == REF_CUSTOM else None,
        "universe_keys": universe,
    })
    return res, info


def parse_metabolite_list(text: str) -> list:
    """
    Pasted metabolite list -> unique labels. One per line (or tab / semicolon separated); commas separate
    only lines made purely of IDs (e.g. "HMDB0000094, C00158"), since names such as
    "Fructose 1,6-bisphosphate" contain commas.
    """
    items = []
    for line in re.split(r"[\r\n\t;]+", str(text or "")):
        parts = [p.strip() for p in line.split(",")]
        ids_only = len(parts) > 1 and all(libs._identifier_in_name(p)[1] for p in parts if p)
        items.extend(parts if ids_only else [line])
    items = [t.strip().strip('"').strip("'").strip() for t in items]
    items = [t for t in items if t and t.lower() not in ("x", "metabolite", "compound", "hmdb", "hmdb_id", "name")]
    return list(dict.fromkeys(items))


def run_msea_ora_on_list(metabolites, library_name: str = None, library=None, reference: str = REF_LIBRARY,
                         custom_reference_keys=None, min_size: int = 2, max_size: int = None,
                         test: str = TEST_HYPERGEOM, adjust: str = ADJUST_ALL):
    """
    MSEA on an arbitrary metabolite list (names, synonyms, HMDB/KEGG/ChEBI/PubChem IDs), bypassing the
    Statistics-tab selection -- the same pipeline from ID standardization onward (validation runs). The
    list itself is the 'detected' background. Defaults = plain compound-list ORA (whole-library
    background, correction over all library sets).
    Returns (results, info, id_map).
    """
    idm = standardize_for_option_a(list(metabolites))
    lib = library if library is not None else libs.get_library(library_name or libs.LIB_KEGG)
    res, info = run_msea_ora(idm, idm, lib, reference, custom_reference_keys, min_size, max_size, test, adjust)
    return res, info, idm


MG_COLUMNS = ["Metabolite", "HMDB_ID", "HMDB_Name", "Gene", "Protein", "UniProt", "Gene_Label"]
UG_COLUMNS = ["Gene", "Symbol_Used_For_ORA", "Protein", "UniProt", "Associated_Metabolites", "N_Metabolites",
              "In_Gene_Set_Collection"]


def build_metabolite_gene_tables(mapping_df: pd.DataFrame, collection_genes=None):
    """
    Option B, step 3: Metabolite -> HMDB -> Gene, from the REAL HMDB association table.

    metabolite_gene_df: one row per HMDB metabolite-protein association of every Mapped
        metabolite -- exactly the rows HMDB lists (e.g. 70 for HMDB0004952) -- with columns
        Metabolite, HMDB_ID, HMDB_Name, Gene, Protein, UniProt, Gene_Label. Gene_Label flags the
        few HMDB entries whose "gene" is not a gene symbol (GENE_LABEL_OTHER); they stay in this
        table but are excluded from the Unique Gene List / ORA input.
    unique_gene_df: one row per distinct HMDB gene symbol -- Gene (as in HMDB), Symbol_Used_For_ORA
        (the form present in the selected collection: identical unless HMDB uses a retired HGNC
        symbol and the collection the current one, or vice versa; see
        hmdb_gene_mapping.symbol_for_collection), Protein, UniProt, Associated_Metabolites,
        N_Metabolites, In_Gene_Set_Collection.
    no_gene: mapped metabolites whose HMDB associations contain no gene symbol at all.
    """
    coll = {str(g).upper() for g in collection_genes} if collection_genes is not None else None
    frames, no_gene = [], []
    for _, r in mapping_df[mapping_df["Mapping_Status"] == STATUS_MAPPED].iterrows():
        a = hgm.associations(r["HMDB_ID"])
        a.insert(0, "HMDB_Name", hgm.hmdb_name(r["HMDB_ID"]))
        a.insert(0, "HMDB_ID", r["HMDB_ID"])
        a.insert(0, "Metabolite", r["Metabolite"])
        a["Gene_Label"] = np.where(a["Gene"].map(hgm.is_gene_symbol), GENE_LABEL_SYMBOL, GENE_LABEL_OTHER)
        if not (a["Gene_Label"] == GENE_LABEL_SYMBOL).any():
            no_gene.append(r["Metabolite"])
        frames.append(a)
    mg = pd.concat(frames, ignore_index=True)[MG_COLUMNS] if frames else pd.DataFrame(columns=MG_COLUMNS)
    sym = mg[mg["Gene_Label"] == GENE_LABEL_SYMBOL]
    if len(sym):
        grp = sym.groupby("Gene", sort=True)

        def _join(s):
            return LIST_SEP.join(dict.fromkeys(s))
        ug = pd.DataFrame({
            "Protein": grp["Protein"].agg(_join), "UniProt": grp["UniProt"].agg(_join),
            "Associated_Metabolites": grp["Metabolite"].agg(lambda s: LIST_SEP.join(sorted(dict.fromkeys(s)))),
            "N_Metabolites": grp["Metabolite"].nunique(),
        }).reset_index()
        ug.insert(1, "Symbol_Used_For_ORA", ug["Gene"].map(lambda g: hgm.symbol_for_collection(g, coll)))
        ug["In_Gene_Set_Collection"] = ug["Symbol_Used_For_ORA"].map(
            lambda g: "Yes" if coll is None or g in coll else "No")
        ug = ug.sort_values(["N_Metabolites", "Gene"], ascending=[False, True]).reset_index(drop=True)[UG_COLUMNS]
    else:
        ug = pd.DataFrame(columns=UG_COLUMNS)
    return mg, ug, no_gene


def genes_for_mapping(mapping_df: pd.DataFrame, collection_genes=None) -> set:
    """Distinct HMDB gene symbols of every Mapped metabolite (in collection form if a collection is given)."""
    coll = {str(g).upper() for g in collection_genes} if collection_genes is not None else None
    ids = mapping_df.loc[mapping_df["Mapping_Status"] == STATUS_MAPPED, "HMDB_ID"]
    out = set()
    for h in dict.fromkeys(ids):
        for g in hgm.genes(h, symbols_only=True):
            out.add(hgm.symbol_for_collection(g, coll) if coll is not None else str(g).upper())
    return out


def run_gene_ora(selected_map: pd.DataFrame, background_map: pd.DataFrame, gene_sets: dict,
                 set_ids: dict = None, background: str = BG_DATASET, min_size: int = 3,
                 max_size: int = None, test: str = TEST_HYPERGEOM):
    """Legacy (not used by the UI): unique gene list -> ORA against a gene-set collection."""
    gene_sets = {k: {str(g).upper() for g in v} for k, v in gene_sets.items()}
    collection = set().union(*gene_sets.values()) if gene_sets else set()
    selected_genes = genes_for_mapping(selected_map, collection)
    if background == BG_DATASET:
        universe = genes_for_mapping(background_map, collection) & collection
    else:
        universe = collection
    res, info = run_ora(selected_genes, universe, gene_sets, set_ids=set_ids, min_size=min_size,
                        max_size=max_size, test=test, gene_mode=True)
    info["unmapped_genes"] = sorted(selected_genes - collection)
    info["n_selected_genes_total"] = len(selected_genes)
    return res, info


def parse_gmt_text(text: str) -> dict:
    """Standard GMT (term<TAB>description<TAB>gene1<TAB>...) -> {term: [GENES]} (symbols upper-cased)."""
    sets = {}
    for line in text.strip().splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue
        genes = [g.strip().upper() for g in parts[2:] if g.strip()]
        if genes:
            sets[parts[0].strip()] = genes
    return sets


# ===========================================================================
# 4. Dot plot
# ===========================================================================
SORT_OPTIONS = {"P-value": ("P-value", True), "FDR": ("FDR", True), "Enrichment ratio": ("Enrichment_Ratio", False)}


def order_results(results: pd.DataFrame, sort_by: str = "P-value", top_n=None) -> pd.DataFrame:
    """Sort (ties broken by P-value) and keep the top_n rows (None = all)."""
    col, asc = SORT_OPTIONS[sort_by]
    keys, ascs = [col], [asc]
    if col != "P-value":
        keys.append("P-value")
        ascs.append(True)
    out = results.sort_values(keys, ascending=ascs, kind="mergesort")
    return out.head(top_n) if top_n else out


def neglog10(series: pd.Series) -> pd.Series:
    """-log10 with exact zeros floored to the smallest positive float (never +inf)."""
    return -np.log10(series.astype(float).clip(lower=np.finfo(float).tiny))


def r_heat_colors(n: int) -> list:
    """R grDevices::heat.colors(n)-equivalent gradient: red -> yellow -> pale yellow."""
    import colorsys
    n = max(int(n), 1)
    j = n // 4
    i = n - j
    cols = [colorsys.hsv_to_rgb((k / (i - 1) if i > 1 else 0.0) / 6.0, 1.0, 1.0) for k in range(i)]
    if j:
        sat = np.linspace(1 - 1 / (2 * j), 1 / (2 * j), j)
        cols += [colorsys.hsv_to_rgb(1 / 6.0, float(v), 1.0) for v in sat]
    return [matplotlib.colors.to_hex(c) for c in cols]


def enrichment_heat_cmap(n_results: int, top_n: int = 25):
    """Dot-plot colour scale: heat.colors(n_results)[1 : top_n] -- red at P = 0 to the colour of the
    last plotted rank at the largest plotted P."""
    cols = r_heat_colors(n_results)[:max(1, min(top_n, n_results))]
    lo, hi = cols[0], cols[-1] if len(cols) > 1 else "#FFA500"
    return matplotlib.colors.LinearSegmentedColormap.from_list("enrichment_heat", [lo, hi])


def pathway_dotplot(results: pd.DataFrame, x_col: str = "neglog10_p", x_label: str = "-log10 (p-value)",
                    top_n=25, sort_by: str = "P-value", title: str = None, cmap="viridis",
                    size_label: str = "Enrichment Ratio", reference_line: float = None,
                    size_col: str = "Enrichment_Ratio", color_by: str = "pvalue",
                    label_wrap: int = 40, max_label_lines: int = 3, x_from_zero: bool = False):
    """
    Metabolite-set enrichment dot plot (MSEA tab), sized as a compact figure (~half the footprint of a
    full-page plot) so it sits neatly in the results panel.
      Y-axis : metabolite set (ordered by `sort_by`; best at the top). Long names wrap onto up to
               `max_label_lines` lines (never truncated with "..." unless a single name exceeds that many
               lines) and each row is given exactly the vertical space its own line count needs, so a
               multi-line label never overlaps or crowds out a neighbouring row's label.
      X-axis : x_col -- default "neglog10_p" (-log10 P-value)
      Colour : color_by -- default "pvalue" (raw P, 0 .. largest plotted P); "neglog10_fdr" is kept for callers
      Size   : size_col -- default Enrichment_Ratio (hits / expected), with a compact size legend whose dots
               are spaced by their own diameter plus a small fixed gap, so they never overlap.
    The legend column (size key + colour scale) is a fixed height in inches, vertically centered, so it
    stays compact and legible regardless of how many pathways are plotted. Visual conventions follow the
    rest of MetaboAI Pro (sans-serif, bold titles, top/right spines removed).
    """
    if results is None or len(results) == 0:
        raise PathwayAnalysisError("No pathway results to plot.")
    df = order_results(results, sort_by, top_n).copy()
    df["neglog10_fdr"] = neglog10(df["FDR"])
    df["neglog10_p"] = neglog10(df["P-value"])
    n = len(df)

    def _wrap(name: str) -> list:
        lines = textwrap.wrap(str(name), label_wrap) or [str(name)]
        if len(lines) > max_label_lines:
            lines = lines[:max_label_lines]
            lines[-1] = lines[-1].rstrip() + " …"
        return lines

    line_lists = [_wrap(p) for p in df["Pathway"]]
    labels = ["\n".join(ll) for ll in line_lists]
    bands = np.array([len(ll) for ll in line_lists], dtype=float)  # vertical slots needed per row
    total_units = float(bands.sum())
    edges = np.concatenate([[0.0], np.cumsum(bands)])
    centers = (edges[:-1] + edges[1:]) / 2.0
    y = total_units - centers  # first (best) row at the top

    # Compact figure: roughly half the linear size of the original layout (width 8.9in -> 5.5in, and
    # the per-row height budget shrunk to match), with every fixed-inch margin below scaled the same
    # way so the proportions -- and the vertical centering of the legend -- still hold.
    MARGIN_SCALE = 0.60
    title_text = title or "Pathway Enrichment"
    if len(title_text) > 46:
        title_text = "\n".join(textwrap.wrap(title_text, 46)[:2])
    title_lines = title_text.count("\n") + 1
    height = max(2.3, 0.185 * total_units + 1.0) + 0.16 * (title_lines - 1)
    # Left margin from the MEASURED width of the longest (wrapped) term name at the tick-label
    # font, so long names (e.g. HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION) are never cut off
    # at the figure edge; the figure widens instead of squeezing the plot.
    from matplotlib.backends.backend_agg import FigureCanvasAgg as _FCA
    _mfig = matplotlib.figure.Figure()
    _r = _FCA(_mfig).get_renderer()
    _label_w_in = 0.0
    for _ll in line_lists:
        for _l in _ll:
            _t = _mfig.text(0, 0, _l, fontsize=7.6)
            _label_w_in = max(_label_w_in, _t.get_window_extent(renderer=_r).width / _mfig.dpi)
            _t.remove()
    _label_w_in += 0.2                      # tick marks + padding
    fig_w = max(5.5, _label_w_in + 3.05)    # keep >= ~3in for the plot + legend column
    fig = plt.figure(figsize=(fig_w, height))
    # Title is a figure-level suptitle with the plotting grid starting below it, so
    # it can never collide with the legend column headings on the right.
    top = 1.0 - (0.62 * MARGIN_SCALE + 0.15 * (title_lines - 1)) / height
    bottom = (0.62 * MARGIN_SCALE) / height
    left = max(0.25 * 5.5 / fig_w, _label_w_in / fig_w)
    gs = gridspec.GridSpec(1, 2, width_ratios=[5.6, 1.55], wspace=0.06,
                           figure=fig, left=left, right=0.975, top=top, bottom=bottom)
    ax = fig.add_subplot(gs[0, 0])
    ax_leg = fig.add_subplot(gs[0, 1])
    ax_leg.axis("off")
    ax_leg.set_xlim(0, 1)
    ax_leg.set_ylim(0, 1)

    counts = df[size_col].astype(float)
    cmin, cmax = counts.min(), counts.max()

    def _s(v):
        return 90.0 if cmax == cmin else 28.0 + 150.0 * (v - cmin) / (cmax - cmin)

    if color_by == "pvalue":
        cvals = df["P-value"].astype(float)
        vmin, vmax = 0.0, float(cvals.max()) or 1.0
        color_title = "P-value"
    else:
        cvals = df["neglog10_fdr"]
        vmin, vmax = float(cvals.min()), float(cvals.max())
        color_title = "−log10(FDR)"
    if np.isclose(vmin, vmax):
        vmin, vmax = max(0.0, vmin - 0.5), vmax + 0.5
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    cm = matplotlib.colormaps.get_cmap(cmap) if isinstance(cmap, str) else cmap

    for yi in y:
        ax.axhline(yi, color="#EEEEEE", linewidth=0.6, zorder=0)
    ax.scatter(df[x_col], y, s=counts.map(_s), c=cvals, cmap=cm, norm=norm,
               edgecolor="#333333", linewidth=0.5, zorder=3)
    if reference_line is not None:
        ax.axvline(reference_line, color="#808080", linestyle="--", linewidth=0.8, zorder=1)
    xvals = df[x_col].astype(float).values
    xmax = float(np.nanmax(xvals)) if n else 1.0
    if x_from_zero:
        upper = max(xmax, reference_line or 0) * 1.15 or 1.0
        ax.set_xlim(0, upper)
    else:
        xmin = float(np.nanmin(xvals))
        pad = 0.06 * (xmax - xmin) if xmax > xmin else 0.5
        ax.set_xlim(max(0.0, xmin - pad), xmax + pad)
        ax.grid(axis="x", color="#EEEEEE", linewidth=0.6, zorder=0)
    ax.set_ylim(-0.5, total_units + 0.5)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7.6, linespacing=0.95)
    ax.tick_params(axis="x", labelsize=8)
    ax.set_xlabel(x_label or x_col.replace("_", " "), fontsize=9.3)
    fig.suptitle(title_text, fontsize=10, fontweight="bold", linespacing=1.15,
                 x=0.5, y=1.0 - (0.18 * MARGIN_SCALE) / height, va="top")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # ---- Legend column: size key (top) then colour scale (bottom) ----
    # Rows are spaced by a FIXED inch amount (converted to this axis's fraction), so the legend stays
    # compact and legible whether the plot has 5 rows or 49, instead of stretching to fill the height.
    axis_height_in = max((top - bottom) * height, 0.1)

    def rows(inches):
        return inches / axis_height_in

    integer_sizes = bool(np.allclose(counts, np.round(counts)))
    reps = sorted(set(int(round(v)) for v in np.linspace(cmin, cmax, 3))) if integer_sizes else \
        sorted(set(float(f"{v:.2g}") for v in np.linspace(cmin, cmax, 3)))

    # Even, non-overlapping dot spacing: the gap between legend dots is set to the LARGEST legend
    # dot's own diameter (converted from points to inches) plus a small fixed pad, so dots of any
    # enrichment-ratio range sit with a small, even gap and never overlap.
    max_area_pts2 = max(_s(c) for c in reps)
    dot_diam_in = 2.0 * math.sqrt(max_area_pts2 / math.pi) / 72.0
    dot_step_in = dot_diam_in + 0.045
    dot_step = rows(dot_step_in)

    # Total legend block height (size key + colour scale), in the same inch-based units used for
    # spacing, so the whole block can be centred vertically in the legend column rather than
    # anchored at the top.
    size_title_lines = size_label if len(size_label) <= 15 else "\n".join(textwrap.wrap(size_label, 13))
    size_title_h = 0.12 * (size_title_lines.count("\n") + 1) + 0.06
    dots_h = dot_step_in * len(reps)
    rule_h = 0.04 + 0.13
    # Title line (~0.11in) + clearance for the top colour-bar tick label, which is centred on
    # the bar's top edge and so pokes ~half a text line above it.
    color_title_h = 0.27
    cbar_in = 0.55
    block_in = size_title_h + dots_h + rule_h + color_title_h + cbar_in
    y_cursor = min(1.0 - rows(0.06), 0.5 + rows(block_in) / 2.0)

    ax_leg.text(0.02, y_cursor, size_title_lines, fontsize=7, fontweight="bold", va="top", ha="left",
               linespacing=1.15)
    y_cursor -= rows(size_title_h)
    for c in reps:
        ax_leg.scatter([0.16], [y_cursor], s=_s(c), c="#9A9A9A", edgecolor="#333333", linewidth=0.5,
                       clip_on=False, zorder=3)
        label_txt = f"{c}" if integer_sizes else f"{c:g}"
        ax_leg.text(0.40, y_cursor, label_txt, fontsize=6.8, va="center", ha="left")
        y_cursor -= dot_step

    y_cursor -= rows(0.04)
    ax_leg.axhline(y_cursor, xmin=0.0, xmax=0.92, color="#DDDDDD", linewidth=0.6)
    y_cursor -= rows(0.16)

    ax_leg.text(0.02, y_cursor, color_title, fontsize=7.8, fontweight="bold", va="top", ha="left")
    y_cursor -= rows(color_title_h)
    cbar_h = rows(cbar_in)
    cax = ax_leg.inset_axes([0.06, max(y_cursor - cbar_h, 0.0), 0.22, cbar_h])
    sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cm)
    sm.set_array([])
    cb = fig.colorbar(sm, cax=cax)
    cb.ax.tick_params(labelsize=7.8, length=2.5, pad=2)
    cb.ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=4))
    cb.outline.set_linewidth(0.5)
    return fig


def enrichment_dotplot(results: pd.DataFrame, top_n: int = 25, title: str = None, sort_by: str = "P-value"):
    """
    Standard ORA dot plot: top 25 sets by raw P (best at the top), x = -log10(P), colour = raw P on a
    heat-colour gradient from 0 (red) to the largest plotted P, dot size = enrichment ratio (hits /
    expected). Long names wrap onto up to 3 lines rather than being truncated.
    """
    n_top = min(len(results), top_n or len(results))
    return pathway_dotplot(results, x_col="neglog10_p", x_label="-log10 (p-value)", top_n=top_n,
                           sort_by=sort_by, title=title or f"Overview of Enriched Metabolite Sets (Top {n_top})",
                           cmap=enrichment_heat_cmap(len(results), top_n or len(results)),
                           size_label="Enrichment Ratio", size_col="Enrichment_Ratio", color_by="pvalue",
                           x_from_zero=False)


# ===========================================================================
# 5. Streamlit UI
# ===========================================================================
def _csv(df: pd.DataFrame, index_col: str = None) -> bytes:
    return utils.to_download_bytes_csv(df.set_index(index_col) if index_col else df)


def _safe(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(s)).strip("_") or "Comparison"


TAB_MSEA = "Metabolite-Set Enrichment Analysis (MSEA)"
TAB_MGPA = "Metabolite–Gene Pathway Analysis (MGPA)"
MSEA_REFERENCE_LABEL_LIBRARY = "All Metabolites in Selected Library (default)"
MSEA_REFERENCE_LABEL_DETECTED = "All Detected/Measured Metabolites"


def _no_stats_warning(st):
    st.warning("Run a **two-group comparison** or an **ANOVA** in the Statistics tab first — this analysis starts "
               "from those existing results (nothing is recomputed here).")


def render_msea_tab(st, state, render_single_figure):
    """
    Tab "Metabolite-Set Enrichment Analysis (MSEA)" -- metabolite-set over-representation analysis.
    `state` is st.session_state; `render_single_figure` is app.py's shared centered-figure helper.
    """
    st.header(TAB_MSEA)
    st.caption("Metabolite Selection → ID Mapping → Reference Metabolome → Enrichment Analysis → Results.")
    sources = available_stats_sources(state)
    if not sources:
        _no_stats_warning(st)
        return
    ctx = _render_metabolite_selection(st, sources, METHOD_A, kp="msea")
    if ctx is not None:
        _render_msea(st, state, render_single_figure, ctx)


def render_mgpa_tab(st, state, render_single_figure):
    """Tab "Metabolite–Gene Pathway Analysis (MGPA)" -- metabolite -> HMDB -> gene -> gene-set ORA."""
    st.header(TAB_MGPA)
    st.caption("Metabolite Selection → HMDB Mapping → Gene Mapping → Gene-Set ORA → Results → Dot Plot.")
    sources = available_stats_sources(state)
    if not sources:
        _no_stats_warning(st)
        return
    ctx = _render_metabolite_selection(st, sources, METHOD_B)
    if ctx is not None:
        _render_option_b(st, state, render_single_figure, ctx)


def _render_metabolite_selection(st, sources, method, kp="pa"):
    """
    Section 1 (both tabs): selection from the existing Statistics-tab results only. Returns a context
    dict, or None when nothing is selected. `kp` prefixes the widget keys (the MGPA tab keeps the
    original 'pa' keys; the MSEA tab uses 'msea' so both tabs can render in the same run).
    """
    st.subheader("Input Metabolite Selection")
    c1, c2, c3, c4 = st.columns([1.6, 1.4, 1.0, 1.2])
    src_label = c1.selectbox("Statistical results to use", list(sources.keys()), key=f"{kp}_source",
                             help="Only results already computed in the Statistics tab are listed.")
    src = sources[src_label]
    stats_df, p_col, fdr_col = src["df"], src["p_col"], src["fdr_col"]
    mode = c2.radio("Select metabolites", [MODE_P, MODE_FDR, MODE_ALL], key=f"{kp}_mode",
                    format_func=lambda m: m.replace("FDR", fdr_col) if m == MODE_FDR else m)
    if mode == MODE_P:
        thr = c3.number_input("P-value ≤", value=0.05, min_value=0.0001, max_value=1.0, step=0.01,
                              format="%.4f", key=f"{kp}_thr_p")
    elif mode == MODE_FDR:
        thr = c3.number_input(f"{fdr_col} ≤", value=0.25, min_value=0.0001, max_value=1.0, step=0.01,
                              format="%.4f", key=f"{kp}_thr_fdr")
    else:
        thr = None
        c3.caption("No threshold — every tested metabolite is used.")
    direction = "both"
    if src["kind"] == "two_group" and src["fc_col"]:
        g_a, g_b = src["groups"]
        dir_label = c4.selectbox("Direction", ["Both directions", f"Up in {g_b} (Log2FC > 0)",
                                               f"Down in {g_b} (Log2FC < 0)"], key=f"{kp}_direction")
        direction = "both" if dir_label.startswith("Both") else ("up" if dir_label.startswith("Up") else "down")
    else:
        c4.caption("ANOVA is non-directional — no up/down filter.")

    selected = select_metabolites(stats_df, p_col, fdr_col, mode, thr if thr is not None else 1.0,
                                  direction, src["fc_col"])
    m1, m2 = st.columns(2)
    m1.metric("Metabolites tested (Statistics tab)", len(stats_df))
    m2.metric("Selected as ORA input", len(selected))
    crit = ("all tested metabolites" if mode == MODE_ALL else
            f"{p_col if mode == MODE_P else fdr_col} ≤ {thr:g}")
    st.caption(f"Source: **{src_label}** · criterion: **{crit}**"
               + ("" if direction == "both" else f" · direction: **{direction}**") + ".")
    with st.expander(f"Selected metabolites ({len(selected)})", expanded=False):
        st.dataframe(utils.format_df_for_display(selected), width="stretch", height=260)
    if len(selected) == 0:
        st.warning("No metabolites pass this criterion — relax the threshold or choose another option.")
        return None
    return {"src_label": src_label, "src": src, "stats_df": stats_df, "p_col": p_col, "fdr_col": fdr_col,
            "mode": mode, "thr": thr, "direction": direction, "selected": selected}


# ---------------------------------------------------------------------------
# MSEA tab (sections 2-5)
# ---------------------------------------------------------------------------
MSEA_ID_COLUMNS = ["Metabolite", "Matched_Name", "HMDB_ID", "SMPDB_ID", "KEGG_Map_ID", "N_SMPDB_Pathways",
                   "N_KEGG_Maps", "KEGG_ID", "ChEBI_ID", "PubChem_CID", "Match_Method", "Matched_On"]
_EXCEL_ID_CAP = 2500  # SMPDB IDs per cell in the Excel workbook (Excel's 32,767-character cell limit)



def pca_module_palettes():
    from metabolomics_modules import pca_module
    return pca_module.PALETTES.keys()


# ---------------------------------------------------------------------------
# 5c. Metabolite-level (feature-space) UMAP -- one dot per METABOLITE (not per
# sample), positioned by the similarity of its expression profile across all
# samples, colored by which top enriched pathway it belongs to. This is the
# metabolomics analogue of a single-cell "feature UMAP": lots of points,
# clustered by co-expression pattern, colored by group/pathway identity --
# rather than the sample-space ordination above, which has one dot per sample.
# ---------------------------------------------------------------------------
def assign_metabolite_pathway_labels(bg_idmap: pd.DataFrame, in_lib: pd.Series, library,
                                     results: pd.DataFrame, top_n_pathways: int = 8,
                                     rank_by: str = "FDR", min_members: int = 3):
    """
    For every background metabolite, assign it a single pathway label: the
    highest-ranked (by `rank_by`) of the top `top_n_pathways` enriched pathways
    it is a detected member of; metabolites in none of them are labeled
    'Other'. A metabolite in more than one shown pathway is colored by
    whichever of those pathways ranks best (priority = enrichment rank, not
    alphabetical) -- it still counts toward that pathway's dot cluster only.
    Pathways with fewer than `min_members` detected members (after this
    priority assignment) are dropped back into 'Other' so tiny, unreadable
    clusters don't fragment the legend.
    Returns a Series of labels indexed by metabolite NAME (bg_idmap["Metabolite"]),
    so it lines up directly with a display-named log2 matrix.
    """
    labels = pd.Series("Other", index=bg_idmap.index)
    if results is None or len(results) == 0:
        return labels
    ranked = results.sort_values(rank_by, ascending=True).head(top_n_pathways)
    for _, row in ranked.iterrows():
        pname = row["Pathway"]
        key_set = library.sets.get(pname, set())
        if not key_set:
            continue
        mask = in_lib.isin(key_set) & (in_lib != "") & (labels == "Other")
        labels[mask] = pname
    counts = labels.value_counts()
    too_small = counts[(counts.index != "Other") & (counts < min_members)].index
    labels[labels.isin(too_small)] = "Other"
    labels.index = bg_idmap["Metabolite"].values
    return labels


def metabolite_expression_umap(log2_df: pd.DataFrame, labels: pd.Series, palette: str = "Default (tab10)",
                               n_neighbors: int = 15, min_dist: float = 0.3, metric: str = "correlation"):
    """
    UMAP embedding of METABOLITES (rows of log2_df), each represented by its
    z-scored expression profile across all samples (columns). Returns
    (fig, emb_df [metabolite x UMAP1/UMAP2]). Points are colored by `labels`
    (from assign_metabolite_pathway_labels); 'Other' is drawn first, small,
    pale gray, so the highlighted pathway clusters sit on top and stand out --
    the same convention single-cell UMAPs use for an "all other cells" layer.
    """
    from metabolomics_modules import pca_module
    import umap
    z = log2_df.sub(log2_df.mean(axis=1), axis=0).div(log2_df.std(axis=1).replace(0, np.nan), axis=0)
    z = z.fillna(0.0)
    n_neighbors = int(min(n_neighbors, max(2, z.shape[0] - 1)))
    reducer = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist, metric=metric,
                        n_components=2, random_state=42)
    emb = reducer.fit_transform(z.values)
    emb_df = pd.DataFrame(emb, index=z.index, columns=["UMAP1", "UMAP2"])

    fig, ax = plt.subplots(figsize=(6.6, 5.6))
    lab = labels.reindex(emb_df.index).fillna("Other")
    pathway_cats = [c for c in lab.value_counts().index.tolist() if c != "Other"]
    colors = pca_module.PALETTES.get(palette, palette) if isinstance(palette, str) else palette

    if "Other" in lab.values:
        m = lab == "Other"
        ax.scatter(emb_df.loc[m, "UMAP1"], emb_df.loc[m, "UMAP2"], s=14, alpha=0.35,
                  color="#d1d5db", edgecolor="none", label=f"Other (n={m.sum()})", zorder=1)
    for i, cat in enumerate(pathway_cats):
        m = lab == cat
        color = colors[i % len(colors)]
        ax.scatter(emb_df.loc[m, "UMAP1"], emb_df.loc[m, "UMAP2"], s=34, alpha=0.9,
                  color=color, edgecolor="k", linewidth=0.3, label=f"{cat} (n={m.sum()})", zorder=2)

    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title("Metabolite-Level UMAP\n(colored by top enriched pathway)",
                 fontsize=10.5, fontweight="bold", loc="center")
    ax.legend(fontsize=7, markerscale=1.1, loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False)
    fig.tight_layout()
    return fig, emb_df


def _render_metabolite_umap(st, state, render_single_figure, results, bg_idmap, in_lib, library, tag):
    st.subheader("Metabolite-Level UMAP (Feature Space)")
    st.caption(
        "A UMAP-style view like a single-cell feature UMAP, but here every dot is a **metabolite**, not a "
        "sample: each is placed by how similar its expression pattern is across all samples (z-scored log2 "
        "abundance), so metabolites that rise and fall together cluster together. Dots are colored by which "
        "top enriched pathway they belong to; metabolites in none of the shown pathways are pale gray "
        "('Other'). This is a feature-space view -- it complements, and is independent from, the "
        "sample-space PCA/UMAP above."
    )
    log2_data = state.get("log2_data")
    if log2_data is None:
        st.info("Load and normalize data (Data Upload -> Normalization) to compute the metabolite-level UMAP.")
        return
    log2_df = utils.with_display_feature_names(log2_data)
    log2_df = log2_df.loc[[m for m in bg_idmap["Metabolite"] if m in log2_df.index]]
    if log2_df.shape[0] < 15:
        st.info(f"Only {log2_df.shape[0]} background metabolites are available -- too few for a stable "
               "metabolite-level UMAP (need at least ~15).")
        return

    u1, u2, u3, u4 = st.columns(4)
    rank_by = u1.selectbox("Rank pathways by", ["FDR", "P-value", "Holm"], key="mumap_rank")
    max_pw = min(15, len(results)) if results is not None else 0
    top_n_pathways = u2.number_input("Pathways to highlight", min_value=1, max_value=max(1, max_pw),
                                     value=min(8, max(1, max_pw)), step=1, key="mumap_topn")
    min_members = u3.number_input("Min metabolites/pathway", min_value=2, max_value=20, value=3, step=1,
                                  key="mumap_minmembers")
    palette_choice = u4.selectbox("Palette", list(pca_module_palettes()), key="mumap_palette")

    labels = assign_metabolite_pathway_labels(bg_idmap, in_lib, library, results,
                                              top_n_pathways=int(top_n_pathways), rank_by=rank_by,
                                              min_members=int(min_members))
    labels = labels.reindex(log2_df.index).fillna("Other")
    n_highlighted = int((labels != "Other").sum())
    st.caption(f"**{log2_df.shape[0]} metabolites** plotted -- **{n_highlighted}** colored across "
              f"**{labels[labels != 'Other'].nunique()}** highlighted pathways, "
              f"**{log2_df.shape[0] - n_highlighted}** shown as 'Other'.")
    if n_highlighted == 0:
        st.warning("None of the ranked pathways have enough detected members to highlight -- try increasing "
                  "'Pathways to highlight' or lowering 'Min metabolites/pathway'.")
        return

    n1, n2, n3 = st.columns(3)
    n_neighbors = n1.slider("n_neighbors", 2, min(50, max(2, log2_df.shape[0] - 1)),
                            min(15, max(2, log2_df.shape[0] - 1)), 1, key="mumap_neighbors")
    min_dist = n2.slider("min_dist", 0.0, 0.99, 0.3, 0.05, key="mumap_min_dist")
    metric = n3.selectbox("Distance metric", ["correlation", "euclidean", "cosine"], key="mumap_metric")

    try:
        with st.spinner("Computing metabolite-level UMAP..."):
            fig, emb_df = metabolite_expression_umap(log2_df, labels, palette=palette_choice,
                                                      n_neighbors=int(n_neighbors), min_dist=min_dist,
                                                      metric=metric)
        render_single_figure(fig, caption="Each dot is one metabolite; color = top enriched pathway membership.",
                            download_name=f"{tag}_MetaboliteUMAP", key_prefix="msea_metabolite_umap",
                            width_ratio=(1, 4, 1))
    except ImportError:
        st.error("The `umap-learn` package isn't available in this environment -- the metabolite-level UMAP "
                "can't be computed.")
        return

    dl = emb_df.copy()
    dl["Pathway_Label"] = labels.reindex(dl.index)
    st.download_button(f"Download {tag}_Metabolite_UMAP_Coordinates.csv", utils.to_download_bytes_csv(dl),
                       file_name=f"{tag}_Metabolite_UMAP_Coordinates.csv", mime="text/csv",
                       key="msea_dl_metabolite_umap")


def _render_msea(st, state, render_single_figure, ctx):
    src, stats_df, selected, mode = ctx["src"], ctx["stats_df"], ctx["selected"], ctx["mode"]

    # ---------------- 2: ID mapping (HMDB_ID + SMPDB_ID + KEGG_Map_ID) ----------------
    st.subheader("ID Mapping")
    row_annotations = state.get("row_annotations")
    id_cols = libs.detected_id_columns(row_annotations)
    found = [f"`{c}` ({k.upper() if k != 'pubchem' else 'PubChem'})" for k, c in id_cols.items() if c]
    st.caption(f"HMDB_ID matched via {'row annotations (' + ', '.join(found) + ')' if found else 'name/synonym'} "
               f"→ SMPDB_ID / KEGG_Map_ID lookup (up to {libs.ID_LIST_DISPLAY_MAX} shown per row; CSV holds all).")
    bg_idmap = standardize_for_option_a(stats_df.index, row_annotations)
    sel_idmap = bg_idmap.set_index("Metabolite").loc[[str(i) for i in selected.index]].reset_index()
    sel_view = libs.add_hmdb_pathway_ids(sel_idmap)[MSEA_ID_COLUMNS]
    k1, k2, k3, k4, k5 = st.columns(5)
    k1.metric("Selected", len(sel_view))
    k2.metric("With HMDB_ID", int((sel_view["HMDB_ID"] != "").sum()))
    k3.metric("With SMPDB_ID", int((sel_view["N_SMPDB_Pathways"] > 0).sum()))
    k4.metric("With KEGG_Map_ID", int((sel_view["N_KEGG_Maps"] > 0).sum()))
    k5.metric("Unmapped", int((sel_view["Match_Method"] == libs.MATCH_UNMAPPED).sum()))
    st.markdown("**ID Mapping table**")
    st.dataframe(sel_view, width="stretch", height=min(380, 38 + 35 * max(1, len(sel_view))), hide_index=True)

    # ---------------- 3: reference metabolome (all detected/measured metabolites) ----------------
    st.subheader("Reference Metabolome")
    lc1, lc2 = st.columns([1.3, 2.2])
    lib_choice = lc1.selectbox("Metabolite-set library", libs.BUILTIN_LIBRARIES + [libs.LIB_CUSTOM],
                               key="msea_library")
    library = None
    if lib_choice == libs.LIB_CUSTOM:
        up = lc2.file_uploader("Custom metabolite sets (CSV: Set_Name, Metabolite — long or ';'-separated; "
                               "or .gmt)", type=["csv", "txt", "gmt", "tsv"], key="msea_custom_lib")
        if up is None:
            lc2.info("Upload a metabolite-set file to use a custom library.")
        else:
            try:
                library = libs.parse_custom_library(up.getvalue().decode("utf-8", errors="ignore"), up.name)
                lc2.caption(f"Parsed **{len(library)}** sets / **{len(library.universe)}** unique members "
                            f"from `{up.name}`.")
            except ValueError as e:
                lc2.error(f"Could not parse that file: {e}")
    else:
        library = libs.get_library(lib_choice)
        lc2.caption(f"**{library.name}** — {len(library)} sets, {len(library.universe)} unique members.")
    if library is None:
        return
    if len(library) == 0:
        st.error("The selected library contains no sets.")
        return
    in_lib = library_membership(bg_idmap, library)
    detected = {k for k in in_lib if k}
    sel_in_lib = library_membership(sel_idmap, library)
    n_sets_detected = sum(1 for s in library.sets.values() if s & detected)
    restrict_detected = st.checkbox(
        "Restrict reference to detected/measured metabolites only", value=False, key="msea_restrict_detected",
        help="Off (default): background = whole library. On: background = only metabolites detected in this "
             "dataset (more conservative).")
    reference = REF_DETECTED if restrict_detected else REF_LIBRARY
    ref_label = MSEA_REFERENCE_LABEL_DETECTED if restrict_detected else MSEA_REFERENCE_LABEL_LIBRARY
    st.caption(f"Background / reference: **{ref_label}** ({library.name}).")
    r1, r2, r3, r4 = st.columns(4)
    r1.metric("Detected/measured metabolites", len(bg_idmap))
    r2.metric("Detected in library", len(detected))
    r3.metric("Selected in library (n)", len({k for k in sel_in_lib if k}))
    r4.metric("Library sets with ≥ 1 detected member", f"{n_sets_detected} / {len(library)}")
    sel_view = sel_view.assign(In_Selected_Library=np.where(sel_in_lib.values != "", "Yes", "No"))
    not_in_lib = sel_view[sel_view["In_Selected_Library"] == "No"]
    if len(not_in_lib):
        st.warning(f"**{len(not_in_lib)} of {len(sel_view)} selected metabolites are not in the {library.name} "
                   "library** and cannot contribute to this enrichment.")
        with st.expander(f"Selected metabolites not in the library ({len(not_in_lib)})", expanded=False):
            st.dataframe(not_in_lib[["Metabolite", "HMDB_ID", "KEGG_ID", "Match_Method"]], width="stretch",
                         height=min(300, 38 + 35 * len(not_in_lib)), hide_index=True)
    labels = getattr(library, "member_labels", {}) or {}
    key_to_name = {}
    for k, nm in zip(in_lib, bg_idmap["Metabolite"]):
        if k:
            key_to_name.setdefault(k, str(nm))
    universe_keys = detected if restrict_detected else set(library.universe)
    ref_df = pd.DataFrame({"Library_Key": sorted(universe_keys)})
    ref_df["Metabolite"] = ref_df["Library_Key"].map(lambda k: key_to_name.get(k) or labels.get(k, k))
    ref_df["Library_Member_Name"] = ref_df["Library_Key"].map(lambda k: labels.get(k, ""))
    ref_df["Detected_In_Dataset"] = np.where(ref_df["Library_Key"].isin(detected), "Yes", "No")
    with st.expander(f"Reference metabolome used ({len(ref_df)} metabolites)", expanded=False):
        st.dataframe(ref_df, width="stretch", height=min(320, 38 + 35 * max(1, len(ref_df))), hide_index=True)

    # ---------------- 4: enrichment analysis ----------------
    st.subheader("Enrichment Analysis")
    s1, s2, s3 = st.columns(3)
    test = s1.selectbox("Statistical test", [TEST_HYPERGEOM, TEST_FISHER], key="msea_test",
                        help="Both give the same P-value for this 2×2 enrichment framing.")
    min_size = s2.number_input("Min set size (in background)", min_value=1, max_value=500, value=2, step=1,
                               key="msea_min")
    max_size_in = s3.number_input("Max set size (0 = no limit)", min_value=0, max_value=5000, value=0, step=5,
                                  key="msea_max")
    max_size = int(max_size_in) or None
    results, info = run_msea_ora(sel_idmap, bg_idmap, library, reference, None, int(min_size), max_size, test,
                                 ADJUST_ALL)
    st.caption("Correction: **Holm** and **FDR** across all sets in range (untested sets count as P = 1).")
    if mode == MODE_ALL and restrict_detected:
        st.info("Selecting all metabolites with the detected-only reference makes every set P = 1 by "
                "construction — use a significance-based selection instead.")
    e1, e2, e3, e4, e5 = st.columns(5)
    e1.metric("Background size (N)", info["N"])
    e2.metric("Input in background (n)", info["n"])
    e3.metric("Sets with ≥ 1 hit", info["n_sets_tested"])
    e4.metric("Sets in correction (m)", info["n_adjusted"])
    e5.metric("FDR ≤ 0.05", int((results["FDR"] <= 0.05).sum()) if len(results) else 0)
    universe_desc = (f"{info['N']} metabolites detected/tested in this dataset that are members of the "
                     f"{library.name} library" if restrict_detected else
                     f"{info['N']} metabolites in the {library.name} library (default background)")
    st.caption(f"N = {universe_desc}. n = selected metabolites in the universe. K = set members within the "
               f"universe; k = input in the set; Expected = n·K/N; P = P(X ≥ k) ({test.lower()}); Holm and FDR "
               f"across {info['n_adjusted']} sets (size {int(min_size)}–{max_size or '∞'}).")

    # ---------------- 5: results ----------------
    st.subheader("Results")
    st.markdown(f"**Pathway Results — Library: {library.name} · Reference: {ref_label}**")
    if len(results) == 0:
        st.warning("No metabolite set had an overlapping input within the size limits — nothing to report.")
    else:
        st.dataframe(utils.format_df_for_display(results), width="stretch", height=380, hide_index=True)
        st.caption("Enrichment_Ratio = (k/n)/(K/N) = k/Expected. Sorted by P-value.")

    fig = None
    tag = f"{_safe(src['comparison'])}_MSEA_{_safe(library.name)}"
    if len(results):
        st.markdown("**Enrichment Dot Plot**")
        d1, d2 = st.columns(2)
        top_choice = d1.selectbox("Sets shown", ["Top 10", "Top 20", "Top 25", "Top 30", "All"], index=2,
                                  key="msea_top")
        top_n = None if top_choice == "All" else int(top_choice.split()[1])
        sort_by = d2.selectbox("Sort by", list(SORT_OPTIONS.keys()), key="msea_sort")
        shown = min(len(results), top_n or len(results))
        fig = enrichment_dotplot(results, top_n=top_n, sort_by=sort_by)
        render_single_figure(fig, caption=f"{shown} of {len(results)} sets with hits ({library.name}), sorted by "
                                          f"{sort_by}.", download_name=f"{tag}_DotPlot", key_prefix="msea_dotplot",
                             width_ratio=(1, 3, 1))

    # ---------------- 5c: metabolite-level (feature-space) UMAP ----------------
    if len(results):
        _render_metabolite_umap(st, state, render_single_figure, results, bg_idmap, in_lib, library, tag)

    # ---------------- downloads ----------------
    st.markdown("**Downloads**")
    # "ID mapping (CSV)" exports exactly the Section 2 table above (same columns, same capped ID lists) --
    # the fuller, uncapped version with library membership lives only in the Excel workbook, to avoid a CSV
    # that looks like the table but silently differs.
    xls_view = libs.add_hmdb_pathway_ids(sel_idmap, cap=_EXCEL_ID_CAP)[MSEA_ID_COLUMNS].assign(
        In_Selected_Library=sel_view["In_Selected_Library"].values)
    params = {
        "Statistics source": ctx["src_label"], "Selection": mode, "Threshold": "" if ctx["thr"] is None else ctx["thr"],
        "Direction": ctx["direction"], "P-value column": ctx["p_col"], "FDR column (input stats)": ctx["fdr_col"],
        "Method": TAB_MSEA, "Library": library.name, "Library sets": len(library),
        "Library unique members": len(library.universe),
        "Reference metabolome": ref_label,
        "Reference definition": universe_desc,
        "Background size (N)": info["N"], "Input size in background (n)": info["n"],
        "Statistical test": test,
        "Multiple-testing correction": f"Holm and Benjamini-Hochberg across {info['n_adjusted']} sets (all sets "
                                       "in the selected library within the size limits)",
        "Min set size": int(min_size), "Max set size": max_size or "none",
        "Sets with >= 1 hit": info["n_sets_tested"], "Sets skipped (size filter)": info["n_sets_size_filtered"],
        "Sets with no overlap (P = 1, counted in the correction)": info["n_sets_no_overlap"],
        "Selected metabolites": len(sel_view), "Selected in library": info["n_selected_in_library"],
        "ID mapping source": f"HMDB_ID: row annotations / Omnix Metabolomics embedded reference ({len(ref.METABOLITES)} "
                             f"metabolites); SMPDB_ID + KEGG_Map_ID: {libs.HMDB_PATHWAY_SOURCE}",
        "Library sources": f"KEGG = {libs.KEGG_SOURCE}; SMPDB from the {libs.HMDB_PATHWAY_SOURCE}; LIPID MAPS "
                           "curated",
        "Dot plot": "x = -log10(P-value), colour = P-value, size = enrichment ratio "
                   "(PNG/TIFF/SVG/PDF download beneath the figure)",
    }
    params_df = pd.DataFrame({"Value": [str(v) for v in params.values()]},
                             index=pd.Index(list(params.keys()), name="Parameter"))
    b = st.columns([1, 1, 1, 3])
    b[0].download_button("ID mapping (CSV)", _csv(sel_view[MSEA_ID_COLUMNS], "Metabolite"),
                         file_name=f"{tag}_ID_Mapping.csv", mime="text/csv", key="msea_dl_map")
    b[1].download_button("Pathway results (CSV)", _csv(results, "Pathway_ID"), file_name=f"{tag}_Pathway_Results.csv",
                         mime="text/csv", key="msea_dl_res", disabled=len(results) == 0)
    sheets = {"Parameters": params_df, "Input_Selection": selected,
              "ID_Mapping": xls_view.set_index("Metabolite"),
              "Not_In_Library": xls_view[xls_view["In_Selected_Library"] == "No"].set_index("Metabolite"),
              "Reference_Metabolome": ref_df.set_index("Library_Key"),
              "Pathway_Results": results.set_index("Pathway_ID") if len(results) else results}
    b[2].download_button("Complete results (Excel)", utils.to_download_bytes_xlsx(sheets),
                         file_name=f"{tag}_Complete_Results.xlsx",
                         mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", key="msea_dl_xlsx")
    if fig is not None:
        plt.close(fig)


# ---------------------------------------------------------------------------
# MGPA tab -- Metabolite -> HMDB -> Gene -> Unique Gene List -> gene-set ORA
# (Section 4 onward mirrors ProteoAI Pro's "Over-Representation Analysis (ORA)"; engine and dot plot:
#  gene_set_enrichment.run_ora_hypergeometric / ora_dotplot, ported from ProteoAI's enrichment.py)
# ---------------------------------------------------------------------------
def _render_option_b(st, state, render_single_figure, ctx):
    method = METHOD_B
    src_label, src, stats_df = ctx["src_label"], ctx["src"], ctx["stats_df"]
    p_col, fdr_col, mode, thr, direction = ctx["p_col"], ctx["fdr_col"], ctx["mode"], ctx["thr"], ctx["direction"]
    selected = ctx["selected"]

    # ---------------- Step 2: HMDB mapping ----------------
    st.subheader("HMDB Mapping")
    row_annotations = state.get("row_annotations")
    hmdb_col = find_hmdb_column(row_annotations)
    hsum = hgm.summary()
    st.caption(f"HMDB_ID matched via {'row annotations (`' + hmdb_col + '`)' if hmdb_col else 'name/synonym'}. "
               f"Genes are HMDB protein associations ({hsum['associations']:,} associations, "
               f"{hsum['genes']:,} genes).")
    background_map = map_metabolites_to_hmdb(stats_df.index, row_annotations)
    sel_map = background_map.set_index("Metabolite").loc[[str(i) for i in selected.index]].reset_index()
    n_mapped = int((sel_map["Mapping_Status"] == STATUS_MAPPED).sum())
    n_hmdb_only = int((sel_map["Mapping_Status"] == STATUS_NO_GENES).sum())
    n_unmapped = int((sel_map["Mapping_Status"] == STATUS_UNMAPPED).sum())
    k1, k2, k3, k4 = st.columns(4)
    k1.metric("Selected", len(sel_map))
    k2.metric("Mapped (HMDB + genes)", n_mapped)
    k3.metric("HMDB ID, no HMDB genes", n_hmdb_only)
    k4.metric("Unmapped", n_unmapped)
    mapped_view = sel_map[sel_map["Mapping_Status"] != STATUS_UNMAPPED]
    st.dataframe(mapped_view[["Metabolite", "HMDB_ID", "HMDB_Name", "N_HMDB_Proteins", "N_HMDB_Genes",
                              "Mapping_Source", "Mapping_Status", "Note"]],
                 width="stretch", height=min(360, 38 + 35 * max(1, len(mapped_view))), hide_index=True)
    not_used = sel_map[sel_map["Mapping_Status"] != STATUS_MAPPED]
    if len(not_used):
        st.warning(f"**{len(not_used)} of {len(sel_map)} selected metabolites cannot contribute genes** "
                   f"({n_unmapped} with no HMDB ID, {n_hmdb_only} whose HMDB record lists no protein "
                   "associations). They are excluded from the gene list and listed below.")
        st.dataframe(not_used[["Metabolite", "HMDB_ID", "Mapping_Status"]], width="stretch",
                     height=min(250, 38 + 35 * len(not_used)), hide_index=True)
    if n_mapped == 0:
        st.error("None of the selected metabolites have HMDB gene associations — gene-based ORA cannot run. "
                 "Add an `HMDB_ID` column to your row-annotation file (Tab 1) or use recognizable names.")
        return

    # Step 3 is rendered into a placeholder so it appears above the ORA settings, while the
    # Unique Gene List's In_Gene_Set_Collection flag reflects the collection chosen in step 4.
    step3 = st.container()

    # ---------------- Step 4: gene-set ORA settings (mirrors ProteoAI Pro's ORA) ----------------
    st.subheader("Gene-Set ORA Settings")
    st.caption("Hypergeometric test against a chosen gene-set library. **MSigDB (Hallmark)** runs offline; "
               "other databases are fetched live from Enrichr and require internet access.")
    c1, c2 = st.columns(2)
    organism = c1.selectbox("Organism", list(gse.ORGANISM_CATALOG), key="mgpa_ora_organism")
    org_cfg = gse.ORGANISM_CATALOG[organism]
    db_options = gse.ORA_GSEA_DATABASE_COLLECTIONS
    db_labels = {c: c + ("" if org_cfg["databases"][c]["available"] else " — unavailable") for c in db_options}
    collection = c2.selectbox("Gene-set database", db_options, key="mgpa_ora_collection",
                              index=db_options.index(gse.DEFAULT_ORA_COLLECTION), format_func=lambda c: db_labels[c])
    db_info = org_cfg["databases"][collection]
    gene_sets, resolved_library, resolved_version = None, db_info.get("library"), db_info.get("version")
    if not db_info["available"]:
        st.warning(f"**{collection}** is not available for {organism}: {db_info.get('note', '')}")
    else:
        st.caption(f"Database version: **{db_info['version']}**  ·  "
                   f"Gene identifier type: **{org_cfg['gene_id_type']}**")
        try:
            with st.spinner(f"Fetching {collection} for {organism}..."):
                gene_sets, resolved_library, resolved_version = gse.fetch_gene_set_library(organism, collection)
        except gse.EnrichmentError as e:
            st.error(str(e))
    collection_label = collection
    if gene_sets is None:
        # Section 3's In_Gene_Set_Collection flags need a collection: fall back to the bundled Hallmark set
        # (the ORA itself stays disabled until the selected database can be loaded).
        st.warning(f"Section 3's collection flags use the bundled **{gse.DEFAULT_ORA_COLLECTION}** meanwhile; "
                   "the ORA can run once the selected database is available.")
        collection_label = gse.DEFAULT_ORA_COLLECTION
        collection_genes = set().union(*[{str(g).upper() for g in v}
                                         for v in gse.fetch_gene_set_library(organism, collection_label)[0].values()])
    else:
        collection_genes = set().union(*[{str(g).upper() for g in v} for v in gene_sets.values()])

    universe_choice = st.radio(
        "Background / universe", gse.ORA_UNIVERSE_OPTIONS, key="mgpa_ora_universe",
        help="Genes considered testable. Detected-metabolite genes (default) reflect what could actually have "
             "been observed in this experiment.")
    c3, c4 = st.columns(2)
    min_size_ora = c3.number_input("Min gene set size", min_value=1, value=1, step=1, key="mgpa_ora_min_size")
    max_size_ora = c4.number_input("Max gene set size (0 = no limit)", min_value=0, value=0, step=10,
                                   key="mgpa_ora_max_size")

    # ---------------- Step 3 (filled now): Metabolite–Gene mapping + Unique Gene List ----------------
    mg_df, ug_df, no_gene = build_metabolite_gene_tables(sel_map, collection_genes)
    with step3:
        st.subheader("Metabolite–Gene Mapping & Unique Gene List")
        n_in_coll = int((ug_df["In_Gene_Set_Collection"] == "Yes").sum()) if len(ug_df) else 0
        m1, m2, m3, m4 = st.columns(4)
        m1.metric("Metabolites contributing genes", mg_df["Metabolite"].nunique() if len(mg_df) else 0)
        m2.metric("Metabolite–gene associations", len(mg_df))
        m3.metric("Unique genes", len(ug_df))
        m4.metric(f"In {collection_label[:28]}", n_in_coll)
        st.markdown("**Metabolite–Gene Mapping**")
        st.dataframe(mg_df, width="stretch", height=260, hide_index=True)
        n_other = int((mg_df["Gene_Label"] == GENE_LABEL_OTHER).sum()) if len(mg_df) else 0
        if n_other:
            st.caption(f"{n_other} association(s) carry a non-gene HMDB label — kept above, excluded from the "
                       "Unique Gene List / ORA.")
        st.markdown("**Unique Gene List**")
        st.caption(f"{len(ug_df)} unique genes — input list for the ORA below. "
                   f"`In_Gene_Set_Collection = No`: {len(ug_df) - n_in_coll} gene(s) absent from {collection_label}.")
        st.dataframe(ug_df, width="stretch", height=260, hide_index=True)
        if no_gene:
            st.warning(f"{len(no_gene)} mapped metabolite(s) have no gene-symbol association in HMDB "
                       f"({', '.join(no_gene[:10])}{'…' if len(no_gene) > 10 else ''}).")
    if len(ug_df) == 0:
        st.error("The selected metabolites have no associated genes — gene-set ORA cannot run.")
        return

    comparison_genes = list(dict.fromkeys(ug_df["Symbol_Used_For_ORA"]))
    detected_genes = sorted(genes_for_mapping(background_map, collection_genes))
    signature = (tuple(comparison_genes), organism, collection, universe_choice, int(min_size_ora),
                 int(max_size_ora))
    if st.button("Run Over-Representation Analysis", type="primary", icon="🧬", key="mgpa_run_ora",
                 disabled=gene_sets is None):
        universe = detected_genes if universe_choice == gse.UNIVERSE_DETECTED_LABEL else None
        try:
            with st.spinner(f"Running hypergeometric ORA against {len(gene_sets)} gene sets..."):
                result, run_meta = gse.run_ora_hypergeometric(
                    comparison_genes, gene_sets, universe_genes=universe, min_set_size=int(min_size_ora),
                    max_set_size=(int(max_size_ora) if max_size_ora > 0 else None))
            run_meta.update({"Organism": organism, "Gene identifier type": org_cfg["gene_id_type"],
                             "Database name": resolved_library, "Database version": resolved_version})
            state["mgpa_ora_result"] = result
            state["mgpa_ora_method_used"] = f"Over-Representation Analysis ({organism}, {collection})"
            state["mgpa_ora_run_meta"] = run_meta
            state["mgpa_ora_universe_genes"] = universe
            state["mgpa_ora_signature"] = signature
            state["mgpa_dotplot_on"] = False
            if len(result) == 0:
                st.warning("No gene sets could be tested — none of the comparison genes were found in the chosen "
                           "universe. See the metadata panel below for how many genes were mapped.")
            else:
                st.success(f"{len(result)} gene sets tested.")
        except gse.EnrichmentError as e:
            st.error(str(e))

    result = state.get("mgpa_ora_result")
    fig = None
    if result is None:
        st.info("Click **Run Over-Representation Analysis** to test the Unique Gene List (Section 3).")
    else:
        # ---------------- Step 5: results ----------------
        if state.get("mgpa_ora_signature") != signature:
            st.warning("The input selection or ORA settings changed since the last run — click **Run "
                       "Over-Representation Analysis** again to update the results below.")
        ora_tag = f"{_safe(src['comparison'])}_ORA_{_safe(collection)}"
        st.subheader(f"Results — {state.get('mgpa_ora_method_used')}")
        st.dataframe(utils.format_df_for_display(result), width="stretch", height=350)
        st.download_button(f"Download {ora_tag}_Results.csv", utils.to_download_bytes_csv(result),
                           file_name=f"{ora_tag}_Results.csv", mime="text/csv", key="mgpa_dl_ora_results")
        run_meta = state.get("mgpa_ora_run_meta") or {}
        if run_meta:
            with st.expander("📋 Analysis metadata (for reproducibility)"):
                meta_df = pd.DataFrame([{"Field": k, "Value": str(v)} for k, v in run_meta.items()])
                st.dataframe(meta_df, width="stretch", hide_index=True)
                st.download_button("Download analysis metadata (.csv)", utils.to_download_bytes_csv(meta_df),
                                   file_name=f"{ora_tag}_Analysis_Metadata.csv", mime="text/csv", key="mgpa_dl_ora_meta")

        # ---------------- Step 6: pathway dot plot (ProteoAI ora_dotplot) ----------------
        if len(result) > 0:
            st.subheader("Pathway Dot Plot")
            st.caption("X: Fold Enrichment (dashed line = 1.0). Size: Genes in Overlap. Color: -log10(FDR).")
            d1, d2 = st.columns(2)
            mgpa_top_choice = d1.selectbox("Sets shown", ["Top 10", "Top 20", "Top 25", "Top 30", "All"],
                                           index=2, key="mgpa_top")
            dotplot_top_n = None if mgpa_top_choice == "All" else int(mgpa_top_choice.split()[1])
            mgpa_sort_by = d2.selectbox("Sort by", list(gse.ORA_SORT_OPTIONS.keys()), key="mgpa_sort")
            if st.button("Generate Pathway Dot Plot", type="primary", icon="📌", key="mgpa_dotplot_btn"):
                state["mgpa_dotplot_on"] = True
            if state.get("mgpa_dotplot_on"):
                try:
                    shown = min(len(result), dotplot_top_n or len(result))
                    fig = gse.ora_dotplot(result, top_n=dotplot_top_n, sort_by=mgpa_sort_by)
                    render_single_figure(fig, caption=f"{shown} of {len(result)} gene sets, sorted by "
                                                     f"{mgpa_sort_by}.",
                                         download_name="ORA_Pathway_Dotplot", key_prefix="mgpa_dotplot",
                                         width_ratio=(1, 6, 1))
                except gse.EnrichmentError as e:
                    st.error(str(e))

    # ---------------- Step 7: downloads ----------------
    st.subheader("Downloads")
    tag = f"{_safe(src['comparison'])}_MGPA"
    params = {
        "Statistics source": src_label, "Selection": mode, "Threshold": "" if thr is None else thr,
        "Direction": direction, "P-value column": p_col, "FDR column (input stats)": fdr_col,
        "Method": method, "Organism": organism, "Gene-set database": collection,
        "Database name": resolved_library, "Database version": resolved_version,
        "Background / universe": universe_choice, "Min gene set size": int(min_size_ora),
        "Max gene set size": int(max_size_ora) or "none",
        "Selected metabolites": len(sel_map), "Mapped (HMDB + genes)": n_mapped,
        "HMDB ID but no HMDB protein associations": n_hmdb_only, "Unmapped (no HMDB)": n_unmapped,
        "Metabolite-gene associations": len(mg_df), "Unique genes (input list)": len(ug_df),
        "Distinct symbols tested (after collection symbol harmonization)": len(comparison_genes),
        "ORA engine": "gene_set_enrichment.run_ora_hypergeometric (ProteoAI Pro enrichment.run_ora_hypergeometric)",
        "Metabolite-gene source": f"{hsum['source']}: {hsum['associations']} associations, "
                                  f"{hsum['metabolites']} metabolites, {hsum['genes']} gene labels",
    }
    if result is not None:
        params["ORA results from"] = ("the current settings" if state.get("mgpa_ora_signature") == signature
                                      else "an earlier run (settings changed since)")
        params.update({f"ORA: {k}": v for k, v in (state.get("mgpa_ora_run_meta") or {}).items()})
    params_df = pd.DataFrame({"Value": [str(v) for v in params.values()]},
                             index=pd.Index(list(params.keys()), name="Parameter"))
    sheets = {"Parameters": params_df,
              "Input_Selection": selected,
              "HMDB_Mapping": sel_map.set_index("Metabolite"),
              "Unmapped": not_used.set_index("Metabolite"),
              "Metabolite_Gene_Mapping": mg_df.set_index("Metabolite"),
              "Unique_Gene_List": ug_df.set_index("Gene")}
    if result is not None:
        sheets["ORA_Results"] = result.set_index("Gene Set") if len(result) else result
        uni = state.get("mgpa_ora_universe_genes")
        if uni is not None:
            sheets["Universe_Genes"] = pd.DataFrame({"Gene": sorted(set(uni))}).set_index("Gene")
    b = st.columns([1, 1, 1, 1, 2])
    b[0].download_button("HMDB mapping (CSV)", _csv(sel_map, "Metabolite"), file_name=f"{tag}_HMDB_Mapping.csv",
                         mime="text/csv", key="pa_b_dl_map")
    b[1].download_button("Metabolite–Gene mapping (CSV)", _csv(mg_df, "Metabolite"),
                         file_name=f"{tag}_Metabolite_Gene_Mapping.csv", mime="text/csv", key="pa_b_dl_mg")
    b[2].download_button("Unique gene list (CSV)", _csv(ug_df, "Gene"),
                         file_name=f"{tag}_Unique_Gene_List.csv", mime="text/csv", key="pa_b_dl_ug")
    b[3].download_button("Complete results (Excel)", utils.to_download_bytes_xlsx(sheets),
                         file_name=f"{tag}_Complete_Results.xlsx",
                         mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", key="pa_b_dl_xlsx")
    if fig is not None:
        plt.close(fig)
