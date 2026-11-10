"""
gene_set_enrichment.py - Gene-set Over-Representation Analysis (ORA) engine for
the "Metabolite–Gene Pathway Analysis (MGPA)" tab (metabolite -> gene -> gene-set ORA).

The MGPA tab's Section 4 onward uses the ProteoAI Pro mirror at the END of this module
(ORGANISM_CATALOG, fetch_gene_set_library, run_ora_hypergeometric, ora_dotplot -- ProteoAI's
enrichment.py ORA, verbatim). The extended engine below (run_gene_set_ora: hypergeometric /
Fisher / chi-square, fixed-number universe, threshold filter, overlap matrix) is kept as a
library API but is no longer wired into the UI.

Provenance
----------
This is a port of the ORA engine of the sister Proteomics app (ProteoAI Pro,
metabolomics_modules/enrichment.py -- the same code lineage later ported into the RNA-seq
app's metabolomics_modules/enrichment.py), NOT a new implementation:

  * parse_gmt_text()          -- ProteoAI `parse_gmt_text`, verbatim.
  * run_gene_set_ora()        -- ProteoAI `run_ora_hypergeometric`: identical
                                 N/K/n/k definitions, universe handling (user
                                 background list, or the union of all library
                                 genes when none is given), min/max set-size
                                 filter, Fold Enrichment = (k/n)/(K/N),
                                 Benjamini-Hochberg FDR (statsmodels
                                 multipletests 'fdr_bh') across every gene set
                                 actually tested, results sorted by P-value, and
                                 the same reproducibility `run_meta` dict.
  * neglog10_fdr()            -- ProteoAI `_neglog10_fdr`, verbatim.

It is extended -- exactly as specified by the user's reference gene-set
enrichment script (hypergeom / fisher / chi2, fixed-number or gene-list
universe, P-value threshold, Gene x GeneSet overlap matrix) -- with:

  * hypergeometric_test / fisher_test / chi_square_test -- the reference
    script's three tests, with its exact 2x2 contingency framing and variable
    names:
                         In gene set   Not in gene set
        Comparison            a               b
        Not comparison        c               d
        a = overlap, b = comparison_size - overlap,
        c = gene_set_size - overlap, d = universe_size - gene_set_size - b
  * a fixed-number universe (the reference script's `universe_genes_number`,
    GSEA default 45956): N = that number, gene sets and the comparison list
    are used unrestricted (no gene identities are available to intersect with);
  * a gene-list universe (the reference script's `universe_gene_list_file` /
    ProteoAI's `universe_genes`): N = |universe|, and BOTH the gene sets and
    the comparison list are intersected with it before testing;
  * the reference script's result columns (Gene_set, P-value,
    #genes_in_universe, #genes_in_gene_set, #genes_in_comparison,
    #genes_in_overlap, Overlapped_genes -- comma-joined) plus FDR, Fold
    Enrichment and Gene Ratio (app-family conventions);
  * build_overlap_matrix() -- the reference script's binary Gene x GeneSet
    matrix for gene sets with >= 1 overlapping gene passing the threshold.

Everything runs offline. The bundled MSigDB Hallmark collection
(h.all.v7.0.symbols.gmt, the official Broad Institute MSigDB v7.0 Hallmark
GMT, as redistributed in the GSEApy project's test data; MSigDB is licensed
CC BY 4.0) is the same collection ProteoAI serves as "MSigDB (Hallmark)"
(via Enrichr's MSigDB_Hallmark_2020 library).
"""

import os
import re

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


class EnrichmentError(Exception):
    pass


# ---------------------------------------------------------------------------
# Methods (reference-script identifiers -> display labels)
# ---------------------------------------------------------------------------
METHOD_HYPERGEOM = "hypergeom"
METHOD_FISHER = "fisher"
METHOD_CHI2 = "chi2"
METHOD_LABELS = {
    METHOD_HYPERGEOM: "Hypergeometric test (one-sided)",
    METHOD_FISHER: "Fisher's exact test (one-sided)",
    METHOD_CHI2: "Chi-square test (Pearson, no continuity correction)",
}
METHODS = list(METHOD_LABELS)

# ---------------------------------------------------------------------------
# Universe / background definitions
# ---------------------------------------------------------------------------
GSEA_DEFAULT_UNIVERSE_SIZE = 45956  # reference script's "GSEA default" universe_genes_number

UNI_FIXED = "Fixed number of genes"
UNI_COLLECTION = "All genes in the selected gene-set collection"
UNI_DETECTED = "Genes of all detected/measured metabolites"
UNI_CUSTOM = "Upload custom universe gene list"
UNIVERSE_OPTIONS = [UNI_FIXED, UNI_COLLECTION, UNI_DETECTED, UNI_CUSTOM]
UNIVERSE_HELP = {
    UNI_FIXED: "Background N = a fixed number of genes (default 45956, the GSEA/MSigDB default used by the "
               "reference script). Gene sets and the input gene list are used as-is (no gene identities to "
               "intersect with). Genome-scale backgrounds give the smallest p-values.",
    UNI_COLLECTION: "Background = every gene that appears in at least one gene set of the selected collection "
                    "(ProteoAI's default when no background list is supplied). Input genes that are in no set "
                    "are discarded.",
    UNI_DETECTED: "Background = every gene associated with the HMDB-mapped metabolites measured and tested in "
                  "this dataset -- the gene-level analogue of Option A's 'All detected/measured metabolites'. "
                  "Most conservative: only genes your assay could implicate can be hits.",
    UNI_CUSTOM: "Background = your own gene list (first column; one HGNC symbol per row). Gene sets and the "
                "input list are intersected with it before testing (reference script's universe file mode).",
}

# ---------------------------------------------------------------------------
# Gene-set collections
# ---------------------------------------------------------------------------
COLL_HALLMARK = "MSigDB (Hallmark)"
COLL_CURATED = "Curated metabolic pathway gene sets (embedded)"
COLL_CUSTOM = "Upload custom .gmt"
COLLECTION_OPTIONS = [COLL_HALLMARK, COLL_CURATED, COLL_CUSTOM]

_BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HALLMARK_GMT_PATH = os.path.join(_BASE_DIR, "h.all.v7.0.symbols.gmt")
HALLMARK_VERSION = "MSigDB v7.0 Hallmark collection (h.all.v7.0.symbols.gmt; 50 gene sets, HGNC symbols)"

ORA_RESULT_COLUMNS = ["Gene_set", "P-value", "FDR", "#genes_in_universe", "#genes_in_gene_set",
                      "#genes_in_comparison", "#genes_in_overlap", "Fold_Enrichment", "Gene_Ratio",
                      "Overlapped_genes"]


# ---------------------------------------------------------------------------
# Readers
# ---------------------------------------------------------------------------
def parse_gmt_text(text: str) -> dict:
    """
    Parse GMT-format text (the standard MSigDB/Enrichr convention:
    term<TAB>description<TAB>gene1<TAB>gene2<TAB>...) into {term: [gene_symbols]}.
    (Ported verbatim from the Proteomics app's proteomics_modules/enrichment.py.)
    """
    gene_sets = {}
    for line in text.strip().splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue
        term = parts[0].strip()
        genes = [g.strip().upper() for g in parts[2:] if g.strip()]
        if genes:
            gene_sets[term] = genes
    return gene_sets


_HALLMARK_CACHE = {}


def load_hallmark_gene_sets() -> dict:
    """The bundled MSigDB Hallmark collection as {term: [genes]} (parsed once, then cached)."""
    if "sets" not in _HALLMARK_CACHE:
        if not os.path.exists(HALLMARK_GMT_PATH):
            raise EnrichmentError(f"Bundled Hallmark GMT not found at {HALLMARK_GMT_PATH}.")
        with open(HALLMARK_GMT_PATH, "r", encoding="utf-8") as f:
            _HALLMARK_CACHE["sets"] = parse_gmt_text(f.read())
    return {k: list(v) for k, v in _HALLMARK_CACHE["sets"].items()}


_HEADER_TOKENS = {"gene", "genes", "symbol", "gene_symbol", "gene symbol", "genesymbol", "hgnc_symbol",
                  "hgnc symbol", "gene_name", "gene name", "id", "name"}


def read_gene_list_text(text: str, filename: str = "") -> list:
    """
    Gene symbols from the FIRST column of a text/CSV/TSV file (the reference
    script's read_gene_list convention), stripped, upper-cased and de-duplicated
    in file order. Parsed line by line (first field before a tab, comma,
    semicolon or whitespace) rather than with delimiter sniffing, which can
    mis-detect a delimiter inside one-column symbol lists (e.g. 'HLA-A').
    A first-row header such as 'Gene' / 'Symbol' is skipped.
    """
    genes, first = [], True
    for line in text.lstrip("\ufeff").splitlines():
        line = line.strip()
        if not line:
            continue
        g = re.split(r"[\t,;]|\s+", line, maxsplit=1)[0].strip().strip('"').strip("'").strip()
        if first:
            first = False
            if g.lower() in _HEADER_TOKENS:
                continue
        if g:
            genes.append(g.upper())
    return list(dict.fromkeys(genes))


# ---------------------------------------------------------------------------
# The three tests (reference script's functions and 2x2 framing)
# ---------------------------------------------------------------------------
def _contingency(universe_size, gene_set_size, comparison_size, overlap):
    a = overlap
    b = comparison_size - overlap
    c = gene_set_size - overlap
    d = universe_size - gene_set_size - b
    return a, b, c, d


def hypergeometric_test(universe_size, gene_set_size, comparison_size, overlap) -> float:
    """One-sided hypergeometric P(X >= overlap) = hypergeom.sf(overlap-1, N, K, n); 1.0 when overlap <= 0."""
    if overlap <= 0:
        return 1.0
    return float(stats.hypergeom.sf(overlap - 1, universe_size, gene_set_size, comparison_size))


def fisher_test(universe_size, gene_set_size, comparison_size, overlap) -> float:
    """One-sided (alternative='greater') Fisher's exact test on [[a, b], [c, d]]; NaN for an invalid table."""
    a, b, c, d = _contingency(universe_size, gene_set_size, comparison_size, overlap)
    if min(a, b, c, d) < 0:
        return np.nan
    return float(stats.fisher_exact(np.array([[a, b], [c, d]]), alternative="greater")[1])


def chi_square_test(universe_size, gene_set_size, comparison_size, overlap) -> float:
    """
    Pearson chi-square test on [[a, b], [c, d]] with correction=False -- the
    reference script's exact call (scipy.stats.chi2_contingency), whose p-value
    is two-sided (it is returned unchanged, i.e. conservatively, for enriched
    sets). Because this is an OVER-representation analysis, a gene set whose
    observed overlap does not exceed its expectation (a <= (a+b)(a+c)/N, i.e.
    neutral or depleted) is reported as P = 1.0 rather than letting depletion
    masquerade as enrichment. NaN for an invalid/degenerate table.
    """
    a, b, c, d = _contingency(universe_size, gene_set_size, comparison_size, overlap)
    if min(a, b, c, d) < 0:
        return np.nan
    n_tot = a + b + c + d
    if n_tot <= 0:
        return np.nan
    expected_a = (a + b) * (a + c) / n_tot
    try:
        p = float(stats.chi2_contingency(np.array([[a, b], [c, d]]), correction=False)[1])
    except ValueError:
        return np.nan
    if np.isnan(p):
        return np.nan
    return p if a > expected_a else 1.0


_TESTS = {METHOD_HYPERGEOM: hypergeometric_test, METHOD_FISHER: fisher_test, METHOD_CHI2: chi_square_test}


def _bh_fdr(pvals) -> np.ndarray:
    """BH-FDR via statsmodels multipletests('fdr_bh') as in ProteoAI; NaN p-values stay NaN (excluded)."""
    p = np.asarray(pvals, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    if ok.any():
        out[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return out


# ---------------------------------------------------------------------------
# ORA (ProteoAI run_ora_hypergeometric, extended with method + fixed universe)
# ---------------------------------------------------------------------------
def run_gene_set_ora(comparison_genes, gene_sets: dict, method: str = METHOD_HYPERGEOM,
                     universe_genes=None, universe_size: int = None, min_set_size: int = 1,
                     max_set_size: int = None):
    """
    Over-Representation Analysis of `comparison_genes` in every gene set.

    comparison_genes: the FULL input gene list (here: the Unique Gene List of the
        selected metabolites) -- used exactly as given (duplicates collapsed).
    gene_sets: {term: [gene_symbols]}.
    method: 'hypergeom' | 'fisher' | 'chi2'.
    Universe (exactly one of these, checked in this order):
      universe_genes (iterable) -> N = |universe|; K = |gene set ∩ universe|;
          n = |comparison ∩ universe| (reference script's universe-file mode /
          ProteoAI's user-provided background).
      universe_size (int)       -> N = universe_size; K = |gene set|;
          n = |comparison| (reference script's fixed-number mode). N must be at
          least the number of distinct genes in the input list and the collection.
      neither                   -> universe = union of all genes in `gene_sets`
          (ProteoAI's default).
    min_set_size / max_set_size: skip gene sets whose K falls outside this range.

    For each gene set with K > 0 inside the size limits (sets with no
    overlapping gene are tested and reported too, P = 1, as in both ProteoAI and
    the reference script):
        k = |comparison ∩ gene set|, P-value by `method`,
        Fold_Enrichment = (k/n)/(K/N), Gene_Ratio = k/n.
    FDR = Benjamini-Hochberg across every gene set tested.

    Returns (result_df sorted by P-value ascending, run_meta dict of reproducibility
    metadata in ProteoAI's format, info dict of N/n and genes not in the collection/universe).
    """
    if method not in _TESTS:
        raise EnrichmentError("Invalid method. Use 'hypergeom', 'fisher', or 'chi2'.")
    if not gene_sets:
        raise EnrichmentError("No gene sets to test (the selected collection is empty).")
    comparison_raw = [str(g).strip() for g in comparison_genes if str(g).strip()]
    if len(comparison_raw) == 0:
        raise EnrichmentError("Comparison gene list is empty.")
    comparison = set(g.upper() for g in comparison_raw)
    n_duplicates = len(comparison_raw) - len(comparison)

    sets_upper = {term: set(g.upper() for g in genes) for term, genes in gene_sets.items()}
    all_library_genes = set().union(*sets_upper.values())

    restrict = True
    if universe_genes is not None:
        universe = set(str(g).strip().upper() for g in universe_genes if str(g).strip())
        if not universe:
            raise EnrichmentError("The universe/background gene list is empty.")
        universe_source = f"user-provided background gene list ({len(universe)} genes)"
        N = len(universe)
    elif universe_size is not None:
        N = int(universe_size)
        required = len(comparison | all_library_genes)
        if N < required:
            raise EnrichmentError(
                f"Fixed universe size {N} is smaller than the number of distinct genes in the input list "
                f"and the gene-set collection combined ({required}); increase it.")
        universe = None
        restrict = False
        universe_source = f"fixed number of genes (N = {N}; gene sets and input list not restricted)"
    else:
        universe = set(all_library_genes)
        universe_source = f"union of all genes in the selected gene-set collection ({len(universe)} genes)"
        N = len(universe)

    comparison_in_universe = comparison & universe if restrict else set(comparison)
    n = len(comparison_in_universe)
    n_comparison_dropped = len(comparison) - n
    test = _TESTS[method]

    rows = []
    n_sets_skipped_size = n_sets_empty = 0
    if n > 0:
        for term, genes in sets_upper.items():
            gene_set_genes = genes & universe if restrict else genes
            K = len(gene_set_genes)
            if K == 0:
                n_sets_empty += 1
                continue
            if K < min_set_size or (max_set_size is not None and K > max_set_size):
                n_sets_skipped_size += 1
                continue
            overlap = comparison_in_universe & gene_set_genes
            k = len(overlap)
            pval = test(universe_size=N, gene_set_size=K, comparison_size=n, overlap=k)
            rows.append({
                "Gene_set": term, "P-value": pval, "#genes_in_universe": N, "#genes_in_gene_set": K,
                "#genes_in_comparison": n, "#genes_in_overlap": k,
                "Fold_Enrichment": (k / n) / (K / N) if (n > 0 and K > 0 and N > 0) else 0.0,
                "Gene_Ratio": k / n if n > 0 else 0.0,
                "Overlapped_genes": ",".join(sorted(overlap)),
            })

    result_df = pd.DataFrame(rows, columns=[c for c in ORA_RESULT_COLUMNS if c != "FDR"])
    if len(result_df):
        result_df["FDR"] = _bh_fdr(result_df["P-value"].values)
        result_df = (result_df[ORA_RESULT_COLUMNS]
                     .sort_values("P-value", ascending=True, kind="mergesort", na_position="last")
                     .reset_index(drop=True))
    else:
        result_df = pd.DataFrame(columns=ORA_RESULT_COLUMNS)

    run_meta = {
        "Analysis method": "Over-Representation Analysis (ORA), gene-set level",
        "Statistical test": METHOD_LABELS[method],
        "Multiple-testing correction": "Benjamini-Hochberg FDR (across all gene sets tested)",
        "Background/universe definition": universe_source,
        "# genes in universe (N)": N,
        "# genes in input comparison list": len(comparison_raw),
        "# duplicate genes in input (collapsed)": n_duplicates,
        "# genes in comparison list used for testing (n)": n,
        "# genes in comparison list NOT found in universe (discarded)": n_comparison_dropped,
        "# genes in comparison list present in the collection": len(comparison & all_library_genes),
        "# genes in comparison list absent from the collection": len(comparison - all_library_genes),
        "# unique genes in collection": len(all_library_genes),
        "# gene sets in library": len(gene_sets),
        "# gene sets tested (after size filter)": len(result_df),
        "# gene sets skipped (outside min/max size)": n_sets_skipped_size,
        "# gene sets with no gene in the universe": n_sets_empty,
    }
    info = {"N": N, "n": n, "method": method, "universe_source": universe_source,
            "genes_not_in_collection": sorted(comparison - all_library_genes),
            "genes_not_in_universe": sorted(comparison - comparison_in_universe),
            "n_sets_total": len(gene_sets), "n_sets_tested": len(result_df)}
    return result_df, run_meta, info


# ---------------------------------------------------------------------------
# Threshold filter + overlap matrix (reference script sections 10-11)
# ---------------------------------------------------------------------------
def filter_results(result_df: pd.DataFrame, metric: str = "P-value", threshold: float = 1.0,
                   apply: bool = True) -> pd.DataFrame:
    """Rows with `metric` <= threshold (inclusive, NaN never passes) when `apply`, else all rows."""
    if not apply or result_df is None or len(result_df) == 0:
        return result_df
    return result_df[(result_df[metric] <= threshold).fillna(False)]


def build_overlap_matrix(result_df: pd.DataFrame, metric: str = "P-value", threshold: float = 1.0,
                         apply: bool = True) -> pd.DataFrame:
    """
    Binary Gene x GeneSet matrix (1 = gene is an overlapping input gene of that
    set) over gene sets with >= 1 overlapping gene, a non-NaN p-value and
    `metric` <= threshold -- the reference script's `.overlap` output. Rows
    (genes) and columns (gene sets) are sorted; index name 'Gene'. Empty frame
    if no gene set qualifies.
    """
    if result_df is None or len(result_df) == 0:
        return pd.DataFrame(index=pd.Index([], name="Gene"))
    ok = (result_df["#genes_in_overlap"] > 0) & result_df["P-value"].notna()
    if apply:
        ok &= (result_df[metric] <= threshold).fillna(False)
    sub = result_df[ok]
    overlap = {r["Gene_set"]: set(filter(None, str(r["Overlapped_genes"]).split(","))) for _, r in sub.iterrows()}
    if not overlap:
        return pd.DataFrame(index=pd.Index([], name="Gene"))
    genes = sorted(set().union(*overlap.values()))
    mat = pd.DataFrame(0, index=genes, columns=sorted(overlap), dtype=int)
    for term, gs in overlap.items():
        mat.loc[sorted(gs), term] = 1
    mat.index.name = "Gene"
    return mat


def neglog10_fdr(fdr_series: pd.Series) -> pd.Series:
    """-log10(FDR) with exact zeros floored to the smallest positive float (ProteoAI `_neglog10_fdr`)."""
    return -np.log10(fdr_series.astype(float).clip(lower=np.finfo(float).tiny))


def dotplot_view(result_df: pd.DataFrame) -> pd.DataFrame:
    """
    Column adapter so gene-set ORA results feed the tab's shared
    pathway_analysis.pathway_dotplot() unchanged (Y = gene set, X = Gene ratio
    or Fold enrichment, size = overlap count, color = -log10(FDR)). Only sets
    with >= 1 overlapping gene and a valid p-value are plotted.
    """
    df = result_df[(result_df["#genes_in_overlap"] > 0) & result_df["P-value"].notna()]
    return pd.DataFrame({
        "Pathway": df["Gene_set"].values, "P-value": df["P-value"].values, "FDR": df["FDR"].values,
        "Overlap_Count": df["#genes_in_overlap"].values, "Enrichment_Ratio": df["Fold_Enrichment"].values,
        "Gene_Ratio": df["Gene_Ratio"].values,
    })


# ===========================================================================
# ProteoAI Pro "Over-Representation Analysis (ORA)" -- mirrored for the MGPA tab
# ===========================================================================
# Ported from the Proteomics app's proteomics_modules/enrichment.py (ORGANISM_CATALOG, ORA_GSEA_DATABASE_COLLECTIONS,
# get_database_info, get_enrichr_gene_set_library, fetch_gene_set_library, ORA_RESULT_COLUMNS,
# run_ora_hypergeometric, _enrichment_dotplot, ora_dotplot). The statistics, result columns, run
# metadata and dot plot are ProteoAI's verbatim. MetaboAI-specific differences, all deliberate:
#   * only the Human catalog entry -- HMDB metabolite-protein associations are human (HGNC symbols);
#   * "MSigDB (Hallmark)" is served from the bundled offline GMT (h.all.v7.0.symbols.gmt) instead of
#     Enrichr's MSigDB_Hallmark_2020, so the tab works without internet; every other collection is
#     fetched live from Enrichr exactly as ProteoAI does (internet required), cached per session.

import math  # noqa: E402
import textwrap  # noqa: E402
import matplotlib  # noqa: E402
import matplotlib.figure  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.gridspec as gridspec  # noqa: E402
from scipy.stats import hypergeom as _hypergeom  # noqa: E402

matplotlib.use("Agg")

ENRICHR_BASE = "https://maayanlab.cloud/Enrichr"
YEAST_ENRICHR_BASE = "https://maayanlab.cloud/YeastEnrichr"
FISH_ENRICHR_BASE = "https://maayanlab.cloud/FishEnrichr"
ORGANISM_HUMAN = "Human (Homo sapiens)"
ORGANISM_CATALOG = {
    ORGANISM_HUMAN: {
        "string_taxid": 9606,
        "enrichr_base": ENRICHR_BASE,
        "gene_id_type": "Gene Symbol (HGNC)",
        "databases": {
            "GO Biological Process": {"library": "GO_Biological_Process_2026", "version": "2026", "available": True},
            "GO Molecular Function": {"library": "GO_Molecular_Function_2026", "version": "2026", "available": True},
            "GO Cellular Component": {"library": "GO_Cellular_Component_2026", "version": "2026", "available": True},
            "KEGG": {"library": "KEGG_2026", "version": "2026", "available": True},
            "Reactome": {"library": "Reactome_Pathways_2024", "version": "2024 (most recent Enrichr release; no 2026 Reactome library exists as of this writing)", "available": True},
            "WikiPathways": {"library": "WikiPathways_2024_Human", "version": "2024 (most recent Enrichr release; no 2026 release exists)", "available": True},
            "MSigDB (Hallmark)": {"library": "h.all.v7.0.symbols.gmt", "version": HALLMARK_VERSION + " -- bundled offline copy", "available": True, "bundled": True},
        },
    },
    "Mouse (Mus musculus)": {
        "string_taxid": 10090,
        "enrichr_base": ENRICHR_BASE,
        "gene_id_type": "Gene Symbol (MGI)",
        "databases": {
            "GO Biological Process": {"library": None, "version": None, "available": False,
                                       "note": "Enrichr's GO libraries are human-gene-symbol only; no mouse-specific GO release exists on Enrichr."},
            "GO Molecular Function": {"library": None, "version": None, "available": False,
                                       "note": "Enrichr's GO libraries are human-gene-symbol only; no mouse-specific GO release exists on Enrichr."},
            "GO Cellular Component": {"library": None, "version": None, "available": False,
                                       "note": "Enrichr's GO libraries are human-gene-symbol only; no mouse-specific GO release exists on Enrichr."},
            "KEGG": {"library": "KEGG_2019_Mouse", "version": "2019 (most recent mouse-specific KEGG library on Enrichr; no 2026 release exists)", "available": True},
            "Reactome": {"library": None, "version": None, "available": False,
                         "note": "No mouse-specific Reactome library exists on Enrichr."},
            "WikiPathways": {"library": "WikiPathways_2024_Mouse", "version": "2024 (most recent Enrichr release; no 2026 release exists)", "available": True},
            "MSigDB (Hallmark)": {"library": None, "version": None, "available": False,
                                  "note": "Enrichr's MSigDB Hallmark library is human-only."},
        },
    },
    "Rat (Rattus norvegicus)": {
        "string_taxid": 10116,
        "enrichr_base": None,
        "gene_id_type": "Gene Symbol (RGD)",
        "databases": {c: {"library": None, "version": None, "available": False,
                          "note": "Enrichr has no rat-specific gene-set libraries or organism portal."}
                     for c in ["GO Biological Process", "GO Molecular Function", "GO Cellular Component",
                               "KEGG", "Reactome", "WikiPathways", "MSigDB (Hallmark)"]},
    },
    "Yeast (S. cerevisiae)": {
        "string_taxid": 4932,
        "enrichr_base": YEAST_ENRICHR_BASE,
        "gene_id_type": "Gene Symbol (SGD)",
        "databases": {
            "GO Biological Process": {"library_prefix": "GO_Biological_Process", "version": "current YeastEnrichr release (organism-specific Enrichr instance; not year-stamped -- exact name resolved at query time)", "available": True},
            "GO Molecular Function": {"library_prefix": "GO_Molecular_Function", "version": "current YeastEnrichr release (not year-stamped)", "available": True},
            "GO Cellular Component": {"library_prefix": "GO_Cellular_Component", "version": "current YeastEnrichr release (not year-stamped)", "available": True},
            "KEGG": {"library": None, "version": None, "available": False,
                     "note": "No confirmed KEGG library on YeastEnrichr."},
            "Reactome": {"library": None, "version": None, "available": False,
                        "note": "No confirmed Reactome library on YeastEnrichr."},
            "WikiPathways": {"library_prefix": "WikiPathways", "version": "current YeastEnrichr release (not year-stamped)", "available": True},
            "MSigDB (Hallmark)": {"library": None, "version": None, "available": False, "note": "No MSigDB library on YeastEnrichr."},
        },
    },
    "Zebrafish (Danio rerio)": {
        "string_taxid": 7955,
        "enrichr_base": FISH_ENRICHR_BASE,
        "gene_id_type": "Gene Symbol (ZFIN)",
        "databases": {
            "GO Biological Process": {"library_prefix": "GO_Biological_Process", "version": "current FishEnrichr release (organism-specific Enrichr instance; not year-stamped -- exact name resolved at query time)", "available": True},
            "GO Molecular Function": {"library_prefix": "GO_Molecular_Function", "version": "current FishEnrichr release (not year-stamped)", "available": True},
            "GO Cellular Component": {"library_prefix": "GO_Cellular_Component", "version": "current FishEnrichr release (not year-stamped)", "available": True},
            "KEGG": {"library": None, "version": None, "available": False,
                     "note": "No confirmed KEGG library on FishEnrichr."},
            "Reactome": {"library": None, "version": None, "available": False, "note": "No confirmed Reactome library on FishEnrichr."},
            "WikiPathways": {"library_prefix": "WikiPathways", "version": "current FishEnrichr release (not year-stamped)", "available": True},
            "MSigDB (Hallmark)": {"library": None, "version": None, "available": False, "note": "No MSigDB library on FishEnrichr."},
        },
    },
}

# Major plant species: no Enrichr-family portal exists for any plant genome, so ORA is unavailable
# for all of them in this app.
_PLANT_SPECIES = {
    "Arabidopsis thaliana": 3702, "Rice (Oryza sativa)": 39947,
    "Maize (Zea mays)": 4577, "Soybean (Glycine max)": 3847,
    "Tomato (Solanum lycopersicum)": 4081,
}
for _name, _taxid in _PLANT_SPECIES.items():
    ORGANISM_CATALOG[f"{_name} (plant)"] = {
        "string_taxid": _taxid, "enrichr_base": None, "gene_id_type": "Gene Symbol / Locus ID",
        "databases": {c: {"library": None, "version": None, "available": False,
                          "note": "No Enrichr-family portal exists for plant species."}
                     for c in ["GO Biological Process", "GO Molecular Function", "GO Cellular Component",
                               "KEGG", "Reactome", "WikiPathways", "MSigDB (Hallmark)"]},
    }

ORA_GSEA_DATABASE_COLLECTIONS = ["GO Biological Process", "GO Molecular Function", "GO Cellular Component",
                                  "KEGG", "Reactome", "WikiPathways", "MSigDB (Hallmark)"]
DEFAULT_ORA_COLLECTION = "MSigDB (Hallmark)"   # the one collection available offline

UNIVERSE_DETECTED_LABEL = "All genes of detected/measured metabolites in this study (recommended)"
UNIVERSE_LIBRARY_LABEL = "Entire gene-set library"
ORA_UNIVERSE_OPTIONS = [UNIVERSE_DETECTED_LABEL, UNIVERSE_LIBRARY_LABEL]

_LIBRARY_CACHE = {}


def get_database_info(organism: str, collection: str) -> dict:
    """Look up one organism/collection cell of ORGANISM_CATALOG (ProteoAI)."""
    org = ORGANISM_CATALOG.get(organism)
    if org is None:
        raise EnrichmentError(f"Unsupported organism: '{organism}'.")
    info = org["databases"].get(collection)
    if info is None:
        raise EnrichmentError(f"Unsupported database collection: '{collection}'.")
    return info


def get_enrichr_gene_set_library(library: str, enrichr_base: str = ENRICHR_BASE):
    """Fetch a full Enrichr gene-set library in GMT format and parse it into {term: [gene_symbols]} (ProteoAI)."""
    try:
        import requests
    except ImportError as e:  # pragma: no cover
        raise EnrichmentError(f"The 'requests' package is required to fetch Enrichr libraries: {e}")
    try:
        resp = requests.get(f"{enrichr_base}/geneSetLibrary",
                            params={"mode": "text", "libraryName": library}, timeout=60)
        resp.raise_for_status()
    except requests.exceptions.RequestException as e:
        raise EnrichmentError(f"Could not reach {enrichr_base} — check your internet connection. Details: {e}")
    if not resp.text.strip():
        raise EnrichmentError(f"Library '{library}' returned no data from {enrichr_base}.")
    return parse_gmt_text(resp.text)


def fetch_gene_set_library(organism: str, collection: str):
    """
    Resolve and fetch the gene-set library of one (organism, collection) cell (ProteoAI), cached per process.
    Returns (gene_sets {term: [genes]}, resolved_library_name, documented_version).
    """
    info = get_database_info(organism, collection)
    if not info.get("available", False):
        raise EnrichmentError(f"'{collection}' is not available for {organism}: "
                              f"{info.get('note', 'no database available for this combination.')}")
    # Most organisms give an exact, year-stamped "library" name. Yeast/zebrafish libraries live on
    # their own organism-specific Enrichr instance and are not year-stamped there, so only a
    # "library_prefix" is known ahead of time -- it is tried as the literal library name (which is
    # also the un-dated name these portals commonly use).
    library_name = info.get("library") or info.get("library_prefix")
    if library_name is None:
        raise EnrichmentError(f"No gene-set library is configured for '{collection}' ({organism}).")
    key = (organism, collection)
    if key not in _LIBRARY_CACHE:
        if info.get("bundled"):
            gene_sets = load_hallmark_gene_sets()
        else:
            try:
                gene_sets = get_enrichr_gene_set_library(library_name, ORGANISM_CATALOG[organism]["enrichr_base"])
            except EnrichmentError:
                if info.get("library_prefix") and not info.get("library"):
                    raise EnrichmentError(
                        f"Could not fetch '{library_name}' from the {organism} Enrichr portal -- the current "
                        "library name there may differ from the expected prefix.")
                raise
        _LIBRARY_CACHE[key] = gene_sets
    return {k: list(v) for k, v in _LIBRARY_CACHE[key].items()}, library_name, info["version"]


PROTEOAI_ORA_COLUMNS = ["Gene Set", "P-value", "FDR", "# genes in universe", "# Genes in Gene Set",
                        "# genes in comparison", "# Genes in Overlap", "Fold Enrichment", "Overlapped genes"]


def run_ora_hypergeometric(comparison_genes, gene_sets: dict, universe_genes=None,
                           min_set_size: int = 1, max_set_size: int = None):
    """
    Over-Representation Analysis via the hypergeometric test -- ProteoAI Pro's run_ora_hypergeometric,
    verbatim (here the comparison list is the MGPA Unique Gene List).
        N = |universe|, K = |gene set ∩ universe|, n = |comparison ∩ universe|,
        k = |comparison ∩ gene set ∩ universe|, P = hypergeom.sf(k-1, N, K, n),
        Fold Enrichment = (k/n)/(K/N); FDR = Benjamini-Hochberg across every gene set tested.
    universe_genes None -> the union of every gene in `gene_sets`.
    Returns (result_df with PROTEOAI_ORA_COLUMNS, run_meta).
    """
    if not gene_sets:
        raise EnrichmentError("No gene sets to test (the selected library returned nothing).")
    comparison_raw = [str(g).strip() for g in comparison_genes if str(g).strip()]
    if len(comparison_raw) == 0:
        raise EnrichmentError("Comparison gene list is empty.")
    comparison = set(g.upper() for g in comparison_raw)
    n_duplicates = len(comparison_raw) - len(comparison)

    all_library_genes = set()
    for genes in gene_sets.values():
        all_library_genes.update(g.upper() for g in genes)

    if universe_genes is not None:
        universe_raw = [str(g).strip() for g in universe_genes if str(g).strip()]
        universe = set(g.upper() for g in universe_raw)
        universe_source = (f"user-provided background ({len(universe)} genes -- all genes of the metabolites "
                           "detected/measured in this study)")
    else:
        universe = set(all_library_genes)
        universe_source = f"union of all genes in the selected gene-set library ({len(universe)} genes)"

    comparison_in_universe = comparison & universe
    n_comparison_dropped = len(comparison) - len(comparison_in_universe)
    N = len(universe)
    n = len(comparison_in_universe)

    rows = []
    n_sets_skipped_size = 0
    if n > 0:
        for term, genes in gene_sets.items():
            gene_set_genes = set(g.upper() for g in genes) & universe
            K = len(gene_set_genes)
            if K == 0:
                continue
            if K < min_set_size or (max_set_size is not None and K > max_set_size):
                n_sets_skipped_size += 1
                continue
            overlap = comparison_in_universe & gene_set_genes
            k = len(overlap)
            pval = _hypergeom.sf(k - 1, N, K, n) if k > 0 else 1.0
            fold_enrichment = (k / n) / (K / N) if (n > 0 and K > 0 and N > 0) else 0.0
            rows.append({
                "Gene Set": term, "P-value": pval, "# genes in universe": N,
                "# Genes in Gene Set": K, "# genes in comparison": n,
                "# Genes in Overlap": k, "Fold Enrichment": fold_enrichment,
                "Overlapped genes": ";".join(sorted(overlap)),
            })

    result_df = pd.DataFrame(rows, columns=[c for c in PROTEOAI_ORA_COLUMNS if c != "FDR"])
    if len(result_df):
        result_df["FDR"] = multipletests(result_df["P-value"], method="fdr_bh")[1]
        result_df = result_df[PROTEOAI_ORA_COLUMNS].sort_values("P-value").reset_index(drop=True)
    else:
        result_df = pd.DataFrame(columns=PROTEOAI_ORA_COLUMNS)

    run_meta = {
        "Analysis method": "Over-Representation Analysis (ORA)",
        "Statistical test": "Hypergeometric test, one-sided (P(X >= k overlapping genes))",
        "Multiple-testing correction": "Benjamini-Hochberg FDR",
        "Background/universe definition": universe_source,
        "# genes in universe (N)": N,
        "# genes in input comparison list": len(comparison_raw),
        "# duplicate genes in input (collapsed)": n_duplicates,
        "# genes in comparison list mapped to universe (n)": n,
        "# genes in comparison list NOT found in universe (discarded)": n_comparison_dropped,
        "# gene sets in library": len(gene_sets),
        "# gene sets tested (after size filter)": len(result_df),
        "# gene sets skipped (outside min/max size)": n_sets_skipped_size,
    }
    return result_df, run_meta


def _enrichment_dotplot(result_df: pd.DataFrame, term_col: str, x_col: str, size_col: str,
                        color_col: str, x_label: str, title: str, size_legend_label: str = "Count",
                        color_legend_label: str = None, top_n: int = 15, reference_x: float = None,
                        select_col: str = None, select_ascending: bool = True):
    # top_n=None (or 0) means "show all rows" -- mirrors pathway_analysis.order_results' convention.
    """
    Shared enrichment dot-plot renderer (top_n by select_col, ordered by x_col), sized as a compact,
    publication-quality figure. Mirrors pathway_analysis.pathway_dotplot's layout: a GridSpec splits
    the plot from a dedicated legend column (size key on top, colour scale below), the legend block is
    a fixed height in inches and vertically centred regardless of row count, and legend dots are spaced
    by their own diameter (converted from points to inches) plus a small fixed pad, so they sit with an
    even gap and never overlap or crowd together. Long term names wrap onto up to 3 lines rather than
    being truncated or cut off.
    """
    select_col = select_col or color_col
    if result_df is None or len(result_df) == 0:
        raise EnrichmentError("No enrichment results to plot.")
    required = {term_col, x_col, size_col, color_col, select_col}
    missing = required - set(result_df.columns)
    if missing:
        raise EnrichmentError(f"Enrichment results are missing required column(s): {', '.join(sorted(missing))}")

    work = result_df.dropna(subset=[x_col, size_col, color_col, select_col]).copy()
    if len(work) == 0:
        raise EnrichmentError("No enrichment results with complete (non-missing) data to plot.")

    ranked = work.sort_values(select_col, ascending=select_ascending)
    subset = ranked.head(top_n) if top_n else ranked
    subset = subset.sort_values(x_col, ascending=True)
    n = len(subset)
    if n == 0:
        raise EnrichmentError("No gene sets available to plot.")

    def _wrap(name: str) -> list:
        lines = textwrap.wrap(str(name), 40) or [str(name)]
        if len(lines) > 3:
            lines = lines[:3]
            lines[-1] = lines[-1].rstrip() + " …"
        return lines

    line_lists = [_wrap(t) for t in subset[term_col]]
    labels = ["\n".join(ll) for ll in line_lists]
    bands = np.array([len(ll) for ll in line_lists], dtype=float)
    total_units = float(bands.sum())
    edges = np.concatenate([[0.0], np.cumsum(bands)])
    centers = (edges[:-1] + edges[1:]) / 2.0
    y = centers  # subset is already ascending by x_col, so first row (lowest x) sits at the bottom

    MARGIN_SCALE = 0.60
    title_text = title or "Enrichment Dot Plot"
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

    cmap = matplotlib.colormaps.get_cmap("coolwarm")  # blue = low value, red = high value
    cvals = subset[color_col].astype(float)
    vmin, vmax = float(cvals.min()), float(cvals.max())
    if np.isclose(vmin, vmax):
        vmin, vmax = vmin - 0.5, vmax + 0.5
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)

    counts = subset[size_col].astype(float)
    cmin, cmax = counts.min(), counts.max()

    def _s(v):
        return 90.0 if cmax == cmin else 28.0 + 150.0 * (v - cmin) / (cmax - cmin)

    for yi in y:
        ax.axhline(yi, color="#EEEEEE", linewidth=0.6, zorder=0)
    ax.scatter(subset[x_col], y, s=counts.map(_s), c=cvals, cmap=cmap, norm=norm,
               edgecolor="#333333", linewidth=0.5, zorder=3)
    if reference_x is not None:
        ax.axvline(reference_x, color="#808080", linestyle="--", linewidth=0.8, zorder=1)

    xvals = subset[x_col].astype(float).values
    xmin, xmax = float(np.nanmin(xvals)), float(np.nanmax(xvals))
    if reference_x is not None:
        xmin, xmax = min(xmin, reference_x), max(xmax, reference_x)
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

    # ---- Legend column: size key (top) then colour scale (bottom), vertically centred ----
    axis_height_in = max((top - bottom) * height, 0.1)

    def rows(inches):
        return inches / axis_height_in

    integer_sizes = bool(np.allclose(counts, np.round(counts)))
    reps = sorted(set(int(round(v)) for v in np.linspace(cmin, cmax, 3))) if integer_sizes else \
        sorted(set(float(f"{v:.2g}") for v in np.linspace(cmin, cmax, 3)))

    # Even, non-overlapping dot spacing: the gap between legend dots is set to the LARGEST legend
    # dot's own diameter (converted from points to inches) plus a small fixed pad.
    max_area_pts2 = max(_s(c) for c in reps)
    dot_diam_in = 2.0 * math.sqrt(max_area_pts2 / math.pi) / 72.0
    dot_step_in = dot_diam_in + 0.045
    dot_step = rows(dot_step_in)

    size_title_lines = size_legend_label if len(size_legend_label) <= 15 else \
        "\n".join(textwrap.wrap(size_legend_label, 13))
    size_title_h = 0.12 * (size_title_lines.count("\n") + 1) + 0.06
    dots_h = dot_step_in * len(reps)
    rule_h = 0.04 + 0.13
    # Title line + clearance for the top colour-bar tick label (centred on the bar's top edge).
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

    color_title = color_legend_label or color_col
    ax_leg.text(0.02, y_cursor, color_title, fontsize=7.8, fontweight="bold", va="top", ha="left")
    y_cursor -= rows(color_title_h)
    cbar_h = rows(cbar_in)
    cax = ax_leg.inset_axes([0.06, max(y_cursor - cbar_h, 0.0), 0.22, cbar_h])
    sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cb = fig.colorbar(sm, cax=cax)
    cb.ax.tick_params(labelsize=7.8, length=2.5, pad=2)
    cb.ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=4))
    cb.outline.set_linewidth(0.5)
    return fig


#: "Sort by" choices offered for the ORA pathway dot plot, mirroring pathway_analysis.SORT_OPTIONS
#: (label -> (column to rank/select by, ascending)).
ORA_SORT_OPTIONS = {"P-value": ("P-value", True), "FDR": ("FDR", True),
                    "Fold Enrichment": ("Fold Enrichment", False)}


def ora_dotplot(result_df: pd.DataFrame, top_n: int = 15, sort_by: str = "FDR"):
    """
    ProteoAI Pro's ORA dot plot, verbatim: x = Fold Enrichment (dashed line at 1.0), size = Genes in
    Overlap, colour = -log10(FDR); the top_n gene sets ranked by `sort_by` (see ORA_SORT_OPTIONS;
    default FDR, ascending), ordered by Fold Enrichment on the plot itself. top_n=None shows all.
    """
    if result_df is None or len(result_df) == 0:
        raise EnrichmentError("No ORA results to plot.")
    if "FDR" not in result_df.columns:
        raise EnrichmentError("ORA results are missing the FDR column needed for the color scale.")
    df = result_df.copy()
    df["-log10(FDR)"] = neglog10_fdr(df["FDR"])
    select_col, select_ascending = ORA_SORT_OPTIONS.get(sort_by, ORA_SORT_OPTIONS["FDR"])
    return _enrichment_dotplot(
        df, term_col="Gene Set", x_col="Fold Enrichment", size_col="# Genes in Overlap",
        color_col="-log10(FDR)", x_label="Fold Enrichment", title="ORA Enrichment Dot Plot",
        size_legend_label="Genes in Overlap", color_legend_label="-log10(FDR)",
        top_n=top_n, reference_x=1.0, select_col=select_col, select_ascending=select_ascending,
    )
