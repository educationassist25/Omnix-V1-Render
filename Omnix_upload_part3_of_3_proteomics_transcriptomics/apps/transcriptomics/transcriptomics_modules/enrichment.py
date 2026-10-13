"""
enrichment.py - Gene Set Enrichment Analysis for RNA-seq, offering three
statistically-distinct methods:

1. STRING Enrichment Analysis (Szklarczyk D et al. 2023, "The STRING database in
   2023", Nucleic Acids Research) -- hypergeometric-style functional enrichment
   computed SERVER-SIDE by STRING against STRING's own bundled/versioned
   GO/KEGG/Reactome/Pfam/InterPro annotation sets. Needs only a gene list;
   supports every organism in STRING's database (thousands of species via NCBI
   taxonomy ID), including plants. STRING's own database version (queried live via
   `get_string_version()`) governs the annotation vintage here -- NOT Enrichr's
   dated GO/KEGG library names below, which are a separate resource.

2. Over-Representation Analysis (ORA) -- a hypergeometric test computed LOCALLY
   (not via Enrichr's own /enrich endpoint, which does not expose the underlying
   universe/gene-set/comparison/overlap counts this app reports). Gene sets are
   fetched in GMT format from the appropriate Enrichr-family instance for the
   selected organism (see ORGANISM_CATALOG) and tested one at a time:
       N = |universe|, K = |gene set ∩ universe|,
       n = |comparison ∩ universe|, k = |comparison ∩ gene set ∩ universe|
       P-value = P(X >= k) via the hypergeometric survival function
       Fold Enrichment = (k/n) / (K/N)
   BH-FDR is applied across every gene set actually tested. The FULL
   comparison/input gene list is used as given -- no additional cutoff is applied
   inside this module.

3. Gene Set Enrichment Analysis, GSEA (Subramanian A et al. 2005, PNAS "Gene set
   enrichment analysis: a knowledge-based approach..."; gene-set permutation as
   implemented by GSEAPY's prerank mode, Fang Z et al. 2023, Bioinformatics) --
   ranks ALL detected genes by a signed score and computes a weighted
   running-sum enrichment statistic (ES, normalized to NES) per gene set, with
   permutation-based significance and reported leading-edge genes. This is a
   textbook-distinct algorithm from ORA/the hypergeometric test above -- it is
   never used as a substitute for Method 2, and vice versa.

All methods expect gene SYMBOLS (not raw Ensembl/RefSeq accessions) as gene
identifiers -- the convention Enrichr, STRING, and MSigDB/.gmt gene sets all
share. If your gene expression matrix uses a different identifier, convert to
gene symbol (e.g. via a standard ID mapping service, or an annotation file's
"Gene name"/"Symbol" column) before using this tab.

Database versions (verified against the live Enrichr library list; see
ORGANISM_CATALOG for the authoritative, per-organism, per-collection record --
never assume a "2026" label without verifying):
  - GO Biological/Molecular Function/Cellular Component: 2026 releases genuinely
    exist on Enrichr for HUMAN (GO_Biological_Process_2026, etc.) as of this
    writing. No organism-specific 2026 GO release exists for mouse/yeast/
    zebrafish on Enrichr; those use whatever that organism's own portal provides.
  - KEGG: a 2026 release genuinely exists for human (KEGG_2026, no species
    suffix in Enrichr's newer naming). The most recent mouse-specific KEGG
    library on Enrichr is KEGG_2019_Mouse -- explicitly NOT relabeled 2026.
  - Reactome: the most recent Enrichr release is Reactome_Pathways_2024 -- there
    is no 2026 Reactome library on Enrichr as of this writing.
  - WikiPathways: the most recent Enrichr releases are WikiPathways_2024_Human
    and WikiPathways_2024_Mouse -- no 2026 release exists.
  - MSigDB: the most recent Enrichr release is MSigDB_Hallmark_2020 -- no newer
    MSigDB library is available via Enrichr as of this writing.
"""

import io
import numpy as np
import pandas as pd
import requests
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import ticker
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

matplotlib.use("Agg")

ENRICHR_BASE = "https://maayanlab.cloud/Enrichr"
YEAST_ENRICHR_BASE = "https://maayanlab.cloud/YeastEnrichr"
FISH_ENRICHR_BASE = "https://maayanlab.cloud/FishEnrichr"
STRING_BASE = "https://string-db.org/api"


class EnrichmentError(Exception):
    pass


# ---------------------------------------------------------------------------
# Organism / database catalog
#
# Every entry below was checked against the live Enrichr library list (fetched
# through a deployed third-party tool that queries Enrichr's own
# `datasetStatistics` endpoint) rather than assumed. "library" is an exact,
# confirmed-real Enrichr library name. "library_prefix" (used only for the
# organism-specific Enrichr instances, where Enrichr does not publish a
# year-stamped name) is resolved against that instance's live library list at
# query time via resolve_library_name() -- whatever name is actually resolved
# is recorded in the run metadata for reproducibility. "available=False"
# entries are NOT silently substituted with a different database; the UI
# disables them and shows the "note" explaining why, and points to STRING
# Enrichment Analysis instead where that covers the organism.
# ---------------------------------------------------------------------------
ORGANISM_CATALOG = {
    "Human (Homo sapiens)": {
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
            "MSigDB": {"library": "MSigDB_Hallmark_2020", "version": "2020 Hallmark collection (most recent MSigDB library on Enrichr; no newer release exists)", "available": True},
        },
    },
    "Mouse (Mus musculus)": {
        "string_taxid": 10090,
        "enrichr_base": ENRICHR_BASE,
        "gene_id_type": "Gene Symbol (MGI)",
        "databases": {
            "GO Biological Process": {"library": None, "version": None, "available": False,
                                       "note": "Enrichr's GO libraries are human-gene-symbol only; no mouse-specific GO release exists on Enrichr. Use STRING Enrichment Analysis for mouse GO-style annotation."},
            "GO Molecular Function": {"library": None, "version": None, "available": False,
                                       "note": "Enrichr's GO libraries are human-gene-symbol only. Use STRING Enrichment Analysis for mouse."},
            "GO Cellular Component": {"library": None, "version": None, "available": False,
                                       "note": "Enrichr's GO libraries are human-gene-symbol only. Use STRING Enrichment Analysis for mouse."},
            "KEGG": {"library": "KEGG_2019_Mouse", "version": "2019 (most recent mouse-specific KEGG library on Enrichr; no 2026 release exists)", "available": True},
            "Reactome": {"library": None, "version": None, "available": False,
                         "note": "No mouse-specific Reactome library exists on Enrichr. Use STRING Enrichment Analysis for mouse."},
            "WikiPathways": {"library": "WikiPathways_2024_Mouse", "version": "2024 (most recent Enrichr release; no 2026 release exists)", "available": True},
            "MSigDB": {"library": None, "version": None, "available": False,
                       "note": "Enrichr's MSigDB Hallmark library is human-only. Use STRING Enrichment Analysis for mouse."},
        },
    },
    "Rat (Rattus norvegicus)": {
        "string_taxid": 10116,
        "enrichr_base": None,
        "gene_id_type": "Gene Symbol (RGD)",
        "databases": {c: {"library": None, "version": None, "available": False,
                           "note": "Enrichr has no rat-specific gene set libraries or organism portal (there is no 'RatEnrichr'). Use STRING Enrichment Analysis, which supports rat directly."}
                      for c in ["GO Biological Process", "GO Molecular Function", "GO Cellular Component",
                                "KEGG", "Reactome", "WikiPathways", "MSigDB"]},
    },
    "Yeast (S. cerevisiae)": {
        "string_taxid": 4932,
        "enrichr_base": YEAST_ENRICHR_BASE,
        "gene_id_type": "Gene Symbol (SGD)",
        "databases": {
            "GO Biological Process": {"library_prefix": "GO_Biological_Process", "version": "current YeastEnrichr release (organism-specific Enrichr instance; not year-stamped -- exact name resolved and logged at query time)", "available": True},
            "GO Molecular Function": {"library_prefix": "GO_Molecular_Function", "version": "current YeastEnrichr release (not year-stamped)", "available": True},
            "GO Cellular Component": {"library_prefix": "GO_Cellular_Component", "version": "current YeastEnrichr release (not year-stamped)", "available": True},
            "KEGG": {"library": None, "version": None, "available": False,
                     "note": "No confirmed KEGG library on YeastEnrichr. Use STRING Enrichment Analysis for yeast KEGG-style annotation."},
            "Reactome": {"library": None, "version": None, "available": False,
                         "note": "No confirmed Reactome library on YeastEnrichr."},
            "WikiPathways": {"library_prefix": "WikiPathways", "version": "current YeastEnrichr release (not year-stamped)", "available": True},
            "MSigDB": {"library": None, "version": None, "available": False, "note": "No MSigDB library on YeastEnrichr."},
        },
    },
    "Zebrafish (Danio rerio)": {
        "string_taxid": 7955,
        "enrichr_base": FISH_ENRICHR_BASE,
        "gene_id_type": "Gene Symbol (ZFIN)",
        "databases": {
            "GO Biological Process": {"library_prefix": "GO_Biological_Process", "version": "current FishEnrichr release (organism-specific Enrichr instance; not year-stamped -- exact name resolved and logged at query time)", "available": True},
            "GO Molecular Function": {"library_prefix": "GO_Molecular_Function", "version": "current FishEnrichr release (not year-stamped)", "available": True},
            "GO Cellular Component": {"library_prefix": "GO_Cellular_Component", "version": "current FishEnrichr release (not year-stamped)", "available": True},
            "KEGG": {"library": None, "version": None, "available": False,
                     "note": "No confirmed KEGG library on FishEnrichr. Use STRING Enrichment Analysis for zebrafish KEGG-style annotation."},
            "Reactome": {"library": None, "version": None, "available": False, "note": "No confirmed Reactome library on FishEnrichr."},
            "WikiPathways": {"library_prefix": "WikiPathways", "version": "current FishEnrichr release (not year-stamped)", "available": True},
            "MSigDB": {"library": None, "version": None, "available": False, "note": "No MSigDB library on FishEnrichr."},
        },
    },
}

# Major plant species: no Enrichr-family portal exists for any plant genome, so
# ORA/GSEA (Enrichr-backed) are unavailable for all of them -- STRING Enrichment
# Analysis is the supported path, since STRING itself covers these organisms.
_PLANT_SPECIES = {
    "Arabidopsis thaliana": 3702, "Rice (Oryza sativa)": 39947,
    "Maize (Zea mays)": 4577, "Soybean (Glycine max)": 3847,
    "Tomato (Solanum lycopersicum)": 4081,
}
for _name, _taxid in _PLANT_SPECIES.items():
    ORGANISM_CATALOG[f"{_name} (plant)"] = {
        "string_taxid": _taxid, "enrichr_base": None, "gene_id_type": "Gene Symbol / Locus ID",
        "databases": {c: {"library": None, "version": None, "available": False,
                           "note": "No Enrichr-family portal exists for plant species. Use STRING Enrichment Analysis, which supports this organism directly via STRING's own database."}
                      for c in ["GO Biological Process", "GO Molecular Function", "GO Cellular Component",
                                "KEGG", "Reactome", "WikiPathways", "MSigDB"]},
    }

# Kept for the STRING Enrichment Analysis tab's species selector (STRING supports
# every organism above, and many more, directly -- this dict is intentionally the
# same set surfaced for both STRING calls and the ORA/GSEA organism selector).
STRING_SPECIES = {name: cfg["string_taxid"] for name, cfg in ORGANISM_CATALOG.items()}

ORA_GSEA_DATABASE_COLLECTIONS = ["GO Biological Process", "GO Molecular Function", "GO Cellular Component",
                                  "KEGG", "Reactome", "WikiPathways", "MSigDB"]

# Display-only label overrides for gene-set collection names shown to the user
# (e.g. in selectboxes, warnings, spinners). The underlying collection name
# (the dict key used everywhere for ORGANISM_CATALOG lookups) is left
# unchanged so no internal branching logic is affected -- only what is
# rendered on screen differs. "MSigDB" is displayed as "MSigDB (Hallmark)"
# since the only MSigDB collection this app offers is the Hallmark set.
COLLECTION_DISPLAY_LABELS = {"MSigDB": "MSigDB (Hallmark)"}


def display_collection_label(collection: str) -> str:
    """User-facing label for a gene-set collection name (display only)."""
    return COLLECTION_DISPLAY_LABELS.get(collection, collection)


def get_database_info(organism: str, collection: str) -> dict:
    """Look up one organism/collection cell of ORGANISM_CATALOG."""
    org = ORGANISM_CATALOG.get(organism)
    if org is None:
        raise EnrichmentError(f"Unsupported organism: '{organism}'.")
    info = org["databases"].get(collection)
    if info is None:
        raise EnrichmentError(f"Unsupported database collection: '{collection}'.")
    return info


# ---------------------------------------------------------------------------
# Method 1: STRING Enrichment Analysis
# ---------------------------------------------------------------------------
def get_string_version():
    """Current STRING database version + stable API address (live, per-request)."""
    try:
        resp = requests.get(f"{STRING_BASE}/tsv/version", timeout=20)
        resp.raise_for_status()
        df = pd.read_csv(io.StringIO(resp.text), sep="\t")
        if len(df):
            return df.iloc[0].to_dict()
    except Exception:
        pass
    return None


def run_string_enrichment(gene_list, species=9606):
    """
    Functional enrichment via the STRING database API. Returns a DataFrame with
    columns: category (e.g. 'Process', 'KEGG', 'RCTM'), term, description,
    number_of_genes, number_of_genes_in_background, p_value, fdr, inputGenes.
    STRING computes this server-side, hypergeometric-style, against its own
    bundled annotation sets -- see get_string_version() for the exact database
    vintage behind these results.
    """
    genes = [str(g).strip() for g in gene_list if str(g).strip()]
    if len(genes) < 1:
        raise EnrichmentError("Need at least 1 gene for STRING enrichment.")
    try:
        resp = requests.post(
            f"{STRING_BASE}/tsv/enrichment",
            data={"identifiers": "%0d".join(genes), "species": species,
                  "caller_identity": "Omnix_Transcriptomics"},
            timeout=60,
        )
        resp.raise_for_status()
    except requests.exceptions.RequestException as e:
        raise EnrichmentError(
            f"Could not reach STRING (string-db.org) — check your internet connection. "
            f"Details: {e}"
        )
    text = resp.text.strip()
    if not text or text.startswith("Error"):
        return pd.DataFrame(columns=["category", "term", "description", "number_of_genes",
                                      "p_value", "fdr", "inputGenes"])
    return pd.read_csv(io.StringIO(text), sep="\t")


# ---------------------------------------------------------------------------
# Gene set library retrieval (shared by ORA and GSEA)
# ---------------------------------------------------------------------------
def get_available_libraries(enrichr_base: str = ENRICHR_BASE):
    """
    Live list of every gene set library available on a given Enrichr-family
    instance (main Enrichr, YeastEnrichr, FishEnrichr, ...), via that instance's
    `datasetStatistics` endpoint. Used to resolve organism-specific instances'
    un-dated library names (see resolve_library_name()) and to fail loudly
    rather than guess when a requested library doesn't actually exist there.
    """
    try:
        resp = requests.get(f"{enrichr_base}/datasetStatistics", timeout=30)
        resp.raise_for_status()
        data = resp.json()
    except requests.exceptions.RequestException as e:
        raise EnrichmentError(f"Could not reach {enrichr_base} — check your internet connection. Details: {e}")
    except ValueError as e:
        raise EnrichmentError(f"Unexpected response listing libraries from {enrichr_base}: {e}")
    return [entry["libraryName"] for entry in data.get("statistics", [])]


def resolve_library_name(enrichr_base: str, prefix: str) -> str:
    """
    Find the best-matching real library name on `enrichr_base` starting with
    `prefix` (case-insensitive), preferring the lexicographically LAST match
    (year-stamped names sort so the most recent year wins; un-dated names are
    returned as-is when only one match exists). Raises EnrichmentError -- never
    silently falls back to a fabricated name -- if nothing matches.
    """
    available = get_available_libraries(enrichr_base)
    matches = sorted(name for name in available if name.lower().startswith(prefix.lower()))
    if not matches:
        raise EnrichmentError(
            f"No library starting with '{prefix}' was found on {enrichr_base}. "
            f"This resource may not provide that collection for this organism."
        )
    return matches[-1]


def fetch_gene_set_library(organism: str, collection: str):
    """
    Resolve and fetch the GMT gene set library for one (organism, collection)
    cell of ORGANISM_CATALOG, for use by both ORA (Method 2) and GSEA (Method 3).

    Returns (gene_sets: {term: [gene_symbols]}, resolved_library_name: str,
    documented_version: str). Raises EnrichmentError with the catalog's
    documented reason if the combination is unsupported -- never substitutes a
    different organism's/collection's library silently.
    """
    info = get_database_info(organism, collection)
    if not info.get("available", False):
        raise EnrichmentError(
            f"'{collection}' is not available for {organism}: "
            f"{info.get('note', 'no database available for this combination.')}"
        )
    org_cfg = ORGANISM_CATALOG[organism]
    base = org_cfg["enrichr_base"]
    if base is None:
        raise EnrichmentError(f"No Enrichr-family instance configured for {organism}.")

    if info.get("library"):
        library_name = info["library"]
    else:
        library_name = resolve_library_name(base, info["library_prefix"])

    gene_sets = get_enrichr_gene_set_library(library_name, enrichr_base=base)
    return gene_sets, library_name, info["version"]


def get_enrichr_gene_set_library(library: str, enrichr_base: str = ENRICHR_BASE):
    """
    Fetch a full gene set library in GMT format from any Enrichr-family
    instance (main Enrichr by default; pass enrichr_base for YeastEnrichr/
    FishEnrichr) and parse it into {term: [gene_symbols]}.
    """
    try:
        resp = requests.get(
            f"{enrichr_base}/geneSetLibrary",
            params={"mode": "text", "libraryName": library}, timeout=60,
        )
        resp.raise_for_status()
    except requests.exceptions.RequestException as e:
        raise EnrichmentError(
            f"Could not reach {enrichr_base} — check your internet connection. Details: {e}"
        )
    if not resp.text.strip():
        raise EnrichmentError(f"Library '{library}' returned no data from {enrichr_base}.")
    return parse_gmt_text(resp.text)


def parse_gmt_text(text: str) -> dict:
    """
    Parse GMT-format text (the standard MSigDB/Enrichr convention:
    term<TAB>description<TAB>gene1<TAB>gene2<TAB>...) into {term: [gene_symbols]}.
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


# ---------------------------------------------------------------------------
# Method 2: Over-Representation Analysis (ORA) -- hypergeometric test, computed
# locally with an explicit, reported universe. This is a classic contingency-
# table enrichment test and is NOT the same algorithm as GSEA (Method 3): ORA
# takes one fixed comparison list and asks "is this gene set over-represented in
# it?"; GSEA (below) instead evaluates every gene's position in a full ranking.
# The two must never be conflated or substituted for one another.
# ---------------------------------------------------------------------------
ORA_RESULT_COLUMNS = ["Gene Set", "P-value", "FDR", "# genes in universe", "# Genes in Gene Set",
                      "# genes in comparison", "# Genes in Overlap", "Fold Enrichment", "Overlapped genes"]


def run_ora_hypergeometric(comparison_genes, gene_sets: dict, universe_genes=None,
                            min_set_size: int = 1, max_set_size: int = None):
    """
    Over-Representation Analysis via the hypergeometric test.

    comparison_genes: the FULL user-provided comparison/input gene list -- used
        exactly as given, with no additional significance/rank cutoff applied
        inside this function. (Any cutoff used to arrive at this list, e.g. an
        FDR threshold in the Statistics tab, is the caller's own upstream choice,
        made once, before calling this function -- not something this function
        does on its own.)
    gene_sets: {term: [gene_symbols]}, e.g. from fetch_gene_set_library().
    universe_genes: the background/universe gene list. If None, defaults to the
        union of every gene appearing anywhere in `gene_sets` (documented as
        such in the returned metadata). Passing the set of genes actually
        DETECTED in the experiment -- the recommended default wired in the UI --
        is generally the more defensible universe, since it reflects genes that
        could actually have been observed, not the whole genome.
    min_set_size / max_set_size: skip gene sets whose overlap with the universe
        falls outside this range (very small or very large gene sets can be
        excluded from the report; defaults keep everything with >=1 gene).

    For each gene set that clears the size filter:
        N = |universe|
        K = |gene set ∩ universe|
        n = |comparison ∩ universe|
        k = |comparison ∩ gene set ∩ universe|
        P-value = P(X >= k) = hypergeom.sf(k-1, N, K, n)   (one-sided upper tail:
            probability of seeing this many or more overlapping genes by chance)
        Fold Enrichment = (k/n) / (K/N)
    FDR = Benjamini-Hochberg, applied across every gene set actually tested here
    (not across the full library, and not silently limited to "significant" ones).

    Returns (result_df, run_meta):
        result_df has EXACTLY these columns, one row per tested gene set:
        Gene Set, P-value, FDR, # genes in universe, # Genes in Gene Set,
        # genes in comparison, # Genes in Overlap, Fold Enrichment, Overlapped genes
        (Overlapped genes are ';'-delimited.)
        run_meta is a dict of reproducibility metadata (background source and
        size, gene counts before/after universe mapping, statistical test name,
        correction method).
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
        universe_source = f"user-provided background ({len(universe)} genes -- e.g. all genes detected in this study)"
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
            pval = hypergeom.sf(k - 1, N, K, n) if k > 0 else 1.0
            fold_enrichment = (k / n) / (K / N) if (n > 0 and K > 0 and N > 0) else 0.0
            rows.append({
                "Gene Set": term, "P-value": pval, "# genes in universe": N,
                "# Genes in Gene Set": K, "# genes in comparison": n,
                "# Genes in Overlap": k, "Fold Enrichment": fold_enrichment,
                "Overlapped genes": ";".join(sorted(overlap)),
            })

    result_df = pd.DataFrame(rows, columns=[c for c in ORA_RESULT_COLUMNS if c != "FDR"])
    if len(result_df):
        result_df["FDR"] = multipletests(result_df["P-value"], method="fdr_bh")[1]
        result_df = result_df[ORA_RESULT_COLUMNS].sort_values("P-value").reset_index(drop=True)
    else:
        result_df = pd.DataFrame(columns=ORA_RESULT_COLUMNS)

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


# ---------------------------------------------------------------------------
# Method 3: Gene Set Enrichment Analysis (GSEA) -- classic weighted running-sum
# statistic (Subramanian et al. 2005). Statistically distinct from ORA (Method
# 2, above): GSEA never applies a significance cutoff to build its input --
# every detected gene is ranked and used.
# ---------------------------------------------------------------------------
def compute_ranking_score(stats_df: pd.DataFrame, method: str = "signed_neglogp") -> pd.Series:
    """
    Build the per-gene ranking score for GSEA from a Statistics-tab-style
    table (must contain a Log2FC-like column and a PValue-like column). Every
    row of `stats_df` is used -- GSEA ranks the complete detected gene list,
    not a pre-filtered significant subset.

    method:
      'signed_neglogp' : sign(log2FC) * -log10(p-value) -- the standard choice,
                          combining direction and significance (as used
                          throughout the GSEAPY prerank documentation).
      'log2fc'         : log2 fold change alone.
      'statistic'      : the test statistic column, if present (e.g. t-statistic),
                          used as-is (already signed).
    """
    lfc_col = next((c for c in stats_df.columns if c.lower() in ("log2fc", "log2_fold_change")), None)
    p_col = next((c for c in stats_df.columns if c.lower() in ("pvalue", "p_value", "p-value")), None)
    stat_col = next((c for c in stats_df.columns
                     if c.lower() in ("statistic", "stat", "wald_stat", "t-statistic", "t_statistic", "score")), None)

    if method == "statistic":
        if stat_col is None:
            raise EnrichmentError("The 'statistic' ranking needs a test-statistic column (for example the Wald "
                                  "'stat' from a DESeq2 test), and this Statistics result has none. Choose "
                                  "'signed_neglogp' or 'log2fc' instead.")
        score = stats_df[stat_col].astype(float)
    elif method == "log2fc":
        if lfc_col is None:
            raise EnrichmentError("Statistics table has no Log2FC column.")
        score = stats_df[lfc_col].astype(float)
    else:
        if lfc_col is None or p_col is None:
            raise EnrichmentError("Statistics table must contain a Log2FC and a PValue column.")
        lfc = stats_df[lfc_col].astype(float)
        pval = stats_df[p_col].astype(float).clip(lower=1e-300)
        score = np.sign(lfc) * -np.log10(pval)

    score = score.replace([np.inf, -np.inf], 0).fillna(0)
    score.index = score.index.map(lambda g: str(g).upper())
    return score


def _weighted_running_sum(ranked_scores: np.ndarray, hit_mask: np.ndarray, weight: float = 1.0):
    """
    Core Subramanian et al. 2005 weighted Kolmogorov-Smirnov running-sum statistic.
    Returns (ES, running_sum_array).
    """
    n = len(ranked_scores)
    n_hits = int(hit_mask.sum())
    if n_hits == 0:
        return 0.0, np.zeros(n)
    abs_scores = np.abs(ranked_scores) ** weight
    hit_sum = abs_scores[hit_mask].sum()
    if hit_sum == 0:
        hit_step = np.where(hit_mask, 1.0 / n_hits, 0.0)
    else:
        hit_step = np.where(hit_mask, abs_scores / hit_sum, 0.0)
    n_miss = n - n_hits
    miss_step = np.where(~hit_mask, 1.0 / max(n_miss, 1), 0.0)
    running_sum = np.cumsum(hit_step - miss_step)
    peak_idx = int(np.argmax(np.abs(running_sum)))
    es = running_sum[peak_idx]
    return es, running_sum


def run_prerank_gsea(ranked_scores: pd.Series, gene_sets: dict, weight: float = 1.0,
                      n_perm: int = 1000, min_size: int = 15, max_size: int = 500,
                      seed: int = 42):
    """
    Gene Set Enrichment Analysis (Subramanian et al. 2005), preranked mode.
    ranked_scores: Series indexed by gene symbol, EVERY detected gene (from
    compute_ranking_score() -- no significance cutoff applied). gene_sets:
    {term: [genes]}.

    Handles duplicate gene symbols in the ranking by keeping the entry with the
    largest-magnitude score per symbol (the most informative one, if the same
    symbol appears more than once after ID mapping) and drops rows with missing
    (NaN) scores before ranking, both reported in run_meta.

    Significance is estimated via GENE-SET permutation (shuffle which ranked
    positions count as "hits", same set size, recompute ES, repeat n_perm times)
    -- the same approach GSEAPY's prerank mode uses, since the original
    Subramanian et al. phenotype-permutation procedure needs the per-sample
    expression matrix and group labels, which aren't available from a ranked
    list alone. This is a genuinely different statistical test from the
    hypergeometric test used by ORA (Method 2) -- not a relabeling of it.

    Returns (result_df, run_meta). result_df: Term, Size, ES, NES, P-value, FDR
    (BH across tested gene sets), Leading_Edge_Size, Leading_Edge (leading-edge
    genes, ';'-delimited, capped at 50 for display).
    """
    if not gene_sets:
        raise EnrichmentError("No gene sets to test (the selected library returned nothing).")

    n_before = len(ranked_scores)
    ranked_scores = ranked_scores.dropna()
    n_missing = n_before - len(ranked_scores)
    if ranked_scores.index.duplicated().any():
        # Keep the largest-magnitude score per duplicated gene symbol (the most
        # informative entry, if the same symbol appears more than once after ID
        # mapping -- e.g. two transcript isoforms mapping to the same gene).
        tmp = ranked_scores.to_frame(name="score")
        tmp["abs_score"] = tmp["score"].abs()
        tmp = tmp.sort_values("abs_score", ascending=False)
        tmp = tmp[~tmp.index.duplicated(keep="first")]
        ranked_scores = tmp["score"]
    n_duplicates_collapsed = n_before - n_missing - len(ranked_scores)

    if len(ranked_scores) < min_size:
        raise EnrichmentError(
            f"Only {len(ranked_scores)} ranked genes available (after dropping missing/duplicate "
            f"entries) -- need at least {min_size} (the minimum gene set size) for GSEA to run."
        )

    ranked_scores = ranked_scores.sort_values(ascending=False)
    genes_ranked = ranked_scores.index.to_numpy()
    scores_arr = ranked_scores.to_numpy(dtype=float)
    n = len(genes_ranked)
    gene_pos = {g: i for i, g in enumerate(genes_ranked)}

    rng = np.random.default_rng(seed)
    results = []
    n_sets_too_small = 0
    n_sets_too_large = 0
    for term, term_genes in gene_sets.items():
        hit_idx = [gene_pos[g] for g in term_genes if g in gene_pos]
        size = len(hit_idx)
        if size < min_size:
            n_sets_too_small += 1
            continue
        if size > max_size:
            n_sets_too_large += 1
            continue
        hit_mask = np.zeros(n, dtype=bool)
        hit_mask[hit_idx] = True
        es, running_sum = _weighted_running_sum(scores_arr, hit_mask, weight)

        null_es = np.empty(n_perm)
        for p in range(n_perm):
            perm_idx = rng.choice(n, size=size, replace=False)
            perm_mask = np.zeros(n, dtype=bool)
            perm_mask[perm_idx] = True
            null_es[p], _ = _weighted_running_sum(scores_arr, perm_mask, weight)

        if es >= 0:
            pos_null = null_es[null_es >= 0]
            denom = pos_null.mean() if pos_null.size and pos_null.mean() != 0 else 1.0
            nes = es / denom
            pval = ((pos_null >= es).sum() + 1) / (pos_null.size + 1) if pos_null.size else 1.0
        else:
            neg_null = null_es[null_es < 0]
            denom = abs(neg_null.mean()) if neg_null.size and neg_null.mean() != 0 else 1.0
            nes = es / denom
            pval = ((neg_null <= es).sum() + 1) / (neg_null.size + 1) if neg_null.size else 1.0

        peak_idx = int(np.argmax(np.abs(running_sum)))
        if es >= 0:
            le_genes = [genes_ranked[i] for i in hit_idx if i <= peak_idx]
        else:
            le_genes = [genes_ranked[i] for i in hit_idx if i >= peak_idx]

        results.append({
            "Term": term, "Size": size, "ES": es, "NES": nes, "P-value": pval,
            "Leading_Edge_Size": len(le_genes), "Leading_Edge": ";".join(le_genes[:50]),
        })

    if not results:
        result_df = pd.DataFrame(columns=["Term", "Size", "ES", "NES", "P-value", "FDR",
                                           "Leading_Edge_Size", "Leading_Edge"])
    else:
        result_df = pd.DataFrame(results)
        result_df["FDR"] = multipletests(result_df["P-value"], method="fdr_bh")[1]
        result_df = result_df.sort_values("NES", ascending=False).reset_index(drop=True)

    run_meta = {
        "Analysis method": "Gene Set Enrichment Analysis (GSEA), preranked (Subramanian et al. 2005)",
        "Statistical test": "Weighted Kolmogorov-Smirnov running-sum statistic (ES), normalized to NES",
        "Significance estimation": f"Gene-set permutation, {n_perm} permutations per gene set",
        "Multiple-testing correction": "Benjamini-Hochberg FDR",
        "# genes in ranked list (input)": n_before,
        "# genes dropped (missing ranking value)": n_missing,
        "# duplicate gene symbols collapsed (kept largest-magnitude score)": n_duplicates_collapsed,
        "# genes in ranked list (used)": n,
        "# gene sets in library": len(gene_sets),
        "# gene sets tested": len(result_df),
        "# gene sets skipped (below min size)": n_sets_too_small,
        "# gene sets skipped (above max size)": n_sets_too_large,
        "Min/Max gene set size": f"{min_size} / {max_size}",
    }
    return result_df, run_meta




def _wrap_term(label: str, width: int = 38, max_lines: int = 2) -> str:
    """Wrap a long pathway name onto up to `max_lines` lines (the last one ellipsized if needed),
    so long names never run into the plot or get cut off by the figure edge."""
    import textwrap as _tw
    label = str(label)
    if len(label) <= width:
        return label
    # underscores (e.g. HALLMARK_...) are allowed break points, so names are split between words
    lines = [ln.replace("_ ", "_") for ln in _tw.wrap(label.replace("_", "_ "), width)]
    if len(lines) > max_lines:
        lines = lines[:max_lines]
        lines[-1] = lines[-1][: max(0, width - 1)].rstrip() + "…"
    return "\n".join(lines)


def _layout_dotplot(fig, ax, labels, font_size: float = 10.0, legend_col_in: float = 2.2):
    """Size the figure so the (wrapped) y labels fit on the left and a fixed-width legend column
    fits on the right; returns nothing (adjusts fig/ax in place)."""
    longest = max((max(len(part) for part in str(t).split("\n")) for t in labels), default=10)
    label_in = 0.0068 * font_size * longest + 0.3
    width = max(9.5, label_in + 5.2 + legend_col_in)
    fig.set_size_inches(width, fig.get_size_inches()[1])
    fig.subplots_adjust(left=label_in / width, right=1.0 - legend_col_in / width)


def _draw_dot_legends(fig, ax, reps, scale_fn, size_label, cmap, norm, color_label,
                      font_size: float = 10.0, sci_ticks: bool = False):
    """Size key and colour bar stacked in the right-hand legend column, laid out in INCHES from the
    top of the plot, so the two never overlap however tall or short the figure is."""
    import math as _m
    W, H = fig.get_size_inches()
    pos = ax.get_position()
    x0 = pos.x1 + 0.35 / W                       # legend column starts 0.35in right of the plot
    top_in = pos.y1 * H
    bottom_in = pos.y0 * H

    def fy(inches):                              # inches from bottom -> figure fraction
        return inches / H

    y = top_in
    if reps:
        fig.text(x0, fy(y), size_label, fontsize=font_size * 0.95, fontweight="bold", va="top", ha="left")
        y -= 0.30
        max_d_in = max(2.0 * _m.sqrt(scale_fn(c) / _m.pi) / 72.0 for c in reps)
        step = max_d_in + 0.10
        y -= max_d_in / 2.0
        dots_ax = fig.add_axes([0, 0, 1, 1], facecolor="none")
        dots_ax.set_xlim(0, W)
        dots_ax.set_ylim(0, H)
        dots_ax.axis("off")
        last = y
        for c in reps:
            dots_ax.scatter([x0 * W + max_d_in / 2.0], [y], s=scale_fn(c), c="#555555",
                            edgecolor="black", linewidth=0.5, clip_on=False, zorder=3)
            dots_ax.text(x0 * W + max_d_in + 0.12, y, f"{c:g}" if isinstance(c, float) else str(c),
                         fontsize=font_size * 0.9, va="center", ha="left")
            last = y
            y -= step
        y = last - max_d_in / 2.0 - 0.35         # gap between the size key and the colour bar title
    fig.text(x0, fy(y), color_label, fontsize=font_size * 0.95, fontweight="bold", va="top", ha="left")
    y -= 0.42                                    # title line + room for the top tick label
    cbar_h = max(0.6, min(1.5, y - bottom_in))
    cax = fig.add_axes([x0 + 0.02 / W, fy(y - cbar_h), 0.2 / W, fy(cbar_h)])
    sm = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cax)
    cbar.ax.tick_params(labelsize=font_size * 0.8)
    if sci_ticks:
        cbar.formatter = matplotlib.ticker.FormatStrFormatter("%.1e")
        cbar.update_ticks()
    return cbar

def enrichment_dotplot(result_df: pd.DataFrame, top_n: int = 10, sort_by: str = "FDR",
                       selection: str = "Most significant", color_by: str = "FDR",
                       cmap_name: str = "coolwarm", fig_width_in: float = 9.5,
                       fig_height_in: float = None, font_size: float = 10.0,
                       title: str = None):
    """
    Dot plot of the top enriched pathways -- the standard companion to a GSEA/
    ORA/STRING results table (and to clusterProfiler's `dotplot()`), styled to
    match the sister proteomics pipeline's enrichment dot plots: pathways
    ordered along the y-axis by their x-axis value (strongest at the top),
    dot SIZE encoding a gene-count column, dot COLOR encoding significance on a
    blue (low) -> red (high) colormap, subtle alternating row bands for
    readability, and a size legend + colorbar stacked to the right of the plot.

    Works with all three result shapes produced by this module:
      - preranked GSEA: columns Term / NES / P-value / FDR / Size (+ Leading_Edge_Size)
        -- x-axis is NES, with a dashed reference line at 0.
      - ORA:            columns Gene Set / P-value / FDR / Fold Enrichment /
        # Genes in Overlap -- x-axis is Fold Enrichment, with a dashed
        reference line at 1.0 (the "exactly as expected by chance" value).
      - STRING:         columns term / description / p_value / fdr /
        number_of_genes / number_of_genes_in_background (STRING's own,
        lowercase, TSV column names) -- normalized here into the same shape,
        with a computed Gene Ratio as the x-axis (no reference line, since a
        0-1 ratio has no single natural "no enrichment" value).

    selection:
      "Most significant" -- the top_n gene sets by `sort_by` (the usual view)
      "Top up + down"    -- the top_n/2 most positive and top_n/2 most negative
                            NES values, so an enriched-in-both-directions
                            comparison isn't shown as one-sided (GSEA only)
    """
    if result_df is None or len(result_df) == 0:
        raise EnrichmentError("No enrichment results to plot.")

    df = result_df.copy()

    # Normalize STRING's own lowercase TSV column names into the same shape
    # used by the GSEA/ORA result tables above, so this one shared plotting
    # function renders all three methods in an identical visual style rather
    # than needing a separate STRING-specific dot plot.
    if "Term" not in df.columns and "Gene Set" not in df.columns:
        if "description" in df.columns or "term" in df.columns:
            df["Gene Set"] = df["description"] if "description" in df.columns else df["term"]
        if "fdr" in df.columns and "FDR" not in df.columns:
            df["FDR"] = df["fdr"]
        if "p_value" in df.columns and "P-value" not in df.columns:
            df["P-value"] = df["p_value"]
        if ("number_of_genes" in df.columns and "number_of_genes_in_background" in df.columns
                and "Gene Ratio" not in df.columns):
            df["Gene Ratio"] = df["number_of_genes"] / df["number_of_genes_in_background"].replace(0, np.nan)
        if "number_of_genes" in df.columns and "# Genes in Gene Set" not in df.columns:
            df["# Genes in Gene Set"] = df["number_of_genes"]

    term_col = "Term" if "Term" in df.columns else ("Gene Set" if "Gene Set" in df.columns else None)
    if term_col is None:
        raise EnrichmentError("Results table has no 'Term' or 'Gene Set' column to plot.")

    has_nes = "NES" in df.columns
    if "Leading_Edge_Size" in df.columns:
        size_col, size_label = "Leading_Edge_Size", "Leading Edge Size"
    elif "Size" in df.columns:
        size_col, size_label = "Size", "Genes in set"
    elif "# Genes in Overlap" in df.columns:
        size_col, size_label = "# Genes in Overlap", "Genes in Overlap"
    elif "# Genes in Gene Set" in df.columns:
        size_col, size_label = "# Genes in Gene Set", "Gene Count"
    else:
        size_col, size_label = None, None

    if color_by not in df.columns:
        color_by = "FDR" if "FDR" in df.columns else "P-value"
    if sort_by not in df.columns:
        sort_by = "FDR" if "FDR" in df.columns else "P-value"

    # --- pick which gene sets to show ---
    if selection == "Top up + down" and has_nes:
        n_each = max(1, top_n // 2)
        up = df[df["NES"] > 0].nsmallest(n_each, sort_by)
        down = df[df["NES"] < 0].nsmallest(n_each, sort_by)
        sel = pd.concat([up, down])
    else:
        # sort_by is a p-value/FDR here, so smallest = most significant
        sel = df.nsmallest(top_n, sort_by) if sort_by in ("FDR", "P-value") \
            else df.nlargest(top_n, sort_by)

    if len(sel) == 0:
        raise EnrichmentError("No gene sets matched the current selection.")

    # x-axis: NES for GSEA (dashed reference at 0), Fold Enrichment for ORA
    # (dashed reference at 1.0, "exactly as expected by chance"), Gene Ratio
    # for STRING (no single reference value for a 0-1 ratio), else -log10(p).
    reference_x = None
    if has_nes:
        sel = sel.sort_values("NES")
        x = sel["NES"].to_numpy(dtype=float)
        x_label = "Normalized Enrichment Score (NES)"
        reference_x = 0.0
        default_title = "GSEA Enrichment Dot Plot"
    elif "Fold Enrichment" in sel.columns:
        sel = sel.sort_values("Fold Enrichment")
        x = sel["Fold Enrichment"].to_numpy(dtype=float)
        x_label = "Fold Enrichment"
        reference_x = 1.0
        default_title = "ORA Enrichment Dot Plot"
    elif "Gene Ratio" in sel.columns:
        sel = sel.sort_values("Gene Ratio")
        x = sel["Gene Ratio"].to_numpy(dtype=float)
        x_label = "Gene Ratio"
        default_title = "STRING Enrichment Dot Plot"
    else:
        sel = sel.assign(_negp=-np.log10(sel["P-value"].clip(lower=1e-300))).sort_values("_negp")
        x = sel["_negp"].to_numpy(dtype=float)
        x_label = "-log10(P-value)"
        default_title = "Pathway Dot Plot"

    terms = sel[term_col].astype(str).tolist()
    n_terms = len(terms)
    y = np.arange(n_terms)

    # dot sizes scaled into the same readable range used by the proteomics
    # pipeline's dot plots (80-500 points^2).
    if size_col:
        raw_sizes = sel[size_col].to_numpy(dtype=float)
        size_min, size_max = float(np.nanmin(raw_sizes)), float(np.nanmax(raw_sizes))
    else:
        raw_sizes = None
        size_min = size_max = None

    def _scale(v):
        if raw_sizes is None or size_max == size_min:
            return 220.0
        return 80.0 + 420.0 * (v - size_min) / (size_max - size_min)

    sizes = np.full(n_terms, 220.0) if raw_sizes is None else np.array([_scale(v) for v in raw_sizes])
    colors = sel[color_by].to_numpy(dtype=float)

    if fig_height_in is None:
        fig_height_in = max(4.6, 0.5 * n_terms + 1.3)
    fig, ax = plt.subplots(figsize=(fig_width_in, fig_height_in))

    cmap = matplotlib.colormaps.get_cmap(cmap_name)
    vmin, vmax = float(np.nanmin(colors)), float(np.nanmax(colors))
    if vmin == vmax:
        vmax = vmin + 1e-9
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)

    # Subtle alternating row bands behind the dots -- a standard readability
    # aid (clusterProfiler, ggplot2 themes) that makes it easy to trace a
    # given pathway's row across to its y-axis label, especially with many rows.
    for i in range(n_terms):
        if i % 2 == 1:
            ax.axhspan(i - 0.5, i + 0.5, color="#F4F4F4", zorder=0)

    ax.scatter(x, y, s=sizes, c=colors, cmap=cmap, norm=norm,
               edgecolor="black", linewidth=0.5, zorder=3)

    if reference_x is not None:
        ax.axvline(reference_x, color="#444444", linestyle="--", linewidth=1.0, zorder=2)

    ax.set_yticks(y)
    wrapped = [_wrap_term(t) for t in terms]
    ax.set_yticklabels(wrapped, fontsize=font_size)
    ax.set_ylim(-0.7, n_terms - 0.3)
    ax.set_xlabel(x_label, fontsize=font_size * 1.15, fontweight="bold")
    ax.tick_params(axis="x", labelsize=font_size * 0.85)
    ax.set_title(title or default_title, fontsize=font_size * 1.55, fontweight="bold", pad=14)
    ax.grid(True, axis="x", color="#E6E6E6", linewidth=0.6, zorder=1)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    _layout_dotplot(fig, ax, wrapped, font_size=font_size)
    if raw_sizes is not None and size_max > size_min:
        reps = sorted(set(int(round(v)) for v in np.linspace(size_min, size_max, num=3)))
    else:
        reps = []
    if reps:
        _draw_dot_legends(fig, ax, reps, _scale, size_label, cmap, norm, color_by, font_size=font_size,
                          sci_ticks=bool(np.nanmin(colors) < 1e-3))
    else:
        _draw_dot_legends(fig, ax, [], _scale, "", cmap, norm, color_by, font_size=font_size,
                          sci_ticks=bool(np.nanmin(colors) < 1e-3))

    return fig, sel.drop(columns=[c for c in ("_negp",) if c in sel.columns])


def running_score_plot(ranked_scores: pd.Series, gene_sets: dict, term: str, weight: float = 1.0,
                       nes: float = None, pval: float = None, fdr: float = None):
    """
    The classic GSEA 'mountain plot' for one gene set: running enrichment score
    across the ranked list (top), a barcode of hit positions (middle), and the
    ranking metric itself (bottom) -- the same 3-panel layout used by the Broad
    Institute's GSEA desktop tool and reproduced by GSEAPY's plotting functions.

    nes/pval/fdr, when supplied (from the GSEA results table for this term), are
    annotated on the plot alongside the ES so the figure is self-contained and
    interpretable without cross-referencing the results table.
    """
    if term not in gene_sets:
        raise EnrichmentError(f"Gene set '{term}' not found.")
    ranked_scores = ranked_scores.dropna().sort_values(ascending=False)
    genes_ranked = ranked_scores.index.to_numpy()
    scores_arr = ranked_scores.to_numpy(dtype=float)
    n = len(genes_ranked)
    gene_pos = {g: i for i, g in enumerate(genes_ranked)}
    hit_idx = sorted(gene_pos[g] for g in gene_sets[term] if g in gene_pos)
    if not hit_idx:
        raise EnrichmentError(f"None of '{term}'s genes were found in the ranked list.")
    hit_mask = np.zeros(n, dtype=bool)
    hit_mask[hit_idx] = True
    es, running_sum = _weighted_running_sum(scores_arr, hit_mask, weight)
    peak_idx = int(np.argmax(np.abs(running_sum)))

    fig, axes = plt.subplots(3, 1, figsize=(8, 6.4), sharex=True,
                              gridspec_kw={"height_ratios": [3, 0.6, 1.5], "hspace": 0.08})
    axes[0].plot(range(n), running_sum, color="#2E7D32", linewidth=1.6)
    axes[0].axhline(0, color="gray", linewidth=0.6)
    axes[0].axvline(peak_idx, color="black", linestyle="--", linewidth=0.8)
    axes[0].set_ylabel("Running Enrichment Score")

    def _fmt_p(v):
        if v is None or (isinstance(v, float) and np.isnan(v)):
            return None
        return f"{v:.2e}" if v < 1e-3 else f"{v:.3f}"

    # Matches the proteomics pipeline's mountain-plot header exactly: the term
    # name as the figure's overall (bold) title, and a subtitle line of
    # NES / pVal / FDR underneath -- the raw ES is not shown, since NES is the
    # standardized, comparable-across-gene-set-sizes version that a reader
    # should actually judge significance by.
    stat_bits = []
    if nes is not None and not (isinstance(nes, float) and np.isnan(nes)):
        stat_bits.append(f"NES = {nes:.3f}")
    p_txt = _fmt_p(pval)
    if p_txt:
        stat_bits.append(f"pVal = {p_txt}")
    f_txt = _fmt_p(fdr)
    if f_txt:
        stat_bits.append(f"FDR = {f_txt}")
    fig.suptitle(term, fontsize=10.5, fontweight="bold", y=0.985)
    if stat_bits:
        axes[0].set_title("   |   ".join(stat_bits), fontsize=9.5, color="#333333")

    axes[1].vlines(hit_idx, 0, 1, color="black", linewidth=0.6)
    axes[1].set_yticks([])
    axes[1].set_ylabel("Hits", fontsize=8)

    axes[2].fill_between(range(n), scores_arr, color="#4C72B0", alpha=0.7, linewidth=0)
    axes[2].axhline(0, color="gray", linewidth=0.6)
    axes[2].set_ylabel("Ranking Score")
    axes[2].set_xlabel("Rank in ordered gene list")

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    return fig
