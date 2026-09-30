"""
rnaseq_io.py - Data loading and validation for BulkRNAAI Pro.

Expected input formats
-----------------------
Expression matrix (CSV/XLSX), genes as rows, samples as columns:
    Gene, Sample1, Sample2, QC1, QC2, ...
Two input modes are supported (selected on the Data Upload tab):
  - "Raw Counts"   : non-negative integers (or near-integers) straight off an
                     aligner/quantifier (STAR+featureCounts, Salmon+tximport,
                     RSEM, etc.). This is the ONLY mode that supports pyDESeq2's
                     negative-binomial Wald test in the Statistics tab, since
                     DESeq2's model is defined on raw counts (it fits its own
                     size factors and dispersions -- feeding it already-
                     normalized data violates the model's assumptions).
  - "Pre-normalized logCPM" : already log2(CPM+1)-style values (e.g. from
                     edgeR/limma-voom, or exported from another pipeline).
                     pyDESeq2 cannot be used on this input (no counts to model),
                     so the Statistics tab falls back to classical
                     Welch's t-test / one-way ANOVA on the supplied values.

Metadata table (CSV/XLSX):
    Sample, Group, IsQC, Batch, <any other clinical columns>
    - Sample : must match a column name in the expression matrix
    - Group  : primary experimental/clinical group label (any name is fine --
               richer demo data uses "Diagnosis" instead, see
               get_categorical_metadata_columns in utils.py for how any
               categorical column becomes usable as a grouping variable)
    - IsQC   : True/False (or 1/0) flag for technical QC replicates
    - Batch  : optional batch identifier (usable as a DESeq2 covariate)
"""

import numpy as np
import pandas as pd

from transcriptomics_modules.utils import load_table, DataValidationError  # noqa: F401  (re-exported)

REQUIRED_METADATA_COLS = ["Sample", "Group"]


def validate_expression_matrix(df: pd.DataFrame, mode: str = "Raw Counts") -> pd.DataFrame:
    """
    Validate a genes x samples expression matrix. First column must be gene
    identifiers (symbol or Ensembl ID -- used as-is throughout; no ID mapping
    is performed here). Returns a DataFrame indexed by 'Gene' with every
    remaining column coerced to numeric.

    mode == 'Raw Counts': values are rounded and cast to non-negative
    integers (pyDESeq2 requires integer counts; fractional counts, common
    with transcript-level quantifiers like Salmon after gene-level
    summarization, are rounded rather than truncated).
    """
    if df.shape[1] < 2:
        raise DataValidationError(
            "Expression matrix must have a gene column plus at least one sample column."
        )
    first_col = df.columns[0]
    df = df.rename(columns={first_col: "Gene"})
    df["Gene"] = df["Gene"].astype(str)
    if df["Gene"].duplicated().any():
        dupes = df["Gene"][df["Gene"].duplicated()].unique().tolist()
        raise DataValidationError(f"Duplicate gene identifiers found: {dupes[:5]}...")

    sample_cols = df.columns[1:]
    for c in sample_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    out = df.set_index("Gene")

    if mode == "Raw Counts":
        if (out.dropna(how="all").fillna(0).values < 0).any():
            raise DataValidationError(
                "Negative values found in a matrix marked 'Raw Counts' — raw counts "
                "must be non-negative. If this is already-normalized/log data, "
                "reload it as 'Pre-normalized logCPM' instead."
            )
        out = out.round(0)
    return out


def validate_clinical_metadata(meta: pd.DataFrame, sample_cols) -> pd.DataFrame:
    """Validate metadata table and align it to the sample columns present in the expression matrix."""
    missing = [c for c in REQUIRED_METADATA_COLS if c not in meta.columns]
    if missing:
        raise DataValidationError(f"Metadata missing required column(s): {missing}")
    if "IsQC" not in meta.columns:
        meta["IsQC"] = False
    else:
        meta["IsQC"] = meta["IsQC"].astype(str).str.lower().isin(["true", "1", "yes", "qc"])
    if "Batch" not in meta.columns:
        meta["Batch"] = "1"
    meta = meta.set_index("Sample")
    missing_samples = [s for s in sample_cols if s not in meta.index]
    if missing_samples:
        raise DataValidationError(
            f"Samples present in expression matrix but missing from metadata: {missing_samples}"
        )
    return meta.loc[list(sample_cols)]


def split_qc_and_samples(expr_df: pd.DataFrame, meta: pd.DataFrame):
    """Return (qc_columns, sample_columns) based on metadata IsQC flag."""
    qc_cols = meta.index[meta["IsQC"]].tolist()
    sample_cols = meta.index[~meta["IsQC"]].tolist()
    return qc_cols, sample_cols


def library_size_summary(counts_df: pd.DataFrame) -> pd.DataFrame:
    """
    Per-sample library size (total mapped counts) and number of genes detected
    (count > 0). Only meaningful for raw-count input; for logCPM input the
    'library size' concept doesn't apply the same way, so callers should guard
    on mode before using this.
    """
    lib_size = counts_df.sum(axis=0)
    n_detected = (counts_df > 0).sum(axis=0)
    pct_of_median = (lib_size / lib_size.median()) * 100
    return pd.DataFrame({
        "Library_Size": lib_size,
        "Genes_Detected": n_detected,
        "Pct_of_Median_Library": pct_of_median,
    })
