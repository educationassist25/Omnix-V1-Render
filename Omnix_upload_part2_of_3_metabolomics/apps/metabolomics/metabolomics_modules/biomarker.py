"""
biomarker.py - Biomarker discovery via statistical filtering criteria.
"""

import numpy as np
import pandas as pd


def discover_biomarkers(stats_table: pd.DataFrame, criterion: str = "combined",
                         p_cutoff: float = 0.05, fdr_cutoff: float = 0.25,
                         fc_col: str = "Log2FC", fdr_label: str = "FDR"):
    """
    criterion: 'pvalue' | 'fdr' | 'combined'
    Returns a table: Metabolite, p-value, <fdr_label>, Linear_FC, Log2FC, Direction

    fdr_label : str, default "FDR"
        Name of the Benjamini-Hochberg FDR-adjusted p-value column already
        present in `stats_table` (e.g. "BH P Value" for untargeted/unbiased
        assays). Must match the column name `stats_table` was built with.
    """
    df = stats_table.copy()

    if criterion == "pvalue":
        mask = df["p-value"] < p_cutoff
    elif criterion == "fdr":
        mask = df[fdr_label] < fdr_cutoff
    else:  # combined
        mask = (df["p-value"] < p_cutoff) & (df[fdr_label] < fdr_cutoff)

    result = df[mask].copy()
    result["Direction"] = np.where(result[fc_col] > 0, "Up", "Down")

    cols = []
    for c in ["p-value", fdr_label, "Linear_FC", fc_col, "Direction"]:
        if c in result.columns:
            cols.append(c)
    result = result[cols].rename(columns={fc_col: "Log2FC"} if fc_col != "Log2FC" else {})
    return result.sort_values("p-value")
