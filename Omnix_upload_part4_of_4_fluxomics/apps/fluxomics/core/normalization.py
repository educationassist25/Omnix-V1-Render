"""
core.normalization — Step 4: log2 transformation followed by robust IQR normalization.

Two intensity datasets are produced from the (imputed) isotopologue-channel matrix, and both stay available to
every later step through a global switch:

  log2         L[f, s] = log2(area[f, s])                      (areas ≤ 0 → missing)
  log2 + IQR   N[f, s] = (L[f, s] − median_s) / IQR_s × median(IQR) + median(median_s)

median_s and IQR_s (75th − 25th percentile) are taken over the reference channels of sample s (all channels by
default, or only M+0 channels, which do not change with labeling). Re-scaling by the across-sample median IQR and
median keeps the result in log2 units, so differences between groups are still log2 fold changes.

An optional sample-amount adjustment (cell number, protein, user factors) can be applied to the areas before log2.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

DATASETS = {"log2": "log2", "log2_iqr": "log2 + IQR"}


def log2_transform(X: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    V = X.astype(float)
    bad = int(((V <= 0) & V.notna()).values.sum())
    with np.errstate(divide="ignore", invalid="ignore"):
        L = np.log2(V.where(V > 0))
    return L, bad


def sample_factors(X: pd.DataFrame, factors: dict) -> pd.DataFrame:
    """Divide each sample by a user factor (e.g. cell number relative to the median); missing factor = 1."""
    f = pd.Series({s: float(factors.get(s, 1.0) or 1.0) for s in X.columns})
    return X / f


def robust_iqr(L: pd.DataFrame, reference_rows=None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Returns (normalized matrix, per-sample table: median, IQR, scale, shift)."""
    R = L.loc[reference_rows] if reference_rows is not None else L
    med = R.median(axis=0, skipna=True)
    q1, q3 = R.quantile(0.25), R.quantile(0.75)
    iqr = (q3 - q1).replace(0, np.nan)
    target_iqr = float(np.nanmedian(iqr))
    target_med = float(np.nanmedian(med))
    N = (L - med) / iqr * target_iqr + target_med
    info = pd.DataFrame({"Median (log2)": med, "IQR (log2)": iqr, "Scale": target_iqr / iqr,
                         "Shift (log2)": target_med - med * target_iqr / iqr,
                         "Reference channels": R.notna().sum()})
    info.index.name = "Sample"
    return N, info


def distribution_long(M: pd.DataFrame, label: str) -> pd.DataFrame:
    """Long table for per-sample box plots."""
    d = M.melt(var_name="Sample", value_name="Value", ignore_index=True).dropna()
    d["Dataset"] = label
    return d


def check(L: pd.DataFrame, N: pd.DataFrame) -> dict:
    """Spread of sample medians / IQRs before and after (should shrink to ~0 after)."""
    def spread(M):
        return float(M.median().std()), float((M.quantile(.75) - M.quantile(.25)).std())
    a, b = spread(L), spread(N)
    return {"SD of sample medians (log2)": a[0], "SD of sample medians after": b[0],
            "SD of sample IQRs (log2)": a[1], "SD of sample IQRs after": b[1]}
