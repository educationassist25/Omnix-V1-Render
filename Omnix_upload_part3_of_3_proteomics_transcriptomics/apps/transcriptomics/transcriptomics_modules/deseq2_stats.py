"""
deseq2_stats.py - "Statistics" module: pyDESeq2 Wald tests on raw counts
(preferred), matching pyDESeq2's output schema onto the same
Log2FC / p-value / FDR / Linear_FC table shape used throughout this app --
so Volcano Plot, Biomarker Discovery, Heatmap, and Boxplot all work
unmodified regardless of whether a result came from DESeq2 or the classical
t-test/ANOVA fallback used for 'Pre-normalized logCPM' input.
"""

import numpy as np
import pandas as pd

from transcriptomics_modules.utils import significance_flag

try:
    from pydeseq2.ds import DeseqStats
    from pydeseq2.default_inference import DefaultInference
    HAS_PYDESEQ2 = True
except Exception:  # pragma: no cover
    HAS_PYDESEQ2 = False


class StatsError(Exception):
    pass


def _standardize_columns(df: pd.DataFrame, level_a: str = None, level_b: str = None,
                          fdr_threshold: float = 0.05) -> pd.DataFrame:
    """Map pyDESeq2's results_df columns onto this app's shared schema.

    `Significant` is (FDR < fdr_threshold) OR (raw p-value < fdr_threshold) --
    the same single, user-adjustable threshold applied to both sides of the
    OR (see utils.significance_flag). Not an AND: requiring both cutoffs at
    once is stricter than either cutoff alone and silently drops genes that
    clear one side but not the other -- not what "flag as significant if
    either FDR or p-value clears the threshold" means.
    """
    out = pd.DataFrame({
        "baseMean": df["baseMean"],
        "Log2FC": df["log2FoldChange"],
        "lfcSE": df.get("lfcSE"),
        "stat": df.get("stat"),
        "p-value": df["pvalue"],
        "FDR": df["padj"],
    }, index=df.index)
    out["Linear_FC"] = 2.0 ** out["Log2FC"]
    out["Significant"] = significance_flag(out["p-value"], out["FDR"], fdr_threshold)
    if level_a is not None and level_b is not None:
        # DESeq2 contrast convention: level_a is the numerator, so a positive
        # Log2FC means higher expression in level_a.
        out.insert(0, "Group_A", level_a)
        out.insert(1, "Group_B", level_b)
        out["Higher_In"] = np.where(out["Log2FC"] > 0, level_a, np.where(out["Log2FC"] < 0, level_b, "Equal"))
    out.index.name = "Gene"
    return out.sort_values("p-value")


def deseq2_contrast(dds, factor: str, level_a: str, level_b: str, n_cpus: int = 1,
                     shrink_lfc: bool = True, fdr_threshold: float = 0.05) -> pd.DataFrame:
    """
    Wald test contrast: level_a vs level_b (level_a is the numerator, i.e. a
    positive Log2FC means higher in level_a). `dds` must already be fitted
    (see normalization_deseq2.fit_deseq2) with `factor` among its design
    factors. Cheap relative to fitting -- re-derives just the requested
    contrast from the already-fitted dispersions/GLM.

    shrink_lfc: pyDESeq2's own recommended follow-up step -- after the Wald
    test, shrink the log2FoldChange with an apeGLM prior (DeseqStats.
    lfc_shrink), which pulls noisy estimates for low-count/high-variance
    genes toward zero without touching the p-value/FDR. This only works when
    `dds` was fit with `level_b` as the reference level (see
    normalization_deseq2.fit_deseq2's reference_level argument) -- the
    design-matrix coefficient for level_a must exist as "<factor>[T.<level_a>]".
    If it doesn't (dds fit with a different reference), this silently falls
    back to the unshrunk MLE result rather than erroring.
    """
    if not HAS_PYDESEQ2:
        raise StatsError(
            "pydeseq2 is not installed. Run `pip install pydeseq2` to enable "
            "DESeq2-based statistics for raw-count data."
        )
    inference = DefaultInference(n_cpus=n_cpus)
    ds = DeseqStats(dds, contrast=[factor, level_a, level_b], inference=inference)
    ds.summary()
    mle_log2fc = ds.results_df["log2FoldChange"].copy()

    shrink_applied = False
    if shrink_lfc:
        coeff_name = f"{factor}[T.{level_a}]"
        if coeff_name in ds.LFC.columns:
            try:
                ds.lfc_shrink(coeff=coeff_name)
                shrink_applied = True
            except Exception:
                pass

    result = _standardize_columns(ds.results_df, level_a=level_a, level_b=level_b,
                                   fdr_threshold=fdr_threshold)
    if shrink_applied:
        # Log2FC/Linear_FC in the standardized schema are now the SHRUNK values
        # (pyDESeq2's own recommendation for ranking/plotting); the original
        # MLE estimate is kept alongside for transparency.
        result.insert(result.columns.get_loc("Log2FC") + 1, "Log2FC_MLE_Unshrunk",
                       mle_log2fc.reindex(result.index))
        result.insert(result.columns.get_loc("Linear_FC") + 1, "Linear_FC_MLE_Unshrunk",
                       (2.0 ** mle_log2fc).reindex(result.index))
    result.attrs["contrast"] = f"{level_a} vs {level_b}"
    result.attrs["factor"] = factor
    result.attrs["lfc_shrunk"] = shrink_applied
    return result


def all_pairwise_contrasts(dds, factor: str, levels: list, reference: str = None, n_cpus: int = 1,
                           fdr_threshold: float = 0.05) -> dict:
    """
    Convenience wrapper: either every level vs a single `reference` (Dunnett-
    style — the recommended default for a clinical/diagnosis factor with an
    obvious control group), or, if reference is None, every unique pair among
    `levels`. Returns {contrast_label: results_df}.
    """
    from itertools import combinations
    results = {}
    if reference is not None:
        for lvl in levels:
            if lvl == reference:
                continue
            label = f"{lvl} vs {reference}"
            results[label] = deseq2_contrast(dds, factor, lvl, reference, n_cpus=n_cpus,
                                              fdr_threshold=fdr_threshold)
    else:
        for a, b in combinations(levels, 2):
            label = f"{a} vs {b}"
            results[label] = deseq2_contrast(dds, factor, a, b, n_cpus=n_cpus,
                                              fdr_threshold=fdr_threshold)
    return results


def wald_test_logcpm(log_df: pd.DataFrame, meta: pd.DataFrame, factor: str,
                     level_a: str, level_b: str, covariates: list = None,
                     fdr_threshold: float = 0.05) -> pd.DataFrame:
    """
    Wald test for 'Pre-normalized logCPM' input, so BOTH input modes use the same
    test framework (a Wald test) and return the identical column schema.

    How this relates to the raw-counts path
    ---------------------------------------
    A Wald test is "coefficient / standard error, compared against a normal
    distribution" -- it is defined by that construction, not by any one
    likelihood. pyDESeq2 applies it to a negative-binomial GLM fitted on raw
    counts. Here there are no counts to model, so the same Wald construction is
    applied to a Gaussian linear model (OLS) fitted on the supplied logCPM
    values -- which is exactly what limma/voom-style workflows do with
    already-normalized expression values.

    So the two modes are consistent in test statistic and output, but they are
    NOT the same underlying model: raw counts get a negative-binomial GLM with
    DESeq2's dispersion shrinkage, logCPM gets an ordinary linear model with no
    shrinkage. That difference is inherent to the input, not a choice made here
    -- DESeq2's model cannot be fitted to values that are no longer counts.

    Relative to a Student's t-test on the same data, this differs only in
    referencing the normal distribution rather than a t-distribution with
    n-p degrees of freedom, so p-values are slightly anti-conservative in very
    small samples (both are otherwise identical for a simple two-group design).

    level_a is the numerator: a positive Log2FC means higher in level_a.
    Covariates (e.g. Batch) are included as additional design terms and
    controlled for, mirroring the raw-count path's design factors.
    """
    from scipy import stats as sp_stats

    covariates = covariates or []
    samples = [s for s in log_df.columns
               if s in meta.index and meta.loc[s, factor] in (level_a, level_b)]
    if len(samples) < 3:
        raise StatsError(
            f"Need at least 3 samples across '{level_a}' and '{level_b}' to fit a Wald test "
            f"(found {len(samples)})."
        )
    sub = log_df[samples]
    meta_sub = meta.loc[samples]

    # Design matrix: intercept + group indicator (1 for level_a) + covariate dummies.
    design_cols = [np.ones(len(samples)), (meta_sub[factor].values == level_a).astype(float)]
    design_names = ["Intercept", f"{level_a}_vs_{level_b}"]
    for cov in covariates:
        dummies = pd.get_dummies(meta_sub[cov].astype(str), prefix=cov, drop_first=True)
        for c in dummies.columns:
            design_cols.append(dummies[c].values.astype(float))
            design_names.append(c)
    X = np.column_stack(design_cols)
    n, p = X.shape
    if n - p < 1:
        raise StatsError(
            f"Not enough samples ({n}) for the requested design ({p} terms). "
            "Remove a covariate or include more samples."
        )

    Y = sub.to_numpy(dtype=float)                 # genes x samples
    keep = ~np.isnan(Y).any(axis=1)               # OLS below needs complete rows
    Yk = Y[keep]

    XtX_inv = np.linalg.pinv(X.T @ X)
    beta = Yk @ (XtX_inv @ X.T).T                 # genes x p
    resid = Yk - beta @ X.T
    dof = n - p
    sigma2 = (resid ** 2).sum(axis=1) / dof       # per-gene residual variance
    coef_idx = 1                                  # the group contrast term
    se = np.sqrt(sigma2 * XtX_inv[coef_idx, coef_idx])

    log2fc = beta[:, coef_idx]
    with np.errstate(divide="ignore", invalid="ignore"):
        stat = np.where(se > 0, log2fc / se, np.nan)
    pvals = 2.0 * sp_stats.norm.sf(np.abs(stat))  # Wald: normal reference

    out = pd.DataFrame(index=sub.index[keep])
    out.insert(0, "Group_A", level_a)
    out.insert(1, "Group_B", level_b)
    # baseMean on the linear CPM scale, so it is comparable to DESeq2's baseMean.
    out["baseMean"] = np.nanmean(np.power(2.0, Yk) - 1.0, axis=1).clip(min=0)
    out["Log2FC"] = log2fc
    out["lfcSE"] = se
    out["stat"] = stat
    out["p-value"] = pvals
    valid = ~np.isnan(pvals)
    fdr = np.full(len(pvals), np.nan)
    if valid.sum() > 0:
        from statsmodels.stats.multitest import multipletests
        fdr[valid] = multipletests(pvals[valid], method="fdr_bh")[1]
    out["FDR"] = fdr
    out["Linear_FC"] = 2.0 ** out["Log2FC"]
    # (FDR < threshold) OR (raw p-value < threshold), same shared threshold on
    # both sides -- see utils.significance_flag / _standardize_columns
    # docstring above for why this is an OR, not an AND.
    out["Significant"] = significance_flag(out["p-value"], out["FDR"], fdr_threshold)
    out["Higher_In"] = np.where(out["Log2FC"] > 0, level_a, np.where(out["Log2FC"] < 0, level_b, "Equal"))
    out.index.name = "Gene"
    out = out.sort_values("p-value")
    out.attrs["contrast"] = f"{level_a} vs {level_b}"
    out.attrs["factor"] = factor
    return out


def classical_two_group(log_df: pd.DataFrame, group_a_samples, group_b_samples, method: str = "ttest",
                         group_a_label: str = "Group A", group_b_label: str = "Group B",
                         fdr_threshold: float = 0.05) -> pd.DataFrame:
    """
    Fallback for 'Pre-normalized logCPM' input (no raw counts -> no DESeq2
    negative-binomial model possible): classical Welch's t-test / Mann-Whitney
    on the supplied log values, in the exact schema DESeq2 contrasts use, so
    downstream tabs are agnostic to which path produced the table.

    Also used, on the raw-count path, as the optional descriptive
    "Welch's t-test on VST/logCPM" comparison alongside the DESeq2 model-based
    result (kept as a clearly separate table -- see app.py Statistics tab).
    """
    from transcriptomics_modules import stats_analysis
    raw = stats_analysis.two_group_test(log_df, group_a_samples, group_b_samples, method=method,
                                         group_a_label=group_a_label, group_b_label=group_b_label,
                                         fdr_threshold=fdr_threshold)
    out = raw.rename(columns={"Mean_Log2_GroupA": "Mean_A", "Mean_Log2_GroupB": "Mean_B"})
    out.index.name = "Gene"
    return out


def classical_anova(log_df: pd.DataFrame, group_map: pd.Series, posthoc: str = "tukey"):
    """Fallback multi-group test for 'Pre-normalized logCPM' input."""
    from transcriptomics_modules import stats_analysis
    anova_table, posthoc_results = stats_analysis.anova_test(log_df, group_map, posthoc=posthoc)
    anova_table.index.name = "Gene"
    return anova_table, posthoc_results
