"""
core.stats — automatic statistics.

  2 groups  -> Welch's two-sample t-test (two-sided), per feature
               log2FC = mean(B) - mean(A) on the log2 scale;  FC = 2**log2FC
  >=3 groups -> one-way ANOVA per feature, then Tukey-Kramer HSD using the ANOVA error term (MSE) on the same
               (log2) scale.
Benjamini-Hochberg FDR is applied across features (and, for per-time-point contrasts, within each time point).
Tukey p-values are family-wise adjusted within each feature only.

The same functions work on labeling metrics (fractions, incorporation) with log_scale=False: the effect is then
the difference of means (B - A) in the metric's own units, and no fold change is reported.
"""
from __future__ import annotations

import warnings as _w
from itertools import combinations

import numpy as np
import pandas as pd
from scipy import stats as st


def bh_fdr(p) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values; NaN stays NaN."""
    p = np.asarray(p, float)
    out = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    m = ok.sum()
    if m == 0:
        return out
    pv = p[ok]
    order = np.argsort(pv)
    ranked = pv[order] * m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(ranked[::-1])[::-1]
    res = np.empty(m); res[order] = np.minimum(adj, 1.0)
    out[ok] = res
    return out


def _clean(v):
    v = np.asarray(v, float)
    return v[np.isfinite(v)]


def welch_table(X: pd.DataFrame, samples_a: list, samples_b: list, name_a: str = "A", name_b: str = "B",
                log_scale: bool = True) -> pd.DataFrame:
    """Welch t-test for every row of X. Effect = mean(B) - mean(A)."""
    rows = []
    for feat, r in X.iterrows():
        a, b = _clean(r[samples_a]), _clean(r[samples_b])
        ma, mb = (a.mean() if len(a) else np.nan), (b.mean() if len(b) else np.nan)
        sa, sb = (a.std(ddof=1) if len(a) > 1 else np.nan), (b.std(ddof=1) if len(b) > 1 else np.nan)
        warn = []
        t = p = dfw = np.nan
        if len(a) < 2 or len(b) < 2:
            warn.append(f"not tested: needs ≥ 2 values per group (n = {len(a)} vs {len(b)})")
        else:
            if len(a) < 3 or len(b) < 3:
                warn.append(f"low n ({len(a)} vs {len(b)})")
            if sa == 0 or sb == 0:
                warn.append("zero variance in a group")
            va, vb = sa ** 2 / len(a), sb ** 2 / len(b)
            if va + vb == 0:
                warn.append("not tested: no variance in either group")
            else:
                with _w.catch_warnings():
                    _w.simplefilter("ignore")
                    t, p = st.ttest_ind(b, a, equal_var=False)
                dfw = (va + vb) ** 2 / ((va ** 2 / (len(a) - 1)) + (vb ** 2 / (len(b) - 1)))
                if dfw < 2:
                    warn.append(f"very low Welch df ({dfw:.1f}): unreliable variance estimate")
                if min(sa, sb) > 0 and max(sa, sb) / min(sa, sb) > 4:
                    warn.append("variances differ > 16-fold")
        row = {"Feature": feat, f"n {name_a}": len(a), f"n {name_b}": len(b),
               f"Mean {name_a}": ma, f"Mean {name_b}": mb, f"SD {name_a}": sa, f"SD {name_b}": sb}
        if log_scale:
            row["log2FC"] = mb - ma
            row["FC"] = 2 ** (mb - ma) if np.isfinite(mb - ma) else np.nan
        else:
            row["Difference (B−A)"] = mb - ma
        row.update({"t": t, "Welch df": dfw, "p": p, "Warnings": "; ".join(warn)})
        rows.append(row)
    T = pd.DataFrame(rows)
    T.insert(T.columns.get_loc("p") + 1, "FDR (BH)", bh_fdr(T.p))
    return T


from functools import lru_cache


@lru_cache(maxsize=512)
def _qcrit(k: int, df: float) -> float:
    return float(st.studentized_range.ppf(0.95, k, df))


@lru_cache(maxsize=64)
def _sf_table(k: int, df: int):
    from scipy.interpolate import PchipInterpolator
    q = np.concatenate([np.linspace(0, 12, 241), np.linspace(12.25, 40, 112)])
    sf = np.clip(st.studentized_range.sf(q, k, df), 1e-300, 1.0)
    return PchipInterpolator(q, np.log(sf))


def studentized_sf(q, k: int, df: float):
    """Studentized-range survival function; cached interpolation for integer df (Tukey), exact otherwise."""
    q = np.atleast_1d(np.asarray(q, float))
    if float(df).is_integer() and df >= 1:
        f = _sf_table(int(k), int(df))
        out = np.exp(f(np.clip(q, 0, 40)))
        hi = q > 40
        if hi.any():
            out[hi] = st.studentized_range.sf(q[hi], k, df)
        return np.clip(out, 0, 1)
    return np.asarray(st.studentized_range.sf(q, k, df), float)


def tukey_kramer(groups: dict) -> pd.DataFrame:
    """
    Tukey-Kramer HSD for one feature. groups = {name: values}. Uses MSE from the one-way ANOVA (same scale).
    """
    g = {k: _clean(v) for k, v in groups.items()}
    g = {k: v for k, v in g.items() if len(v) >= 1}
    k = len(g)
    N = sum(len(v) for v in g.values())
    dfe = N - k
    if k < 2 or dfe < 1:
        return pd.DataFrame()
    sse = sum(((v - v.mean()) ** 2).sum() for v in g.values())
    mse = sse / dfe
    rows = []
    for (na, va), (nb, vb) in combinations(g.items(), 2):
        diff = vb.mean() - va.mean()
        se = np.sqrt(mse / 2 * (1 / len(va) + 1 / len(vb)))
        if se == 0:
            q = p = np.nan; lo = hi = np.nan
        else:
            q = abs(diff) / se
            p = float(studentized_sf(q, k, dfe)[0])
            qc = _qcrit(k, float(dfe))
            lo, hi = diff - qc * se, diff + qc * se
        rows.append({"Group A": na, "Group B": nb, "Difference (B−A)": diff, "95% CI low": lo,
                     "95% CI high": hi, "q": q, "p (Tukey)": p})
    return pd.DataFrame(rows)


def anova_table(X: pd.DataFrame, groups: dict, log_scale: bool = True):
    """
    One-way ANOVA per feature + Tukey-Kramer. groups = {group name: [samples]}.
    Returns (anova DataFrame, tukey long DataFrame).
    """
    rows, tk = [], []
    for feat, r in X.iterrows():
        vals = {gname: _clean(r[s]) for gname, s in groups.items()}
        usable = {k: v for k, v in vals.items() if len(v) >= 2}
        warn = []
        F = p = np.nan
        if len(usable) < len(vals):
            warn.append("group(s) with < 2 values excluded: " + ", ".join(k for k in vals if k not in usable))
        if len(usable) < 3:
            warn.append("not tested: fewer than 3 groups with ≥ 2 values")
        else:
            if min(len(v) for v in usable.values()) < 3:
                warn.append("low n in some groups")
            with _w.catch_warnings():
                _w.simplefilter("ignore")
                F, p = st.f_oneway(*usable.values())
                try:
                    lev_p = st.levene(*usable.values(), center="median").pvalue
                    if lev_p < 0.05:
                        warn.append(f"unequal variances (Brown-Forsythe p = {lev_p:.3f})")
                except Exception:
                    pass
            t = tukey_kramer(usable)
            if len(t):
                t.insert(0, "Feature", feat)
                tk.append(t)
        row = {"Feature": feat, **{f"Mean {k}": (v.mean() if len(v) else np.nan) for k, v in vals.items()},
               **{f"n {k}": len(v) for k, v in vals.items()}, "F": F, "p (ANOVA)": p, "Warnings": "; ".join(warn)}
        rows.append(row)
    A = pd.DataFrame(rows)
    A.insert(A.columns.get_loc("p (ANOVA)") + 1, "FDR (BH)", bh_fdr(A["p (ANOVA)"]))
    TK = pd.concat(tk, ignore_index=True) if tk else pd.DataFrame()
    if log_scale and len(TK):
        TK["log2FC (B−A)"] = TK["Difference (B−A)"]
        TK["FC"] = 2 ** TK["log2FC (B−A)"]
    return A, TK


def per_time_contrasts(X: pd.DataFrame, design_samples: pd.DataFrame, group_a: str, group_b: str,
                       log_scale: bool = True) -> pd.DataFrame:
    """Welch B vs A within each time point; FDR applied separately within each time point."""
    out = []
    for h, sub in design_samples.groupby("Hours"):
        sa = sub.loc[sub.Group == group_a, "Sample"].tolist()
        sb = sub.loc[sub.Group == group_b, "Sample"].tolist()
        if not sa or not sb:
            continue
        T = welch_table(X, sa, sb, group_a, group_b, log_scale)
        T.insert(1, "Time (h)", h)
        out.append(T)
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


def direction_text(group_a: str, group_b: str, log_scale: bool = True) -> str:
    if log_scale:
        return (f"{group_b} vs {group_a}: log2FC = mean({group_b}) − mean({group_a}); "
                f"positive = higher in {group_b}")
    return f"{group_b} vs {group_a}: difference = {group_b} − {group_a}; positive = higher in {group_b}"


# ---------------------------------------------------------------------------
# v2 additions: multiple-testing methods, contrasts, post hoc options, two-way ANOVA
# ---------------------------------------------------------------------------
ADJUST = {"bh": "Benjamini–Hochberg FDR", "by": "Benjamini–Yekutieli FDR (any dependence)",
          "holm": "Holm (family-wise)", "bonferroni": "Bonferroni (family-wise)", "none": "None (raw p)"}


def p_adjust(p, method: str = "bh") -> np.ndarray:
    """Adjusted p-values; NaN stays NaN. Matches statsmodels.multipletests."""
    p = np.asarray(p, float)
    out = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    m = int(ok.sum())
    if m == 0:
        return out
    pv = p[ok]
    if method == "none":
        out[ok] = pv
        return out
    if method == "bh":
        return bh_fdr(p)
    order = np.argsort(pv)
    ps = pv[order]
    if method == "by":
        c = np.sum(1.0 / np.arange(1, m + 1))
        adj = np.minimum.accumulate((ps * m * c / np.arange(1, m + 1))[::-1])[::-1]
    elif method == "holm":
        adj = np.maximum.accumulate(ps * (m - np.arange(m)))
    elif method == "bonferroni":
        adj = ps * m
    else:
        raise ValueError(method)
    res = np.empty(m)
    res[order] = np.minimum(adj, 1.0)
    out[ok] = res
    return out


def run_contrasts(X: pd.DataFrame, contrasts: list, log_scale: bool = True, adjust: str = "bh",
                  family: str = "contrast") -> pd.DataFrame:
    """
    Welch's t-test for every feature and every contrast.
    contrasts: list of dict(name, a, b, samples_a, samples_b, time).
    family = 'contrast' (adjust within each contrast across features) or 'all' (across every test together).
    Columns: Contrast, Time, Feature, Mean A, Mean B, effect (log2FC or Difference), t, df, p, q, Warnings.
    """
    out = []
    for c in contrasts:
        T = welch_table(X, c["samples_a"], c["samples_b"], "A", "B", log_scale)
        T = T.drop(columns=["FDR (BH)"])
        T.insert(0, "Contrast", c["name"])
        T.insert(1, "A", c["a"])
        T.insert(2, "B", c["b"])
        T.insert(3, "Time", c.get("time") or "")
        out.append(T)
    if not out:
        return pd.DataFrame()
    R = pd.concat(out, ignore_index=True)
    if family == "all":
        R["q"] = p_adjust(R.p, adjust)
    else:
        R["q"] = np.concatenate([p_adjust(g.p.to_numpy(), adjust) for _, g in R.groupby("Contrast", sort=False)]) \
            if len(R) else []
        # groupby with sort=False keeps contrast blocks in order, matching the concatenation order
    R = R.rename(columns={"q": f"q ({adjust.upper()})" if adjust != "none" else "q (none)"})
    return R


def games_howell(groups: dict) -> pd.DataFrame:
    """Games–Howell post hoc test (unequal variances, unequal n) for one feature."""
    g = {k: _clean(v) for k, v in groups.items()}
    g = {k: v for k, v in g.items() if len(v) >= 2}
    k = len(g)
    rows = []
    for (na, va), (nb, vb) in combinations(g.items(), 2):
        diff = vb.mean() - va.mean()
        sa, sb = va.var(ddof=1) / len(va), vb.var(ddof=1) / len(vb)
        se = np.sqrt(sa + sb)
        if se == 0 or k < 2:
            rows.append({"Group A": na, "Group B": nb, "Difference (B−A)": diff, "p (Games-Howell)": np.nan})
            continue
        df = (sa + sb) ** 2 / (sa ** 2 / (len(va) - 1) + sb ** 2 / (len(vb) - 1))
        q = abs(diff) / se * np.sqrt(2)
        p = float(st.studentized_range.sf(q, k, df))
        qc = _qcrit(k, round(float(df), 1)) / np.sqrt(2)      # CI: df rounded to 0.1 (p-value is exact)
        rows.append({"Group A": na, "Group B": nb, "Difference (B−A)": diff, "95% CI low": diff - qc * se,
                     "95% CI high": diff + qc * se, "q": q, "df": df, "p (Games-Howell)": min(p, 1.0)})
    return pd.DataFrame(rows)


def pairwise_welch(groups: dict, adjust: str = "holm") -> pd.DataFrame:
    """All pairwise Welch t-tests for one feature, adjusted within the feature."""
    g = {k: _clean(v) for k, v in groups.items()}
    rows = []
    for (na, va), (nb, vb) in combinations(g.items(), 2):
        if len(va) < 2 or len(vb) < 2 or (va.var(ddof=1) + vb.var(ddof=1)) == 0:
            t = p = np.nan
        else:
            with _w.catch_warnings():
                _w.simplefilter("ignore")
                t, p = st.ttest_ind(vb, va, equal_var=False)
        rows.append({"Group A": na, "Group B": nb, "Difference (B−A)": vb.mean() - va.mean() if len(va) and len(vb)
                     else np.nan, "t": t, "p": p})
    T = pd.DataFrame(rows)
    if len(T):
        T[f"p ({adjust})"] = p_adjust(T.p, adjust)
    return T


POSTHOC = {"tukey": "Tukey–Kramer HSD (equal variances)", "games_howell": "Games–Howell (unequal variances)",
           "pairwise_welch": "Pairwise Welch t-tests (Holm within feature)", "none": "None"}


def anova_v2(X: pd.DataFrame, groups: dict, log_scale: bool = True, adjust: str = "bh", posthoc: str = "tukey",
             welch_anova: bool = False, posthoc_alpha: float | None = 0.05):
    """
    One-way ANOVA (or Welch's ANOVA) per feature + chosen post hoc test. Returns (table, post hoc long).
    Post hoc tests are run for features whose adjusted ANOVA p < posthoc_alpha (None = every feature).
    """
    rows, ph, ph_jobs = [], [], []
    for feat, r in X.iterrows():
        vals = {gname: _clean(r[s]) for gname, s in groups.items()}
        usable = {k: v for k, v in vals.items() if len(v) >= 2}
        warn, F, p = [], np.nan, np.nan
        if len(usable) < len(vals):
            warn.append("group(s) with < 2 values excluded")
        if len(usable) >= 2 and sum(np.var(v) for v in usable.values()) > 0:
            with _w.catch_warnings():
                _w.simplefilter("ignore")
                if welch_anova:
                    F, p = _welch_anova(list(usable.values()))
                else:
                    F, p = st.f_oneway(*usable.values())
                try:
                    if st.levene(*usable.values(), center="median").pvalue < 0.05 and not welch_anova:
                        warn.append("unequal variances (Brown–Forsythe p < 0.05): consider Welch ANOVA / Games–Howell")
                except Exception:
                    pass
            if posthoc != "none":
                ph_jobs.append((feat, usable))
        else:
            warn.append("not tested: fewer than 2 groups with ≥ 2 values or no variance")
        rows.append({"Feature": feat, **{f"Mean {k}": (v.mean() if len(v) else np.nan) for k, v in vals.items()},
                     **{f"n {k}": len(v) for k, v in vals.items()}, "F": F, "p": p, "Warnings": "; ".join(warn)})
    A = pd.DataFrame(rows)
    A.insert(A.columns.get_loc("p") + 1, f"q ({adjust.upper()})", p_adjust(A.p, adjust))
    qmap = dict(zip(A.Feature, A[f"q ({adjust.upper()})"]))
    fn = {"tukey": tukey_kramer, "games_howell": games_howell, "pairwise_welch": pairwise_welch}.get(posthoc)
    for feat, usable in ph_jobs:
        if posthoc_alpha is not None and not (qmap.get(feat, np.nan) < posthoc_alpha):
            continue
        t = fn(usable)
        if len(t):
            t.insert(0, "Feature", feat)
            ph.append(t)
    P = pd.concat(ph, ignore_index=True) if ph else pd.DataFrame()
    if log_scale and len(P):
        P["FC"] = 2 ** P["Difference (B−A)"]
    return A, P


def _welch_anova(groups: list):
    k = len(groups)
    n = np.array([len(g) for g in groups], float)
    m = np.array([g.mean() for g in groups])
    v = np.array([g.var(ddof=1) for g in groups])
    if (v == 0).any():
        return np.nan, np.nan
    w = n / v
    mw = (w * m).sum() / w.sum()
    a = (w * (m - mw) ** 2).sum() / (k - 1)
    tmp = ((1 - w / w.sum()) ** 2 / (n - 1)).sum()
    b = 1 + 2 * (k - 2) / (k ** 2 - 1) * tmp
    F = a / b
    df2 = (k ** 2 - 1) / (3 * tmp)
    return F, float(st.f.sf(F, k - 1, df2))


def _rss(y, Xd):
    beta, *_ = np.linalg.lstsq(Xd, y, rcond=None)
    r = y - Xd @ beta
    return float(r @ r), np.linalg.matrix_rank(Xd)


def _dummies(levels, values):
    return np.column_stack([(values == l).astype(float) for l in levels[1:]]) if len(levels) > 1 \
        else np.zeros((len(values), 0))


def two_way_anova(X: pd.DataFrame, design: pd.DataFrame, factor_a: str = "Group", factor_b: str = "TimeLabel",
                  adjust: str = "bh") -> pd.DataFrame:
    """
    Two-way ANOVA with interaction (Type II sums of squares) per feature: y ~ A + B + A:B.
    design: rows = samples (column 'Sample') with the two factor columns.
    p-values are adjusted across features separately for A, B and A×B.
    """
    d = design.set_index("Sample")
    rows = []
    for feat, r in X.iterrows():
        y = r[d.index].astype(float)
        ok = y.notna().to_numpy()
        yv = y.to_numpy()[ok]
        a = d[factor_a].astype(str).to_numpy()[ok]
        b = d[factor_b].astype(str).to_numpy()[ok]
        la, lb_ = list(dict.fromkeys(a)), list(dict.fromkeys(b))
        res = {"Feature": feat, "n": int(ok.sum())}
        cells = pd.Series(list(zip(a, b))).value_counts()
        if len(la) < 2 or len(lb_) < 2 or len(cells) < len(la) * len(lb_) or len(yv) <= len(la) * len(lb_):
            res.update({"Warnings": "not tested: empty cells or no residual df"})
            rows.append(res)
            continue
        one = np.ones((len(yv), 1))
        DA, DB = _dummies(la, a), _dummies(lb_, b)
        DAB = np.column_stack([DA[:, i] * DB[:, j] for i in range(DA.shape[1]) for j in range(DB.shape[1])])
        full = np.hstack([one, DA, DB, DAB])
        rss_full, rk_full = _rss(yv, full)
        rss_ab, _ = _rss(yv, np.hstack([one, DA, DB]))
        rss_a, _ = _rss(yv, np.hstack([one, DA]))
        rss_b, _ = _rss(yv, np.hstack([one, DB]))
        dfe = len(yv) - rk_full
        mse = rss_full / dfe if dfe > 0 else np.nan
        dfa, dfb = len(la) - 1, len(lb_) - 1
        dfi = dfa * dfb
        ssa, ssb, ssi = rss_b - rss_ab, rss_a - rss_ab, rss_ab - rss_full
        for nm, ss_, df_ in (("A", ssa, dfa), ("B", ssb, dfb), ("A×B", ssi, dfi)):
            Fv = (ss_ / df_) / mse if mse and mse > 0 else np.nan
            res[f"F {nm}"] = Fv
            res[f"p {nm}"] = float(st.f.sf(Fv, df_, dfe)) if np.isfinite(Fv) else np.nan
        res["Residual df"] = dfe
        res["Warnings"] = ""
        rows.append(res)
    T = pd.DataFrame(rows)
    for nm in ("A", "B", "A×B"):
        if f"p {nm}" in T:
            T[f"q {nm}"] = p_adjust(T[f"p {nm}"], adjust)
    T = T.rename(columns=lambda c: c.replace(" A×B", f" {factor_a}×{factor_b}").replace(" A", f" {factor_a}")
                 .replace(" B", f" {factor_b}") if c.split(" ")[0] in ("F", "p", "q") else c)
    return T


# ---------------------------------------------------------------------------
# v2.1: two-group tests and ANOVA family with selectable method; p and FDR always reported
# ---------------------------------------------------------------------------
TTESTS = {"welch": "Welch's t-test (unequal variances)", "student": "Student's t-test (equal variances)",
          "paired": "Paired t-test (pairs matched by replicate number)",
          "mannwhitney": "Mann–Whitney U test (non-parametric)"}
ANOVAS = {"oneway": "One-way ANOVA (equal variances)", "welch": "Welch's ANOVA (unequal variances)",
          "kruskal": "Kruskal–Wallis test (non-parametric)", "twoway": "Two-way ANOVA (group × time)"}
POSTHOC_BY = {"oneway": ["tukey", "pairwise_welch", "none"], "welch": ["games_howell", "pairwise_welch", "none"],
              "kruskal": ["pairwise_mw", "none"]}
POSTHOC["pairwise_mw"] = "Pairwise Mann–Whitney (Holm within feature)"


def two_group(X: pd.DataFrame, samples_a: list, samples_b: list, method: str = "welch", log_scale: bool = True,
              pairs: list | None = None) -> pd.DataFrame:
    """One two-group test per feature. Effect = mean(B) − mean(A) (log2FC on log2 data). Adds FDR (BH)."""
    rows = []
    for feat, r in X.iterrows():
        if method == "paired":
            pr = [(a, b) for a, b in (pairs or []) if np.isfinite(r[a]) and np.isfinite(r[b])]
            a = np.array([r[x] for x, _ in pr], float); b = np.array([r[y] for _, y in pr], float)
        else:
            a, b = _clean(r[samples_a]), _clean(r[samples_b])
        warn, stat, p = [], np.nan, np.nan
        if len(a) < 2 or len(b) < 2:
            warn.append(f"not tested: needs ≥ 2 values per group (n = {len(a)} vs {len(b)})")
        else:
            with _w.catch_warnings():
                _w.simplefilter("ignore")
                try:
                    if method == "welch":
                        stat, p = st.ttest_ind(b, a, equal_var=False)
                    elif method == "student":
                        stat, p = st.ttest_ind(b, a, equal_var=True)
                    elif method == "paired":
                        stat, p = st.ttest_rel(b, a)
                    elif method == "mannwhitney":
                        stat, p = st.mannwhitneyu(b, a, alternative="two-sided")
                except Exception as e:  # pragma: no cover
                    warn.append(str(e))
            if np.var(a) + np.var(b) == 0:
                warn.append("no variance"); p = np.nan
            if min(len(a), len(b)) < 3:
                warn.append("low n")
        ma, mb = (a.mean() if len(a) else np.nan), (b.mean() if len(b) else np.nan)
        row = {"Feature": feat, "n A": len(a), "n B": len(b), "Mean A": ma, "Mean B": mb,
               "SD A": a.std(ddof=1) if len(a) > 1 else np.nan, "SD B": b.std(ddof=1) if len(b) > 1 else np.nan}
        if log_scale:
            row["log2FC (B−A)"] = mb - ma
            row["FC"] = 2 ** (mb - ma) if np.isfinite(mb - ma) else np.nan
        else:
            row["Difference (B−A)"] = mb - ma
        row.update({"Statistic": stat, "p": p, "Warnings": "; ".join(warn)})
        rows.append(row)
    T = pd.DataFrame(rows)
    T.insert(T.columns.get_loc("p") + 1, "FDR", bh_fdr(T.p))
    return T


def pairwise_mw(groups: dict, adjust: str = "holm") -> pd.DataFrame:
    g = {k: _clean(v) for k, v in groups.items()}
    rows = []
    for (na, va), (nb, vb) in combinations(g.items(), 2):
        p = np.nan
        if len(va) >= 2 and len(vb) >= 2:
            with _w.catch_warnings():
                _w.simplefilter("ignore")
                p = st.mannwhitneyu(vb, va, alternative="two-sided").pvalue
        rows.append({"Group A": na, "Group B": nb, "Difference of medians (B−A)": np.median(vb) - np.median(va)
                     if len(va) and len(vb) else np.nan, "p": p})
    T = pd.DataFrame(rows)
    if len(T):
        T[f"p ({adjust})"] = p_adjust(T.p, adjust)
    return T


def anova_method(X: pd.DataFrame, groups: dict, method: str = "oneway", posthoc: str = "tukey",
                 log_scale: bool = True, posthoc_alpha: float | None = 0.05, posthoc_on: str = "FDR"):
    """
    One-way family per feature: 'oneway' (F test), 'welch' (Welch's ANOVA), 'kruskal' (Kruskal–Wallis).
    Returns (table with p and FDR, post hoc long table for features with p or FDR < posthoc_alpha).
    """
    rows, jobs = [], []
    for feat, r in X.iterrows():
        vals = {k: _clean(r[s]) for k, s in groups.items()}
        usable = {k: v for k, v in vals.items() if len(v) >= 2}
        warn, stat, p = [], np.nan, np.nan
        if len(usable) < len(vals):
            warn.append("group(s) with < 2 values excluded")
        if len(usable) >= 2 and sum(np.var(v) for v in usable.values()) > 0:
            with _w.catch_warnings():
                _w.simplefilter("ignore")
                if method == "oneway":
                    stat, p = st.f_oneway(*usable.values())
                elif method == "welch":
                    stat, p = _welch_anova(list(usable.values()))
                else:
                    stat, p = st.kruskal(*usable.values())
            jobs.append((feat, usable))
        else:
            warn.append("not tested: fewer than 2 groups with ≥ 2 values or no variance")
        rows.append({"Feature": feat, **{f"Mean {k}": (v.mean() if len(v) else np.nan) for k, v in vals.items()},
                     **{f"n {k}": len(v) for k, v in vals.items()}, "Statistic": stat, "p": p,
                     "Warnings": "; ".join(warn)})
    A = pd.DataFrame(rows)
    A.insert(A.columns.get_loc("p") + 1, "FDR", bh_fdr(A.p))
    sigcol = dict(zip(A.Feature, A["FDR" if posthoc_on == "FDR" else "p"]))
    fn = {"tukey": tukey_kramer, "games_howell": games_howell, "pairwise_welch": pairwise_welch,
          "pairwise_mw": pairwise_mw}.get(posthoc)
    ph = []
    if fn is not None:
        for feat, usable in jobs:
            if posthoc_alpha is not None and not (sigcol.get(feat, np.nan) < posthoc_alpha):
                continue
            t = fn(usable)
            if len(t):
                t.insert(0, "Feature", feat)
                ph.append(t)
    P = pd.concat(ph, ignore_index=True) if ph else pd.DataFrame()
    return A, P
