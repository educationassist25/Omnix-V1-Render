"""
core.kinetics — Step 7: isotope-label arrival and APPARENT labeling kinetics.

All quantities here describe how fast label appears in a metabolite pool. They are NOT absolute pathway
fluxes: labeling speed also depends on pool size, compartmentation, exchange and the label supply.

  * Summary per metabolite, group and time: n, mean, SD, SEM, 95 % CI (t distribution).
  * Apparent slope between consecutive measured time points: Δ metric / Δ time (e.g. % points per hour).
  * Optional first-order model  E(t) = E∞ (1 − exp(−k t)),  t½ = ln 2 / k
      fitted only when ≥ 3 distinct time points exist AND the group means do not fall over time;
      reports estimates, 95 % CIs, R², RMSE, AICc vs a straight line through the origin, convergence,
      and flags plateau extrapolation.
  * Speed ranking only among supported fits (never from a single time point).
"""
from __future__ import annotations

import warnings as _w

import numpy as np
import pandas as pd
from scipy import stats as st
from scipy.optimize import curve_fit

LN2 = np.log(2)


def summarize(values: pd.DataFrame, value_col: str = "Value", extra_keys: tuple = ()) -> pd.DataFrame:
    """values: Compound, Group, Hours, Sample, Value -> mean/SD/SEM/CI per Compound, Group, Hours (+ extra keys)."""
    keys = ["Compound", *extra_keys, "Group", "Hours"]
    g = values.dropna(subset=[value_col]).groupby(keys)[value_col]
    S = g.agg(n="count", mean="mean", sd="std").reset_index()
    S["sem"] = S["sd"] / np.sqrt(S["n"])
    tcrit = st.t.ppf(0.975, np.maximum(S["n"] - 1, 1))
    S["ci95"] = np.where(S["n"] > 1, tcrit * S["sem"], np.nan)
    return S


def apparent_slopes(summary: pd.DataFrame, unit: str = "% points/h", assume_zero_at_t0: bool = False) -> pd.DataFrame:
    """Slope between consecutive measured time points (and optionally from t = 0 assuming no label)."""
    rows = []
    for (c, g), d in summary.groupby(["Compound", "Group"]):
        d = d.sort_values("Hours")
        pts = list(zip(d.Hours, d["mean"], d["sem"].fillna(0)))
        if assume_zero_at_t0 and len(pts) and pts[0][0] > 0:
            pts = [(0.0, 0.0, 0.0)] + pts
        for (t1, m1, s1), (t2, m2, s2) in zip(pts[:-1], pts[1:]):
            dt = t2 - t1
            if dt <= 0:
                continue
            rows.append({"Compound": c, "Group": g, "From (h)": t1, "To (h)": t2,
                         f"Apparent slope ({unit})": (m2 - m1) / dt, "SE": np.sqrt(s1 ** 2 + s2 ** 2) / dt,
                         "Assumes no label at t = 0": bool(assume_zero_at_t0 and t1 == 0)})
    return pd.DataFrame(rows)


def first_order(t, e_inf, k):
    return e_inf * (1 - np.exp(-k * t))


def _aicc(rss, n, p):
    if n - p - 1 <= 0 or rss <= 0:
        return np.nan
    return n * np.log(rss / n) + 2 * p + 2 * p * (p + 1) / (n - p - 1)


def fit_one(t: np.ndarray, y: np.ndarray, means: pd.Series, ymax_bound: float = 100.0) -> dict:
    """Fit E(t) = E∞(1 − e^{−kt}) to replicate-level points. `means` = group mean per distinct time."""
    out = dict(status="", reason="", n_points=len(t), n_times=int(len(np.unique(t))),
               E_inf=np.nan, E_inf_lo=np.nan, E_inf_hi=np.nan, k=np.nan, k_lo=np.nan, k_hi=np.nan,
               t_half=np.nan, t_half_lo=np.nan, t_half_hi=np.nan, R2=np.nan, RMSE=np.nan,
               dAICc_vs_linear=np.nan, converged=False, flags="")
    ok = np.isfinite(t) & np.isfinite(y)
    t, y = t[ok], y[ok]
    nt = len(np.unique(t))
    out["n_points"], out["n_times"] = len(t), nt
    if nt < 3:
        out.update(status="Not fitted", reason=f"only {nt} distinct time point(s); ≥ 3 required")
        return out
    m = means.sort_index()
    drops = [(m.index[i], m.index[i + 1]) for i in range(len(m) - 1)
             if m.iloc[i + 1] < m.iloc[i] - max(0.05 * abs(m.iloc[i]), 0.5)]
    if drops:
        out.update(status="Not fitted",
                   reason="labeling falls between " + ", ".join(f"{a:g}→{b:g} h" for a, b in drops) +
                          " — a first-order rise to plateau does not describe these data")
        return out
    tmax = t.max()
    p0 = [min(max(y.max() * 1.2, 1e-3), ymax_bound), LN2 / np.median(t[t > 0]) if (t > 0).any() else 0.1]
    try:
        with _w.catch_warnings():
            _w.simplefilter("ignore")
            popt, pcov = curve_fit(first_order, t, y, p0=p0, bounds=([0, 1e-6], [ymax_bound, 50.0]), maxfev=20000)
        out["converged"] = True
    except Exception as e:
        out.update(status="Not fitted", reason=f"did not converge ({type(e).__name__})")
        return out
    e_inf, k = popt
    resid = y - first_order(t, *popt)
    rss = float((resid ** 2).sum())
    tss = float(((y - y.mean()) ** 2).sum())
    dof = max(len(t) - 2, 1)
    se = np.sqrt(np.clip(np.diag(pcov), 0, np.inf)) if np.all(np.isfinite(pcov)) else [np.nan, np.nan]
    tc = st.t.ppf(0.975, dof)
    out.update(E_inf=e_inf, k=k, t_half=LN2 / k,
               E_inf_lo=e_inf - tc * se[0], E_inf_hi=e_inf + tc * se[0],
               k_lo=max(k - tc * se[1], 0), k_hi=k + tc * se[1],
               R2=1 - rss / tss if tss > 0 else np.nan, RMSE=np.sqrt(rss / len(t)))
    out["t_half_lo"] = LN2 / out["k_hi"] if out["k_hi"] > 0 else np.nan
    out["t_half_hi"] = LN2 / out["k_lo"] if out["k_lo"] > 0 else np.inf
    # straight line through the origin as the simpler alternative
    a = float((t * y).sum() / (t * t).sum()) if (t * t).sum() > 0 else 0.0
    rss_lin = float(((y - a * t) ** 2).sum())
    out["dAICc_vs_linear"] = _aicc(rss, len(t), 2) - _aicc(rss_lin, len(t), 1)
    flags = []
    if nt == 3:
        flags.append("only 3 time points: 1 degree of freedom per curve beyond the 2 parameters")
    if e_inf > 1.15 * m.max() or out["t_half"] > tmax:
        flags.append("plateau E∞ extrapolated beyond the observed data")
    if not np.isfinite(se).all() or (np.isfinite(out["t_half_hi"]) and out["t_half_hi"] / max(out["t_half"], 1e-9) > 5) \
            or not np.isfinite(out["t_half_hi"]):
        flags.append("parameters poorly constrained (wide CI)")
    if e_inf >= ymax_bound * 0.999 or k >= 49.9:
        flags.append("estimate at parameter bound")
    if np.isfinite(out["R2"]) and out["R2"] < 0.6:
        flags.append(f"poor fit (R² = {out['R2']:.2f})")
    if np.isfinite(out["dAICc_vs_linear"]) and out["dAICc_vs_linear"] > 2:
        flags.append("straight line fits better (ΔAICc > 2)")
    out["flags"] = "; ".join(flags)
    serious = [f for f in flags if not f.startswith("only 3")]
    out["status"] = "Fitted" if not serious else "Fitted with warnings"
    return out


def fit_all(values: pd.DataFrame, value_col: str = "Value", ymax_bound: float = 100.0) -> pd.DataFrame:
    """Fit every Compound × Group. values: Compound, Group, Hours, Sample, Value (replicate level)."""
    rows = []
    for (c, g), d in values.dropna(subset=["Hours", value_col]).groupby(["Compound", "Group"]):
        means = d.groupby("Hours")[value_col].mean()
        r = fit_one(d.Hours.to_numpy(float), d[value_col].to_numpy(float), means, ymax_bound)
        r.update(Compound=c, Group=g, n_replicates=int(d.groupby("Hours").size().min()))
        rows.append(r)
    cols = ["Compound", "Group", "status", "reason", "flags", "E_inf", "E_inf_lo", "E_inf_hi", "k", "k_lo", "k_hi",
            "t_half", "t_half_lo", "t_half_hi", "R2", "RMSE", "dAICc_vs_linear", "converged", "n_times",
            "n_points", "n_replicates"]
    return pd.DataFrame(rows)[cols] if rows else pd.DataFrame(columns=cols)


def speed_ranking(fits: pd.DataFrame) -> tuple[pd.DataFrame, str]:
    """Rank by apparent t½ only among fits with status 'Fitted' (no serious flags)."""
    ok = fits[(fits.status == "Fitted") & np.isfinite(fits.t_half)]
    if ok.empty:
        return ok, ("No ranking: no metabolite has a supported first-order fit. Ranking by speed is not "
                    "defensible from these data (too few time points, non-monotonic labeling or poorly "
                    "constrained parameters).")
    R = ok.sort_values("t_half")[["Compound", "Group", "t_half", "t_half_lo", "t_half_hi", "E_inf", "R2"]].copy()
    R.insert(0, "Rank", R.groupby("Group").cumcount() + 1)
    return R, f"Ranked {len(R)} metabolite × group curves with supported first-order fits by apparent t½."
