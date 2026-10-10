"""
core.labeling — Step 4: isotopologue fractions and isotope-incorporation metrics.

For every metabolite (one kept peak) and sample, using ALL isotopologues M+0..M+n (n = atoms of the tracer
element in the formula), before any display filter:

  Fraction(M+i)      = Area(M+i) / Σ Area(M+0..M+n)                 (sums to 1 when the total is > 0)
  Incorporation (%)  = (1 − Fraction(M+0)) × 100                       share of molecules with ≥ 1 label
  Mean label number  = Σ i × Fraction(M+i)                             average tracer-derived labels per molecule
  Atom enrichment (%) = Mean label number / n × 100                     share of tracer-element atoms labeled

Incorporation % is NOT the percentage of carbon atoms labeled (that is the atom enrichment). Mean label number
and atom enrichment are only interpretable as tracer-derived labels after natural-abundance correction.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .correction import correct_mid, parse_formula, TRACER_ELEMENT
from .io_reorganize import RA_COLS, MAX_ISO

METRICS = {
    "incorporation": ("Total incorporation (%)", "%"),
    "mean_labels": ("Mean label number", "labels / molecule"),
    "enrichment": ("Atom enrichment (%)", "%"),
}


def n_label_atoms(formula: str, tracer: str = "13C") -> int:
    el = TRACER_ELEMENT.get(tracer, "C")
    return int(parse_formula(formula).get(el, 0))


def compute_labeling(kept_rows: pd.DataFrame, samples: list, tracer: str = "13C", correction: str = "already",
                     purity: float = 0.99, resolution: str = "high"):
    """
    Returns (fractions_long, summary_long).
      fractions_long: Compound, Method, Formula, n, Sample, i, Area, Fraction_raw, Fraction   (Fraction = corrected
                      when correction == 'apply', otherwise as provided)
      summary_long:   Compound, Method, Sample, Total, Incorporation, MeanLabels, Enrichment, SumCheck,
                      Incorporation_raw (uncorrected, for inspection)
    """
    k = kept_rows[kept_rows.Sample.isin(samples)]
    frows, srows = [], []
    for _, r in k.iterrows():
        n = min(n_label_atoms(r.Formula, tracer) or int(r.nC), MAX_ISO)
        areas = np.array([r[RA_COLS[i]] for i in range(n + 1)], float)
        areas = np.where(np.isfinite(areas), areas, 0.0)
        tot = areas.sum()
        if tot <= 0:
            continue
        raw = areas / tot
        frac = correct_mid(areas, r.Formula, tracer, purity, resolution) if correction == "apply" else raw
        idx = np.arange(n + 1)
        for i in idx:
            frows.append((r.Name, r.Mode, r.Formula, n, r.Sample, int(i), areas[i], raw[i], frac[i]))
        ml = float((idx * frac).sum())
        srows.append((r.Name, r.Mode, r.Sample, tot, (1 - frac[0]) * 100, ml, ml / n * 100 if n else np.nan,
                      float(frac.sum()), (1 - raw[0]) * 100))
    F = pd.DataFrame(frows, columns=["Compound", "Method", "Formula", "n", "Sample", "i", "Area",
                                     "Fraction_raw", "Fraction"])
    Sm = pd.DataFrame(srows, columns=["Compound", "Method", "Sample", "Total", "Incorporation", "MeanLabels",
                                      "Enrichment", "SumCheck", "Incorporation_raw"])
    return F, Sm


def metric_matrix(summary: pd.DataFrame, metric: str, samples: list) -> pd.DataFrame:
    col = {"incorporation": "Incorporation", "mean_labels": "MeanLabels", "enrichment": "Enrichment"}[metric]
    return summary.pivot(index="Compound", columns="Sample", values=col).reindex(columns=samples)


def fraction_matrix(fractions: pd.DataFrame, samples: list, channels: set | None = None,
                    max_i: int | None = None) -> pd.DataFrame:
    """Rows 'Compound | M+i'. channels: optional set of (Compound, 'M+i') to keep (e.g. the displayed ones)."""
    F = fractions.copy()
    F["Channel"] = F.Compound + " | M+" + F.i.astype(str)
    if channels is not None:
        F = F[[(c, f"M+{i}") in channels for c, i in zip(F.Compound, F.i)]]
    if max_i is not None:
        F = F[F.i <= max_i]
    return F.pivot(index="Channel", columns="Sample", values="Fraction").reindex(columns=samples)


def group_mid(fractions: pd.DataFrame, design: pd.DataFrame, compound: str) -> pd.DataFrame:
    """Mean ± SD fraction per isotopologue and condition for one metabolite."""
    F = fractions[fractions.Compound == compound].merge(design[["Sample", "Group", "Hours", "Condition"]], on="Sample")
    g = F.groupby(["Condition", "Group", "Hours", "i"]).Fraction
    return g.agg(["mean", "std", "count"]).reset_index()


def labeled_metabolites(summary: pd.DataFrame, design: pd.DataFrame, min_incorp: float = 2.0) -> list:
    """Metabolites whose mean incorporation reaches `min_incorp` % in at least one condition."""
    S = summary.merge(design[["Sample", "Condition"]], on="Sample")
    m = S.groupby(["Compound", "Condition"]).Incorporation.mean().groupby("Compound").max()
    return m[m >= min_incorp].index.tolist()


def reference_m1(kept_rows: pd.DataFrame, ref_sample: str, tracer: str = "13C") -> tuple[dict, dict]:
    """Observed M+1/M+0 ratio per metabolite in the unlabeled reference (for the correction-status check)."""
    k = kept_rows[kept_rows.Sample == ref_sample]
    obs, forms = {}, {}
    for _, r in k.iterrows():
        a0, a1 = r[RA_COLS[0]], r[RA_COLS[1]]
        if pd.notna(a0) and a0 > 0 and pd.notna(a1):
            obs[r.Name] = float(a1) / float(a0)
            forms[r.Name] = r.Formula
    return obs, forms


# ---------------------------------------------------------------------------
# v2: computation from the channel table (core.dataset), vectorized per metabolite
# ---------------------------------------------------------------------------
def compute_from_ct(CT: pd.DataFrame, samples: list, tracer: str = "13C", correction: str = "already",
                    purity: float = 0.99, resolution: str = "high", method: str | None = None,
                    resolving_power: float = 1e5, ref_mz: float = 200.0, instrument: str = "orbitrap"):
    """
    Same outputs as compute_labeling, from the channel table (rows Compound × i, columns samples).
    Missing isotopologue areas (NaN) count as 0 in the denominator: a channel that is not reported is below detection.
    Called with Dataset.labeling_table (every isotopologue M+0…M+n), so Fraction = Area(M+i) / Σ Area(M+0…M+n) —
    Compound Discoverer's "Exchange Rate in Percent".
    Returns (F long, S summary long).
    """
    samples = [s for s in samples if s in CT.columns]
    fparts, sparts = [], []
    for comp, g in CT.groupby("Compound", sort=False):
        g = g.sort_values("i")
        n = int(g.n.iloc[0])
        A = g[samples].to_numpy(float)
        A = np.where(np.isfinite(A), A, 0.0)
        A = np.clip(A, 0, None)
        tot = A.sum(axis=0)
        ok = tot > 0
        if not ok.any():
            continue
        A, tot, smp = A[:, ok], tot[ok], [s for s, o in zip(samples, ok) if o]
        raw = A / tot
        if correction == "apply" and g.Formula.iloc[0]:
            from .correction import correction_matrix, solve
            M = correction_matrix(g.Formula.iloc[0], len(g), tracer, purity, resolution, method,
                                  resolving_power, ref_mz, instrument)
            frac = np.column_stack([solve(M, raw[:, j], method) for j in range(raw.shape[1])])
        else:
            frac = raw
        idx = g.i.to_numpy()
        k = len(smp)
        fparts.append(pd.DataFrame({
            "Compound": comp, "Method": g.Method.iloc[0], "Formula": g.Formula.iloc[0], "n": n,
            "Sample": np.tile(smp, len(idx)), "i": np.repeat(idx, k), "Area": A.ravel(),
            "Fraction_raw": raw.ravel(), "Fraction": frac.ravel()}))
        ml = (idx[:, None] * frac).sum(axis=0)
        sparts.append(pd.DataFrame({
            "Compound": comp, "Method": g.Method.iloc[0], "Sample": smp, "FractionDenominator": tot,
            "Incorporation": (1 - frac[0]) * 100, "MeanLabels": ml, "Enrichment": ml / n * 100 if n else np.nan,
            "SumCheck": frac.sum(axis=0), "Incorporation_raw": (1 - raw[0]) * 100}))
    F = pd.concat(fparts, ignore_index=True) if fparts else pd.DataFrame(
        columns=["Compound", "Method", "Formula", "n", "Sample", "i", "Area", "Fraction_raw", "Fraction"])
    S = pd.concat(sparts, ignore_index=True) if sparts else pd.DataFrame(
        columns=["Compound", "Method", "Sample", "FractionDenominator", "Incorporation", "MeanLabels", "Enrichment", "SumCheck",
                 "Incorporation_raw"])
    return F, S


def reference_m1_ct(CT: pd.DataFrame, ref_sample: str) -> tuple[dict, dict]:
    """Observed M+1/M+0 ratio in the unlabeled reference, from the channel table."""
    obs, forms = {}, {}
    if ref_sample not in CT.columns:
        return obs, forms
    for comp, g in CT[CT.i <= 1].groupby("Compound", sort=False):
        v = g.set_index("i")[ref_sample]
        if 0 in v.index and 1 in v.index and pd.notna(v[0]) and v[0] > 0 and pd.notna(v[1]) and g.Formula.iloc[0]:
            obs[comp] = float(v[1]) / float(v[0])
            forms[comp] = g.Formula.iloc[0]
    return obs, forms


# ---------------------------------------------------------------------------
# complete isotope incorporation table (v2.1)
# ---------------------------------------------------------------------------
def incorporation_table(F: pd.DataFrame, S: pd.DataFrame, design: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Every metabolite: each isotopologue fraction (%), incorporation (%), fractional enrichment (%) and mean label
    number. Returns (by_sample, by_condition):
      by_sample     rows = Compound × Measure, columns = analysed samples (ordered by group, time)
      by_condition  rows = Compound × Measure, columns = <condition> mean / SD / n
    """
    an = design[(design.SampleType == "13C labeled") & design.Include.astype(bool)].copy()
    order = an.sort_values(["Group", "Hours"], key=lambda x: x.map({g: i for i, g in enumerate(dict.fromkeys(an.Group))})
                           if x.name == "Group" else x).Sample.tolist()
    Fa = F[F.Sample.isin(order)].assign(Measure=lambda d: "M+" + d.i.astype(str) + " (%)", Value=lambda d: d.Fraction * 100,
                                        _o=lambda d: d.i)
    parts = [Fa[["Compound", "Formula", "Method", "Measure", "_o", "Sample", "Value"]]]
    meta = F.drop_duplicates("Compound").set_index("Compound")[["Formula", "Method"]]
    for k, (col, name) in enumerate((("Incorporation", "Incorporation (%)"), ("Enrichment", "Fractional enrichment (%)"),
                                     ("MeanLabels", "Mean label number"))):
        X = S[S.Sample.isin(order)][["Compound", "Sample", col]].rename(columns={col: "Value"})
        X["Measure"], X["_o"] = name, 100 + k
        X = X.join(meta, on="Compound")
        parts.append(X[["Compound", "Formula", "Method", "Measure", "_o", "Sample", "Value"]])
    L = pd.concat(parts, ignore_index=True)
    corder = {c: i for i, c in enumerate(dict.fromkeys(F.Compound))}
    by_s = L.pivot_table(index=["Compound", "Formula", "Method", "Measure", "_o"], columns="Sample", values="Value",
                         aggfunc="first").reindex(columns=order).reset_index()
    by_s.columns.name = None
    by_s = by_s.assign(_c=by_s.Compound.map(corder)).sort_values(["_c", "_o"]).drop(columns=["_c", "_o"])
    L = L.merge(an[["Sample", "Condition"]], on="Sample")
    conds = list(dict.fromkeys(an.set_index("Sample").loc[order, "Condition"]))
    g = L.groupby(["Compound", "Formula", "Method", "Measure", "_o", "Condition"]).Value
    st = g.agg(["mean", "std", "count"]).unstack("Condition")
    cols = []
    for c in conds:
        for m, lab in (("mean", "mean"), ("std", "SD"), ("count", "n")):
            if (m, c) in st.columns:
                cols.append(((m, c), f"{c} {lab}"))
    by_c = st[[a for a, _ in cols]]
    by_c.columns = [b for _, b in cols]
    by_c = by_c.reset_index()
    by_c = by_c.assign(_c=by_c.Compound.map(corder)).sort_values(["_c", "_o"]).drop(columns=["_c", "_o"])
    return by_s.reset_index(drop=True), by_c.reset_index(drop=True)
