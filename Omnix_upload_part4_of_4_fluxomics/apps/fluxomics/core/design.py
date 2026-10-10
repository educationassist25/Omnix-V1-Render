"""
core.design — Step 2: experimental design and tracer settings.
The design table is the single source of truth for which samples are analysed and how they are grouped.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict

import numpy as np
import pandas as pd

from .io_reorganize import TYPE_LAB, TYPE_REF, TYPE_BLANK

TO_HOURS = {"min": 1 / 60, "h": 1.0, "d": 24.0}


@dataclass
class Tracer:
    name: str = "[U-¹³C₆]-glucose"
    element: str = "13C"          # 13C, 15N or 2H
    labeled_atoms: int = 6        # number of labeled atoms in the tracer molecule
    purity: float = 0.99          # isotopic purity (atom fraction of the label at labeled positions)
    enrichment: float = 1.0       # fraction of the substrate pool supplied as tracer (1 = 100 % replacement)

    def to_dict(self):
        return asdict(self)


@dataclass
class CorrectionSettings:
    mode: str = "already"         # "already" (data already corrected), "apply" (correct here), "none"
    resolution: str = "high"      # "high": carbon natural abundance only; "low": all elements (nominal mass)
    method: str = "isocor_hr"     # core.correction.METHODS
    instrument: str = "orbitrap"  # AccuCor: orbitrap / fticr / tof
    resolving_power: float = 100000.0   # AccuCor: resolving power (FWHM) at ref_mz
    ref_mz: float = 200.0

    def to_dict(self):
        return asdict(self)


def default_design(S: pd.DataFrame) -> pd.DataFrame:
    """Editable design table pre-filled from the parsed sample names."""
    cols = ["Sample", "SampleType", "Group"] + (["Treatment"] if "Treatment" in S else []) + ["TimeValue", "TimeUnit"]
    D = S[cols].copy()
    if "Treatment" not in D:
        D.insert(3, "Treatment", "")
    D["Treatment"] = D.Treatment.fillna("").astype(str)
    D["Include"] = (D.SampleType == TYPE_LAB) & (S["Include"].astype(bool).to_numpy() if "Include" in S else True)
    D["Replicate"] = 0
    lab = D.SampleType == TYPE_LAB
    D.loc[lab, "Replicate"] = D[lab].groupby(["Group", "Treatment", "TimeValue"], dropna=False).cumcount() + 1
    D["TimeUnit"] = D.TimeUnit.replace("", np.nan)
    return D.reset_index(drop=True)


def time_hours(D: pd.DataFrame) -> pd.Series:
    return D.TimeValue.astype(float) * D.TimeUnit.map(TO_HOURS).fillna(1.0)


def time_label(v, unit) -> str:
    if pd.isna(v):
        return "—"
    return f"{float(v):g} {unit or 'h'}"


def finalize_design(D: pd.DataFrame) -> pd.DataFrame:
    """Add derived columns: Hours, TimeLabel, Condition; validate."""
    D = D.copy()
    if "Treatment" not in D:
        D["Treatment"] = ""
    D["Treatment"] = D.Treatment.fillna("").astype(str).str.strip()
    D["BaseGroup"] = D.Group.astype(str)
    has_t = D.Treatment != ""
    D.loc[has_t, "Group"] = D.loc[has_t, "BaseGroup"] + " + " + D.loc[has_t, "Treatment"]   # arm = group × treatment
    D["Hours"] = time_hours(D)
    D["TimeLabel"] = [time_label(v, u) for v, u in zip(D.TimeValue, D.TimeUnit)]
    D["Condition"] = np.where(D.SampleType == TYPE_LAB, D.Group.astype(str) + " · " + D.TimeLabel, D.SampleType)
    return D


def analysed(D: pd.DataFrame) -> pd.DataFrame:
    """Labeled samples that are switched on."""
    return D[(D.SampleType == TYPE_LAB) & D.Include.astype(bool)]


def groups_of(D: pd.DataFrame) -> list:
    return list(dict.fromkeys(analysed(D).Group))


def times_of(D: pd.DataFrame) -> list:
    a = analysed(D)
    return sorted(a.Hours.dropna().unique().tolist())


def samples_of(D: pd.DataFrame, group=None, hours=None) -> list:
    a = analysed(D)
    if group is not None:
        a = a[a.Group == group]
    if hours is not None:
        a = a[np.isclose(a.Hours, hours)]
    return a.Sample.tolist()


def conditions_of(D: pd.DataFrame) -> list:
    a = analysed(D).sort_values(["Group", "Hours"], key=lambda s: s.map({g: i for i, g in enumerate(groups_of(D))})
                                if s.name == "Group" else s)
    return list(dict.fromkeys(a.Condition))


def reference_samples(D: pd.DataFrame) -> list:
    return D.loc[D.SampleType == TYPE_REF, "Sample"].tolist()


def blank_samples(D: pd.DataFrame) -> list:
    return D.loc[D.SampleType == TYPE_BLANK, "Sample"].tolist()


def design_warnings(D: pd.DataFrame) -> list:
    """Plain-language checks shown under the design table."""
    w = []
    a = analysed(D)
    if a.empty:
        return ["No labeled samples are included."]
    n = a.groupby(["Group", "Hours"], dropna=False).size()
    for (g, h), k in n.items():
        if k < 3:
            w.append(f"{g} at {h:g} h has only {k} replicate(s); statistics will be unreliable (≥ 3 recommended).")
    if a.Hours.isna().any():
        w.append("Some included samples have no time point; kinetics will ignore them.")
    if len(times_of(D)) < 3:
        w.append("Fewer than 3 time points: first-order labeling kinetics will not be fitted.")
    if not reference_samples(D):
        w.append("No unlabeled (12C) reference sample: correction status cannot be checked.")
    return w


# ---------------------------------------------------------------------------
# contrasts (v2)
# ---------------------------------------------------------------------------
CONTRAST_KINDS = {
    "arms_per_time": "Arm B vs arm A at each time point",
    "all_arms_per_time": "Every pair of arms at each time point",
    "time_vs_first": "Each time point vs the first, within each arm",
    "time_consecutive": "Consecutive time points, within each arm",
    "custom": "Custom condition pairs",
}


def _ct(name, a, b, sa, sb, time=None):
    return dict(name=name, a=a, b=b, samples_a=list(sa), samples_b=list(sb), time=time)


def build_contrasts(D: pd.DataFrame, kind: str, arm_a=None, arm_b=None, pairs=None) -> list:
    """Return a list of contrasts (B vs A). D must be finalized."""
    an = analysed(D)
    arms, times = groups_of(D), times_of(D)
    lab = {h: time_label(v, u) for h, v, u in zip(an.Hours, an.TimeValue, an.TimeUnit)}
    out = []
    if kind in ("arms_per_time", "all_arms_per_time"):
        from itertools import combinations
        prs = [(arm_a, arm_b)] if kind == "arms_per_time" else list(combinations(arms, 2))
        for a, b in prs:
            for h in times:
                sa, sb = samples_of(D, a, h), samples_of(D, b, h)
                if sa and sb:
                    out.append(_ct(f"{b} vs {a} · {lab[h]}", a, b, sa, sb, lab[h]))
    elif kind in ("time_vs_first", "time_consecutive"):
        for g in arms:
            tt = [h for h in times if samples_of(D, g, h)]
            for k in range(1, len(tt)):
                h0 = tt[0] if kind == "time_vs_first" else tt[k - 1]
                h1 = tt[k]
                out.append(_ct(f"{g}: {lab[h1]} vs {lab[h0]}", f"{g} · {lab[h0]}", f"{g} · {lab[h1]}",
                               samples_of(D, g, h0), samples_of(D, g, h1), None))
    elif kind == "custom":
        for ca, cb in pairs or []:
            sa = an.loc[an.Condition == ca, "Sample"]
            sb = an.loc[an.Condition == cb, "Sample"]
            if len(sa) and len(sb) and ca != cb:
                out.append(_ct(f"{cb} vs {ca}", ca, cb, sa, sb, None))
    return out


def condition_samples(D: pd.DataFrame, conditions=None) -> dict:
    an = analysed(D)
    conds = conditions or conditions_of(D)
    return {c: an.loc[an.Condition == c, "Sample"].tolist() for c in conds}


def arms_samples(D: pd.DataFrame) -> dict:
    an = analysed(D)
    return {g: an.loc[an.Group == g, "Sample"].tolist() for g in groups_of(D)}
