"""
core.io_reorganize — Step 1: data input and organization.

Reads one or more Compound Discoverer SIL exports (one row per compound x raw file, isotopologue areas in
"Relative Area 0..25") and builds the final reorganized data set with ONE peak per metabolite:
  * isotopologue rows M+0..M+n (n = carbons), only rows non-zero in every labeled sample are displayed
  * peaks whose M+0 is 0 in any labeled sample are removed (unreliable)
  * candidates (all peaks / names / methods of the same metabolite) are ranked:
      1 not a likely in-source fragment (co-elutes with larger compound = this + neutral loss)
      2 detected in more labeled samples   3 more labeled isotopologues
      4 higher signal (>= SIGNAL_FOLD x)   5 fewer CD warnings   6 higher signal
  * v2.2: a metabolite whose kept peak has ONLY the M+0 isotopologue row is excluded from the combined table
No values are changed: every number in the workbook is copied from the source files.

v2.2: sample naming, name matching and selection follow sil_flux_pipeline.py stage 1 EXACTLY (FILE_PATTERN,
GROUP_ORDER, TIME_ORDER, SHORT_GROUP, ALIASES only — no extra library synonyms, no second zero rule). The generic
file-name parser is used only when FILE_PATTERN does not recognise the raw-file names.
"""
from __future__ import annotations

import io
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter
from scipy.optimize import linear_sum_assignment

MAX_ISO = 25
RA_COLS = [f"Relative Area {i}" for i in range(MAX_ISO + 1)]
REQUIRED_COLS = ["Name", "Formula", "Molecular Weight", "RT in min", "Area", "Status",
                 "Study File ID", "Filename"] + RA_COLS

TYPE_REF, TYPE_BLANK, TYPE_LAB = "12C reference", "Blank", "13C labeled"
STEREO_TOKENS = {"l", "d", "dl", "r", "s", "rs"}
MODE_PATTERN = re.compile(r"(RP|HILIC|C18|ZIC|IC|GC)[_\- ]?(Pos|Neg)", re.I)
TIME_TOKEN = re.compile(r"^(\d+(?:\.\d+)?)\s*(h|hr|hrs|hour|hours|min|mins|minute|minutes|m|d|day|days)$", re.I)
REF_TOKEN = re.compile(r"^(12C|C12|unlab|unlabeled|unlabelled|natural|nolabel|nonlabeled|U12C)$", re.I)
BLANK_TOKEN = re.compile(r"^(blank|blk|solvent)", re.I)
CONTROL_HINT = re.compile(r"(ctrl|control|vehicle|wt|wildtype|scramble|shc|sicontrol|mock)", re.I)

# Names that differ between methods but are the same metabolite (name -> name used in the other file)
ALIASES = {
    "Inosinic acid": "Inosine-5'-monophosphate (IMP)",
}


def _raw_key(name: str) -> str:
    tokens = [t for t in re.split(r"[^0-9a-zA-Zα-ωΑ-Ω]+", str(name).lower()) if t and t not in STEREO_TOKENS]
    return " ".join(sorted(tokens))


def _library_keymap() -> dict:
    """key(synonym) -> key(standard name) for every curated library entry, so synonyms count as one metabolite."""
    f = Path(__file__).resolve().parent.parent / "data" / "metabolite_library.csv"
    out = {}
    if f.exists():
        L = pd.read_csv(f, dtype=str).fillna("")
        for _, r in L.iterrows():
            std = _raw_key(r.StandardName)
            for nm in [r.StandardName] + [x for x in r.Synonyms.split("|") if x]:
                out.setdefault(_raw_key(nm), std)
    return out


# v2.2: the curated-library synonym map is NOT used for peak selection (sil_flux_pipeline.py uses ALIASES only).
LIBRARY_KEYMAP: dict = {}
# Force the method for specific metabolites, e.g. {"L-Glutamic acid": "RP_Pos"}
PREFER: dict = {}
SIGNAL_FOLD = 1.5  # signal must differ by at least this factor to decide at rule 4


def mode_of(name: str) -> str:
    """Acquisition method from a file name, e.g. 'RP_Pos', 'HILIC_Neg'."""
    m = MODE_PATTERN.search(name)
    return f"{m.group(1).upper()}_{m.group(2).capitalize()}" if m else Path(name).stem


# ---------------------------------------------------------------------------
# Sample naming — EXACT copy of sil_flux_pipeline.py
# ---------------------------------------------------------------------------
# Filenames look like:
#   08302026_RP_Pos_CellLine_GeneKD_1001_12Hr.raw   -> labeled sample
#   08302026_RP_Pos_CellLine_shControl_1009_12C.raw -> 12C (unlabeled) reference
#   08302026_RP_Pos_Blank_3.raw                     -> blank
# Edit FILE_PATTERN / GROUP_ORDER / TIME_ORDER / SHORT_GROUP for other studies.
FILE_PATTERN = re.compile(r"_(?P<group>GeneKD|shControl)_(?P<num>\d+)_(?P<time>\w+?)\.raw$", re.I)
GROUP_ORDER = ["shControl", "GeneKD"]
TIME_ORDER = ["12 h", "18 h", "40 h"]
SHORT_GROUP = {"shControl": "shCtrl", "GeneKD": "GeneKD"}


def parse_filename(fn: str) -> dict:
    """Return group / time / sample number / type / short sample name for one raw file."""
    if "blank" in fn.lower():
        tag = re.search(r"Blank_?(\w+)\.raw", fn, re.I)
        name = "Blank" + (f"_{tag.group(1)}" if tag else "")
        return dict(Group="Blank", Time="—", SampleNo="—", SampleType=TYPE_BLANK, Sample=name)
    m = FILE_PATTERN.search(fn)
    if not m:  # unknown file: keep it, put it at the end
        stem = Path(fn).stem
        return dict(Group="Other", Time="—", SampleNo="—", SampleType="Other", Sample=stem)
    group, num, time = m.group("group"), m.group("num"), m.group("time")
    if time.upper() == "12C":
        return dict(Group=group, Time="12C (unlabeled)", SampleNo=num,
                    SampleType=TYPE_REF, Sample=f"12C_Ref_{num}")
    time = re.sub(r"(?i)hr?s?$", " h", time).strip()  # "12Hr" -> "12 h"
    short = SHORT_GROUP.get(group, group)
    return dict(Group=group, Time=time, SampleNo=num, SampleType=TYPE_LAB,
                Sample=f"{short}_{time.replace(' ', '')}_{num}")


def _ref_sort_key(r):
    type_rank = {TYPE_REF: 0, TYPE_BLANK: 1, TYPE_LAB: 2, "Other": 3}
    g = GROUP_ORDER.index(r.Group) if r.Group in GROUP_ORDER else 99
    t = TIME_ORDER.index(r.Time) if r.Time in TIME_ORDER else 99
    return (type_rank.get(r.SampleType, 3), g, t, r.SampleNo)


def _add_time_columns(S: pd.DataFrame) -> pd.DataFrame:
    """TimeValue / TimeUnit for the design table, read from the reference Time label ('12 h' -> 12, 'h')."""
    tv, tu = [], []
    for t in S.Time.astype(str):
        m = re.fullmatch(r"(\d+(?:\.\d+)?)\s*(h|min|d)", t.strip())
        tv.append(float(m.group(1)) if m else np.nan)
        tu.append(m.group(2) if m else "")
    return S.assign(TimeValue=tv, TimeUnit=tu)


def reference_recognises(filenames) -> bool:
    """True when every non-blank raw file matches FILE_PATTERN (the reference naming applies)."""
    files = [f for f in dict.fromkeys(filenames) if "blank" not in str(f).lower()]
    return bool(files) and all(FILE_PATTERN.search(str(f)) for f in files)


def build_sample_table_reference(df: pd.DataFrame) -> pd.DataFrame:
    """sil_flux_pipeline.build_sample_table, unchanged (plus TimeValue / TimeUnit for the design step)."""
    files = df[["Filename", "Study File ID"]].drop_duplicates()
    rows = [dict(Filename=f, StudyFileID=i, **parse_filename(f)) for f, i in files.values]
    S = pd.DataFrame(rows)
    S = S.iloc[sorted(range(len(S)), key=lambda i: _ref_sort_key(S.iloc[i]))].reset_index(drop=True)
    if (S.SampleType == TYPE_BLANK).sum() == 1:  # single blank -> just "Blank"
        S.loc[S.SampleType == TYPE_BLANK, "Sample"] = "Blank"
    S["Block"] = np.where(S.SampleType == TYPE_LAB, S.Group + " · " + S.Time, S.SampleType)
    if S.Sample.duplicated().any():
        raise ValueError("Two raw files map to the same sample name; check FILE_PATTERN.")
    return _add_time_columns(S)


# ---------------------------------------------------------------------------
# Generic sample-name parsing (fallback for file names FILE_PATTERN does not recognise)
# ---------------------------------------------------------------------------
def _tokens(fn: str):
    stem = re.sub(r"\.(raw|d|mzml|mzxml|wiff)$", "", Path(fn).name, flags=re.I)
    return [t for t in re.split(r"[_\s]+", stem) if t]


def _unit_norm(u: str) -> str:
    u = u.lower()
    if u.startswith("h"):
        return "h"
    if u.startswith("d"):
        return "d"
    return "min"


def parse_filenames(filenames) -> pd.DataFrame:
    """
    Parse raw-file names into Group / Time / SampleNo / SampleType / Sample without study-specific patterns.
    Tokens identical in every non-blank file of a method (date, column, polarity, cell line) are ignored;
    a token like 12Hr / 30min / 2d is the time; a token like 12C / unlabeled marks the unlabeled reference;
    a numeric token is the sample number; the remaining tokens form the group name.
    """
    files = list(dict.fromkeys(filenames))
    toks = {f: _tokens(f) for f in files}
    nonblank = [f for f in files if not any(BLANK_TOKEN.match(t) for t in toks[f])]
    common = set.intersection(*[set(toks[f]) for f in nonblank]) if len(nonblank) > 1 else set()
    rows = []
    for k, f in enumerate(files):
        t = toks[f]
        if any(BLANK_TOKEN.match(x) for x in t):
            tag = next((x for x in reversed(t) if x.isdigit()), "")
            rows.append(dict(Filename=f, Group="Blank", Time="—", TimeValue=np.nan, TimeUnit="", SampleNo=tag,
                             SampleType=TYPE_BLANK, Sample="Blank" + (f"_{tag}" if tag else "")))
            continue
        rest = [x for x in t if x not in common]
        is_ref = any(REF_TOKEN.match(x) for x in rest)
        tm = next((TIME_TOKEN.match(x) for x in rest if TIME_TOKEN.match(x)), None)
        nums = [x for x in rest if x.isdigit()]
        group_tokens = [x for x in rest if not x.isdigit() and not REF_TOKEN.match(x) and not TIME_TOKEN.match(x)]
        group = "_".join(group_tokens) or "Sample"
        num = nums[-1] if nums else str(k + 1)
        if is_ref:
            rows.append(dict(Filename=f, Group=group, Time="Unlabeled", TimeValue=np.nan, TimeUnit="", SampleNo=num,
                             SampleType=TYPE_REF, Sample=f"12C_Ref_{num}"))
            continue
        if tm:
            val = float(tm.group(1)); unit = _unit_norm(tm.group(2))
            tlabel = f"{val:g} {unit}"
        else:
            val, unit, tlabel = np.nan, "", "—"
        rows.append(dict(Filename=f, Group=group, Time=tlabel, TimeValue=val, TimeUnit=unit, SampleNo=num,
                         SampleType=TYPE_LAB, Sample=f"{group}_{tlabel.replace(' ', '')}_{num}"))
    return pd.DataFrame(rows)


def order_samples(S: pd.DataFrame) -> pd.DataFrame:
    groups = list(dict.fromkeys(S.loc[S.SampleType == TYPE_LAB, "Group"]))
    groups.sort(key=lambda g: (0 if CONTROL_HINT.search(g) else 1, g))
    type_rank = {TYPE_REF: 0, TYPE_BLANK: 1, TYPE_LAB: 2}

    def key(r):
        tv = r.TimeValue if pd.notna(r.TimeValue) else 1e9
        mult = {"min": 1 / 60, "h": 1, "d": 24}.get(r.TimeUnit, 1)
        return (type_rank.get(r.SampleType, 3), groups.index(r.Group) if r.Group in groups else 99, tv * mult,
                str(r.SampleNo).zfill(8))
    S = S.iloc[sorted(range(len(S)), key=lambda i: key(S.iloc[i]))].reset_index(drop=True)
    if (S.SampleType == TYPE_BLANK).sum() == 1:
        S.loc[S.SampleType == TYPE_BLANK, "Sample"] = "Blank"
    S["Block"] = np.where(S.SampleType == TYPE_LAB, S.Group + " · " + S.Time, S.SampleType)
    return S


def build_sample_table(df: pd.DataFrame) -> pd.DataFrame:
    if reference_recognises(df.Filename):
        return build_sample_table_reference(df)
    files = df[["Filename", "Study File ID"]].drop_duplicates()
    S = parse_filenames(files.Filename)
    S["StudyFileID"] = S.Filename.map(dict(zip(files.Filename, files["Study File ID"])))
    if S.Sample.duplicated().any():
        dup = S.Sample[S.Sample.duplicated()].tolist()
        raise ValueError(f"Two raw files map to the same sample name {dup}; rename the files or edit the design.")
    return order_samples(S)


# ---------------------------------------------------------------------------
# Compound handling
# ---------------------------------------------------------------------------
def carbon_count(formula: str) -> int:
    """Number of carbons in a space-separated formula like 'C4 H9 N3 O2'."""
    m = re.search(r"(?:^|\s)C(\d*)(?=\s|$)", str(formula))
    return int(m.group(1) or 1) if m else MAX_ISO


def assign_peaks(df: pd.DataFrame) -> pd.DataFrame:
    """
    Give each row a peak number within its compound name, matching peaks across
    samples by retention time, so no row is ever merged or dropped.
    """
    df = df.copy()
    df["Peak"] = 0
    for _, g in df.groupby("Name", sort=False):
        n_peaks = g.groupby("Sample").size().max()
        # reference RTs = median RT of each peak position over samples that have all peaks
        full = [np.sort(s["RT in min"].values) for _, s in g.groupby("Sample") if len(s) == n_peaks]
        ref_rt = np.median(np.vstack(full), axis=0)
        for _, s in g.groupby("Sample"):
            cost = np.abs(s["RT in min"].values[:, None] - ref_rt[None, :])
            r, c = linear_sum_assignment(cost)
            df.loc[s.index[r], "Peak"] = c + 1
    assert not df.duplicated(["Name", "Peak", "Sample"]).any()
    df["nPeaks"] = df.groupby("Name")["Peak"].transform("max")
    df["PeakLabel"] = df["Peak"].astype(str) + " of " + df["nPeaks"].astype(str)
    df["RT_mean"] = df.groupby(["Name", "Peak"])["RT in min"].transform("mean").round(2)
    return df


def build_wide_table(df, S, mode, require_all=True, min_area: float = 1000.0):
    samples = list(S.Sample)
    lab = list(S.loc[S.SampleType == TYPE_LAB, "Sample"])
    fixed = ["Method", "Compound", "Formula", "MW", "# C", "Peak", "RT (min)", "Isotopologue"]
    rows = []
    for name in dict.fromkeys(df.Name):  # keep source compound order
        g = df[df.Name == name]
        for p in sorted(g.Peak.unique()):
            gp = g[g.Peak == p].set_index("Sample")
            first = gp.iloc[0]
            n = min(int(first.nC), MAX_ISO)
            for i in range(n + 1):
                vals = [gp.at[s, RA_COLS[i]] if s in gp.index else np.nan for s in samples]
                rows.append([mode, name, first.Formula, first["Molecular Weight"], int(first.nC),
                             first.PeakLabel, first.RT_mean, f"M+{i}", *vals])
    W = pd.DataFrame(rows, columns=fixed + samples)

    # a peak whose M+0 is 0 in any labeled sample is not a reliable peak -> drop the whole peak
    m0 = W[W.Isotopologue == "M+0"]
    bad = m0[(m0[lab] == 0).any(axis=1)][["Compound", "Peak"]]
    bad_peaks = set(map(tuple, bad.values))
    dropped_peaks = [f"{mode}: {c} (peak {p})" for c, p in sorted(bad_peaks)]
    W = W[[(c, p) not in bad_peaks for c, p in zip(W.Compound, W.Peak)]].reset_index(drop=True)

    is_m0 = W.Isotopologue == "M+0"
    if require_all:
        keep = ~(W[lab] == 0).any(axis=1) | is_m0           # zero in any labeled sample -> drop
        rule = ("An isotopologue row is kept only if it has a value in every labeled sample. "
                "If it is 0 in any one labeled sample, the row was removed")
    else:
        keep = (W[lab].fillna(0) != 0).any(axis=1) | is_m0  # zero in all labeled samples -> drop
        rule = "Isotopologue rows that are 0 in all labeled samples were removed"
    rule += (" (12C reference and blank are not considered). If M+0 itself is 0 in any labeled sample, "
             "the whole peak is removed (unreliable peak). "
             "A sample where the whole compound peak was not reported (empty) does not count as 0.")
    removed = int((~keep).sum())
    W = W[keep].reset_index(drop=True)

    # v2.3 low-value rule (after the zero rule of sil_flux_pipeline.py): a value below min_area in ANY labeled
    # sample is not reliable. M+0 below min_area -> the whole peak is removed; any other isotopologue below
    # min_area -> that isotopologue row is removed. Empty cells (peak not reported) do not count.
    low = []
    if min_area and min_area > 0:
        is_low = (W[lab] < min_area).any(axis=1)          # NaN < x is False
        m0 = W.Isotopologue == "M+0"
        low_peaks = set(zip(W.Compound[is_low & m0], W.Peak[is_low & m0]))
        low += [f"{mode}: {c} (peak {p}) — M+0 < {min_area:,g}" for c, p in sorted(low_peaks)]
        in_low_peak = np.array([(c, p) in low_peaks for c, p in zip(W.Compound, W.Peak)], dtype=bool)
        n_rows = int((is_low & ~m0 & ~in_low_peak).sum())
        W = W[~in_low_peak & ~(is_low & ~m0)].reset_index(drop=True)
        removed += n_rows
        rule += (f" Low-value rule: an isotopologue row with a value below {min_area:,g} in any labeled sample was "
                 f"removed; if M+0 is below {min_area:,g} in any labeled sample, the whole peak was removed "
                 "(12C reference and blank are not considered; empty cells do not count).")
        low.append(f"{mode}: low-value rule removed {len(low_peaks)} peak(s) and {n_rows} isotopologue row(s)")
    return W, fixed, rule, removed, dropped_peaks + low



# ---------------------------------------------------------------------------
# One peak per metabolite
# ---------------------------------------------------------------------------
def name_key(name: str, aliases: dict | None = None) -> str:
    """Name used for matching across methods: alias, lower case, no stereo prefix, sorted words
    (identical to sil_flux_pipeline.name_key)."""
    al = aliases if aliases is not None else ALIASES
    return _raw_key(al.get(name, name))


RT_TOL = 0.05  # min; peaks closer than this to a larger annotated compound count as co-eluting
# common in-source neutral losses: larger compound = candidate + one of these -> candidate may be a fragment
NEUTRAL_LOSSES = {
    "H2O": "water", "NH3": "ammonia", "CO2": "CO2", "CH2O2": "formic acid", "H3PO4": "phosphoric acid",
    "HPO3": "metaphosphate", "C5H8O4": "ribose", "C5H8O3": "deoxyribose", "C6H10O5": "hexose",
    "C5H9O7P": "ribose-phosphate", "C6H11O8P": "hexose-phosphate", "H4P2O7": "pyrophosphate",
    "C2H5NO2": "glycine", "C3H5NO2": "dehydroalanine", "C5H7NO3": "pyroglutamate", "C3H9N": "trimethylamine",
    "H2O2": "two waters", "H5NO": "water + ammonia",
}


def formula_dict(f):
    d = {}
    for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", str(f).replace(" ", "")):
        d[el] = d.get(el, 0) + int(n or 1)
    return d


def neutral_loss(big, small):
    """Name of the neutral loss if formula(big) - formula(small) is a known loss, else ''."""
    a, b = formula_dict(big), formula_dict(small)
    diff = {el: a.get(el, 0) - b.get(el, 0) for el in set(a) | set(b)}
    if any(v < 0 for v in diff.values()):
        return ""
    for f, name in NEUTRAL_LOSSES.items():
        fd = formula_dict(f)
        if all(diff.get(el, 0) == fd.get(el, 0) for el in set(diff) | set(fd)):
            return name
    return ""


def candidate_peaks(df, W, S):
    """One row per (method, compound, peak) still in the table, with the metrics used to choose."""
    lab = set(S.loc[S.SampleType == TYPE_LAB, "Sample"])
    n_lab = len(lab)
    in_table = set(map(tuple, W[["Method", "Compound", "Peak"]].drop_duplicates().values))
    iso = W[W.Isotopologue != "M+0"].groupby(["Method", "Compound", "Peak"]).size()
    rows = []
    for (mode, name, pk), gp in df[df.Sample.isin(lab)].groupby(["Mode", "Name", "PeakLabel"], sort=False):
        if (mode, name, pk) not in in_table:
            continue
        rows.append(dict(Method=mode, Compound=name, Peak=pk, Formula=gp.Formula.iloc[0],
                         MW=float(gp["Molecular Weight"].iloc[0]), RT=round(float(gp["RT in min"].mean()), 2),
                         Detected=gp.Sample.nunique(), nLab=n_lab, LabeledIso=int(iso.get((mode, name, pk), 0)),
                         MedianArea=float(gp.Area.median()), Warnings=int((gp.Status != "NoWarnings").sum())))
    E = pd.DataFrame(rows)
    E["Key"] = E.Formula.str.replace(" ", "") + "|" + E.Compound.map(name_key)

    # co-elution with a LARGER annotated compound in the same method -> possible in-source fragment
    best_area = E.groupby(["Method", "Key"]).MedianArea.transform("max")
    major = E[E.MedianArea >= 0.5 * best_area]          # main peaks of each compound
    notes = []
    for _, r in E.iterrows():
        m = major[(major.Method == r.Method) & (major.Key != r.Key) & (major.MW > r.MW + 1.0) &
                  ((major.RT - r.RT).abs() <= RT_TOL)]
        hits = []
        for _, x in m.iterrows():
            loss = neutral_loss(x.Formula, r.Formula)
            if loss:
                hits.append(f"{x.Compound} (RT {x.RT}; −{loss})")
        notes.append("; ".join(hits))
    E["CoElutesWith"] = notes
    return E


def compare(a, b):
    """Return (winner, loser, reason) between two candidate peaks."""
    n = a.nLab
    if bool(a.CoElutesWith) != bool(b.CoElutesWith):
        w, l = (b, a) if a.CoElutesWith else (a, b)
        return w, l, f"other peak co-elutes with {l.CoElutesWith.split(' (')[0]} (possible in-source fragment)"
    if a.Detected != b.Detected:
        w, l = (a, b) if a.Detected > b.Detected else (b, a)
        return w, l, f"detected in more samples ({w.Detected}/{n} vs {l.Detected}/{n})"
    if a.LabeledIso != b.LabeledIso:
        w, l = (a, b) if a.LabeledIso > b.LabeledIso else (b, a)
        return w, l, f"more labeled isotopologues ({w.LabeledIso} vs {l.LabeledIso})"
    hi, lo = (a, b) if a.MedianArea >= b.MedianArea else (b, a)
    fold = hi.MedianArea / lo.MedianArea if lo.MedianArea > 0 else np.inf
    if fold >= SIGNAL_FOLD:
        return hi, lo, f"higher signal ({fold:.1f}x median area)"
    if a.Warnings != b.Warnings:
        w, l = (a, b) if a.Warnings < b.Warnings else (b, a)
        return w, l, f"fewer peak warnings ({w.Warnings} vs {l.Warnings}); signal similar ({fold:.1f}x)"
    return hi, lo, f"signal slightly higher ({fold:.2f}x); all other checks equal"


def desc(r):
    return f"{r.Method} {r.Compound} peak {r.Peak} @ {r.RT} min"


def resolve_unique(W, df, S, prefer: dict | None = None):
    PREFER_ = PREFER if prefer is None else prefer
    """
    Keep exactly one peak per metabolite. Candidates = every peak of every name that maps to the same
    metabolite (same formula + same normalized name / alias), in every method.
    Returns (W_kept, W_removed, decisions table, kept set of (method, compound, peak)).
    """
    E = candidate_peaks(df, W, S)
    kept, note, rows = set(), {}, []
    for key, g in E.groupby("Key", sort=False):
        cands = [r for _, r in g.iterrows()]
        forced = [PREFER_[c] for c in g.Compound if c in PREFER_]
        if forced:
            pool = [r for r in cands if r.Method == forced[0]] or cands
        else:
            pool = cands
        win = pool[0]
        for c in pool[1:]:
            win, _, _ = compare(win, c)
        kept.add((win.Method, win.Compound, win.Peak))
        if len(cands) == 1:
            continue
        others = [r for r in cands if r is not win and not (r.Method == win.Method and r.Compound == win.Compound and r.Peak == win.Peak)]
        reasons = {}
        for r in others:
            if forced and r.Method != forced[0]:
                reasons[id(r)] = "manual choice (PREFER)"
            else:
                reasons[id(r)] = compare(win, r)[2]
        runner = others[0] if not others else max(others, key=lambda r: (not r.CoElutesWith, r.Detected, r.LabeledIso, r.MedianArea))
        what = []
        if g.Method.nunique() > 1: what.append(f"{g.Method.nunique()} methods")
        if g.Compound.nunique() > 1: what.append(f"{g.Compound.nunique()} names")
        if len(cands) > g.Compound.nunique() * 1 and (g.groupby(["Method", "Compound"]).size() > 1).any(): what.append("several peaks")
        note[(win.Method, win.Compound, win.Peak)] = f"Kept: best of {len(cands)} candidates ({', '.join(what)}) — {reasons[id(runner)]}"
        for r in others:
            note[(r.Method, r.Compound, r.Peak)] = f"Removed: {desc(win)} kept ({reasons[id(r)]})"
        for r in [win] + others:
            rows.append(dict(Metabolite=win.Compound, Kept="Yes" if r is win else "No", Method=r.Method, Name=r.Compound,
                             Peak=r.Peak, RT=r.RT, Detected=f"{r.Detected}/{r.nLab}", LabeledIso=r.LabeledIso,
                             MedianArea=r.MedianArea, Warnings=r.Warnings, CoElutesWith=r.CoElutesWith,
                             Reason=("chosen" if r is win else reasons[id(r)])))
    W = W.copy()
    trip = list(zip(W.Method, W.Compound, W.Peak))
    W.insert(1, "Selection", [note.get(t, "") for t in trip])
    keep = np.array([t in kept for t in trip], dtype=bool)
    return W[keep].reset_index(drop=True), W[~keep].reset_index(drop=True), pd.DataFrame(rows), kept


def exclude_m0_only(W: pd.DataFrame, W_removed: pd.DataFrame, kept: set):
    """
    v2.2 criterion. After one peak per metabolite is chosen, a metabolite whose kept peak has only the M+0
    isotopologue row (every labeled isotopologue M+1…M+n was removed by the zero rule) carries no labeling
    information and is excluded from the combined table. Its rows move to Removed_peaks.
    Returns (W, W_removed, kept, excluded compound names).
    """
    key = list(zip(W.Method, W.Compound, W.Peak))
    n_lab_iso = W[W.Isotopologue != "M+0"].groupby(["Method", "Compound", "Peak"]).size()
    only = {k for k in dict.fromkeys(key) if int(n_lab_iso.get(k, 0)) == 0}
    if not only:
        return W, W_removed, kept, []
    mask = np.array([k in only for k in key], dtype=bool)
    moved = W[mask].copy()
    moved["Selection"] = "Removed: only M+0 isotopologue kept (no labeled isotopologue) — excluded from combined table"
    W_removed = pd.concat([W_removed, moved], ignore_index=True)
    names = list(dict.fromkeys(moved.Compound))
    return W[~mask].reset_index(drop=True), W_removed, kept - only, names


# ---------------------------------------------------------------------------
# Excel writing
# ---------------------------------------------------------------------------
FONT = "Arial"
N = Font(name=FONT, size=10)
B = Font(name=FONT, size=10, bold=True)
HW = Font(name=FONT, size=10, bold=True, color="FFFFFF")
NAVY = "1F3864"
HEAD_FILL = {TYPE_REF: "7F7F7F", TYPE_BLANK: "A6A6A6", "shControl": "2F75B5", "GeneKD": "C65911", "Other": "595959"}
LIGHT_FILL = {TYPE_REF: "EDEDED", TYPE_BLANK: "F2F2F2", "shControl": "DDEBF7", "GeneKD": "FCE4D6", "Other": "F2F2F2"}
STATUS_FILL = {"Gapped": "FFEB9C", "Contaminated": "FFC7CE", "LowFit": "FFC7CE"}
BLOCK_EDGE = Side(style="medium", color="595959")
ROW_EDGE = Side(style="thin", color="808080")


def fill(hex_):
    return PatternFill("solid", fgColor=hex_)


def family(row):
    f = row.SampleType if row.SampleType != TYPE_LAB else row.Group
    return f if f in HEAD_FILL else "Other"


def clean(v):
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return None
    return v.item() if isinstance(v, np.generic) else v


def header_row(ws, columns, row=1):
    for j, c in enumerate(columns, 1):
        x = ws.cell(row, j, c)
        x.font, x.fill = HW, fill(NAVY)
        x.alignment = Alignment(wrap_text=True, vertical="center")


def write_readme(wb, sources, df, W, S, rule, W_removed=None, D=None, dropped_peaks=(), m0_only=()):
    ws = wb.active
    ws.title = "README"
    n_lab = int((S.SampleType == TYPE_LAB).sum())
    lines = [
        ("13C tracing — reorganized isotopologue data (no calculations)", True),
        *[(f"Source ({m}): {f}", False) for m, f in sources],
        ("All numbers are copied unchanged from the source files.", False),
        ("", False),
        ("SHEETS", True),
        ('Isotopologue_Areas — one row per compound × isotopologue, one column per sample. '
         'Value = "Relative Area i" from the source for M+i. The Method column says which acquisition '
         '(e.g. RP_Pos, HILIC_Neg) the row comes from; methods are stacked one after the other. '
         'Each metabolite appears ONCE, as a single peak from a single method (see Selection column).', False),
        ("Selection_Decisions — every metabolite that had more than one candidate peak (several peaks, several "
         "names for the same compound, or both methods): metrics of every candidate, which one was kept and why.", False),
        ("Removed_peaks — the isotopologue rows of all candidates that were not kept (same layout).", False),
        ("Samples — what each raw file is: group, time point, sample number, type.", False),
        ("Peak_Details — one row per compound peak per sample: RT, total area and Compound Discoverer "
         "status (NoWarnings / Gapped / Contaminated / LowFit).", False),
        ("", False),
        ("HOW THE SOURCE WAS RESHAPED", True),
        ('Source columns "Relative Area 0 … 25" became rows M+0 … M+25, limited to M+n '
         "(n = number of carbons in the formula).", False),
        ("Sample columns are grouped: 12C reference, Blank, then each group × time point "
         "(replicates labeled by sample number).", False),
        ('Compound names reported for more than one chromatographic peak in the same sample are kept as '
         'separate blocks ("Peak 1 of 2", "Peak 2 of 2") with their mean RT. Peaks were matched across '
         "samples by retention time.", False),
        (rule, False),
        *[(f"Peak removed (M+0 = 0 or below the minimum value in a labeled sample) — {x}", False) for x in dropped_peaks],
        ("Empty cell = the compound peak was not reported for that sample. 0 = reported, area is 0.", False),
        ("", False),
        ("ONE PEAK PER METABOLITE", True),
        ("Candidates: every peak of every name that is the same metabolite — same formula and same name after "
         "ignoring L-/D-/DL- prefixes, punctuation and word order, plus known synonyms (e.g. Inosinic acid = IMP) — "
         "in every method. Same-formula compounds with genuinely different names (possible isomers) are NOT merged.", False),
        ("Candidates are ranked in this order: (1) does NOT co-elute (within "
         f"{RT_TOL} min) with a larger annotated compound in the same method — such a peak is likely an in-source "
         "fragment, e.g. hypoxanthine at the inosine RT; (2) detected in more labeled samples; (3) more labeled "
         f"isotopologues kept (M+1 and up); (4) higher signal, median total area, if at least {SIGNAL_FOLD}x different; "
         "(5) fewer Compound Discoverer warnings; (6) higher signal.", False),
        ("Rule (1) only applies when the larger co-eluting compound's formula equals this compound plus a common "
         "in-source neutral loss (water, ammonia, CO2, formic acid, phosphate, ribose, deoxyribose, hexose, "
         "ribose-/hexose-phosphate, pyrophosphate, glycine, pyroglutamate, trimethylamine).", False),
        ("All other candidates move to Removed_peaks. To force a method, use 'Preferred method' in Step 1 (Advanced).", False),
        ("", False),
        ("METABOLITES WITH ONLY M+0", True),
        ("After one peak per metabolite is chosen, a metabolite whose kept peak has only the M+0 isotopologue row "
         "(no labeled isotopologue survived the rules above) is excluded from Isotopologue_Areas and listed in "
         f"Removed_peaks. Excluded: {len(m0_only)}" + (" — " + ", ".join(m0_only) if len(m0_only) else "") + ".", False),
        ("", False),
        ("COUNTS", True),
        *[(f"{m}: {d.Name.nunique()} compound names in source · "
           f"{W[W.Method == m].Compound.nunique()} kept in final table "
           f"({W[W.Method == m][['Compound', 'Peak']].drop_duplicates().shape[0]} peaks, "
           f"{(W.Method == m).sum()} isotopologue rows) · {len(d)} source rows.", False)
          for m, d in df.groupby("Mode", sort=False)],
        (f"Metabolites that had several candidates: {0 if D is None or D.empty else D.Metabolite.nunique()}; "
         f"{0 if D is None or D.empty else int((D.Kept == 'No').sum())} candidate peaks moved to Removed_peaks.", False),
        (f"{len(S)} samples ({n_lab} labeled).", False),
    ]
    for i, (text, bold) in enumerate(lines, 1):
        c = ws.cell(i, 1, text)
        c.font = Font(name=FONT, size=13 if i == 1 else 10.5, bold=bold, color=NAVY if bold else "000000")
        c.alignment = Alignment(wrap_text=True, vertical="top")
    ws.column_dimensions["A"].width = 130


def write_isotopologues(wb, W, fixed, S, title="Isotopologue_Areas"):
    ws = wb.create_sheet(title)
    nf = len(fixed)
    samples = list(S.Sample)
    srow = S.set_index("Sample")

    # header row 1 = fixed column names (merged over 2 rows) + sample blocks; row 2 = sample names
    for j, c in enumerate(fixed, 1):
        ws.merge_cells(start_row=1, end_row=2, start_column=j, end_column=j)
        x = ws.cell(1, j, c)
        x.font, x.fill = HW, fill(NAVY)
        x.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
    block_start, j = set(), nf + 1
    for block, grp in S.groupby("Block", sort=False):
        fam = family(grp.iloc[0])
        block_start.add(j)
        if len(grp) > 1:
            ws.merge_cells(start_row=1, end_row=1, start_column=j, end_column=j + len(grp) - 1)
        x = ws.cell(1, j, block)
        x.font, x.fill, x.alignment = HW, fill(HEAD_FILL[fam]), Alignment(horizontal="center")
        for k, s in enumerate(grp.Sample):
            y = ws.cell(2, j + k, s)
            y.font, y.fill = B, fill(LIGHT_FILL[fam])
            y.alignment = Alignment(horizontal="center", wrap_text=True)
            y.border = Border(left=BLOCK_EDGE if k == 0 else None)
            ws.column_dimensions[get_column_letter(j + k)].width = 15
        j += len(grp)

    # data rows, banded per compound peak
    band, prev, prev_mode = False, None, None
    for i, row in enumerate(W.itertuples(index=False), 3):
        if row[0] != prev_mode and prev_mode is not None:
            for jj in range(1, nf + len(samples) + 1):
                ws.cell(i - 1, jj).border = Border(bottom=Side(style="thick", color="1F3864"))
        prev_mode = row[0]
        key = (row[0], row[fixed.index("Compound")], row[fixed.index("Peak")])
        top = key != prev
        if top:
            band, prev = not band, key
        for j, v in enumerate(row, 1):
            c = ws.cell(i, j, clean(v))
            c.font = B if j in (fixed.index("Compound") + 1, nf) else N
            if j > nf:
                c.number_format = "#,##0"
            elif j == fixed.index("MW") + 1:
                c.number_format = "0.0000"
            if band and j <= nf:
                c.fill = fill("F2F2F2")
            c.border = Border(top=ROW_EDGE if top else None, left=BLOCK_EDGE if j in block_start else None)

    widths = {"Method": 11, "Selection": 34, "Compound": 34, "Formula": 15, "MW": 10, "# C": 5,
              "Peak": 8, "RT (min)": 8, "Isotopologue": 11}
    for j, c in enumerate(fixed, 1):
        ws.column_dimensions[get_column_letter(j)].width = widths.get(c, 12)
    if "Selection" in fixed:
        jd = fixed.index("Selection") + 1
        for i in range(3, len(W) + 3):
            v = ws.cell(i, jd).value or ""
            if v.startswith("Kept"):
                ws.cell(i, jd).fill = fill("E2EFDA")
            elif v.startswith("Removed"):
                ws.cell(i, jd).fill = fill("FCE4D6")
    ws.row_dimensions[2].height = 30
    ws.freeze_panes = ws.cell(3, nf + 1)
    ws.auto_filter.ref = f"A2:{get_column_letter(nf + len(samples))}{len(W) + 2}"


def write_decisions(wb, D):
    ws = wb.create_sheet("Selection_Decisions")
    if D.empty:
        ws.cell(1, 1, "Every metabolite had a single candidate peak.").font = B
        return
    cols = ["Metabolite", "Kept", "Method", "Name", "Peak", "RT", "Detected", "LabeledIso", "MedianArea",
            "Warnings", "CoElutesWith", "Reason"]
    header_row(ws, ["Metabolite (kept name)", "Kept", "Method", "Name in file", "Peak", "RT (min)", "Detected",
                    "Labeled isotopologues", "Median total area", "Warnings", "Possible in-source fragment of", "Reason"])
    prev = None
    for i, r in enumerate(D[cols].itertuples(index=False), 2):
        top = r[0] != prev
        prev = r[0]
        for j, v in enumerate(r, 1):
            c = ws.cell(i, j, clean(v))
            c.font = B if (j == 1 or (j == 2 and v == "Yes")) else N
            if j == 9:
                c.number_format = "#,##0"
            if top:
                c.border = Border(top=ROW_EDGE)
        ws.cell(i, 2).fill = fill("E2EFDA" if r[1] == "Yes" else "FCE4D6")
        if r[10]:
            ws.cell(i, 11).fill = fill("FFEB9C")
    for j, w in enumerate([30, 7, 11, 30, 8, 8, 9, 11, 14, 9, 36, 60], 1):
        ws.column_dimensions[get_column_letter(j)].width = w
    ws.row_dimensions[1].height = 32
    ws.freeze_panes = "B2"
    ws.auto_filter.ref = f"A1:{get_column_letter(len(cols))}{len(D) + 1}"


def write_samples(wb, S, modes):
    ws = wb.create_sheet("Samples")
    file_cols = [f"File {m}" for m in modes]
    out = S[["Sample", "SampleType", "Group", "Time", "SampleNo"] + file_cols]
    header_row(ws, ["Sample", "Sample Type", "Group", "Time", "Sample No."] + [f"Raw file ({m})" for m in modes])
    for i, r in enumerate(S.itertuples(index=False), 2):
        f = fill(LIGHT_FILL[family(r)])
        for j, col in enumerate(out.columns, 1):
            c = ws.cell(i, j, clean(r[S.columns.get_loc(col)]))
            c.font, c.fill = N, f
    for j, w in enumerate([20, 15, 12, 16, 11] + [52] * len(modes), 1):
        ws.column_dimensions[get_column_letter(j)].width = w
    ws.freeze_panes = "B2"


def write_peak_details(wb, df, S, kept=frozenset()):
    ws = wb.create_sheet("Peak_Details")
    sidx = {s: i for i, s in enumerate(S.Sample)}
    cidx = {n: i for i, n in enumerate(dict.fromkeys(df.Name))}
    P = df.merge(S[["Sample", "SampleType", "Group", "Time"]], on="Sample")
    midx = {m: i for i, m in enumerate(dict.fromkeys(df.Mode))}
    P = P.assign(_m=P.Mode.map(midx), _c=P.Name.map(cidx), _s=P.Sample.map(sidx)).sort_values(["_m", "_c", "Peak", "_s"])
    P["InFinal"] = ["Yes" if (m, n, p) in kept else "No" for m, n, p in zip(P.Mode, P.Name, P.PeakLabel)]
    cols = ["Mode", "InFinal", "Name", "Formula", "PeakLabel", "Sample", "SampleType", "Group", "Time",
            "RT in min", "Area", "Status", "Study File ID", "Filename"]
    header_row(ws, ["Method", "In final table", "Compound", "Formula", "Peak", "Sample", "Sample Type", "Group", "Time", "RT (min)",
                    "Total Area (sum of M+0…M+25)", "Status", "Study File ID", "Filename"])
    for i, r in enumerate(P[cols].itertuples(index=False), 2):
        for j, v in enumerate(r, 1):
            c = ws.cell(i, j, clean(v))
            c.font = N
            if j == 11:
                c.number_format = "#,##0"
            if j == 12 and v in STATUS_FILL:
                c.fill = fill(STATUS_FILL[v])
    for j, w in enumerate([11, 14, 34, 15, 8, 18, 14, 11, 15, 9, 18, 13, 9, 50], 1):
        ws.column_dimensions[get_column_letter(j)].width = w
    ws.row_dimensions[1].height = 30
    ws.freeze_panes = "B2"
    ws.auto_filter.ref = f"A1:{get_column_letter(len(cols))}{len(P) + 1}"


# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------
@dataclass
class ReorgResult:
    W: pd.DataFrame                 # final wide table (one peak per metabolite), display-filtered isotopologues
    W_removed: pd.DataFrame         # candidates not kept
    decisions: pd.DataFrame         # Selection_Decisions
    df: pd.DataFrame                # all source rows (long), with Mode / Sample / Peak / nC
    S: pd.DataFrame                 # sample table
    kept: set                       # {(method, compound, peak label)}
    fixed: list
    rule: str
    dropped_peaks: list
    sources: list                   # [(method, file name)]
    log: list = field(default_factory=list)
    m0_only: list = field(default_factory=list)   # metabolites excluded: only M+0 left (v2.2)
    naming: str = "reference"                      # "reference" (FILE_PATTERN) or "generic"

    def kept_rows(self) -> pd.DataFrame:
        """Source rows (all isotopologues M+0..M+25) of the kept peaks only — used for fractions."""
        k = self.kept
        m = [(a, b, c) in k for a, b, c in zip(self.df.Mode, self.df.Name, self.df.PeakLabel)]
        return self.df[m].copy()

    def workbook_bytes(self) -> bytes:
        if getattr(self, "_wb_cache", None):
            return self._wb_cache
        wb = Workbook()
        modes = [m for m, _ in self.sources]
        write_readme(wb, self.sources, self.df, self.W, self.S, self.rule, self.W_removed, self.decisions,
                     self.dropped_peaks, self.m0_only)
        write_isotopologues(wb, self.W, self.fixed, self.S)
        write_decisions(wb, self.decisions)
        write_isotopologues(wb, self.W_removed, self.fixed, self.S, title="Removed_peaks")
        write_samples(wb, self.S, modes)
        write_peak_details(wb, self.df, self.S, self.kept)
        buf = io.BytesIO()
        wb.save(buf)
        self._wb_cache = buf.getvalue()
        return self._wb_cache


def read_table(src, name: str | None = None) -> pd.DataFrame:
    """Read an Excel or CSV export from a path or a file-like object."""
    name = name or str(getattr(src, "name", src))
    if name.lower().endswith((".csv", ".txt", ".tsv")):
        sep = "\t" if name.lower().endswith((".tsv", ".txt")) else ","
        return pd.read_csv(src, sep=sep)
    return pd.read_excel(src)


def process_files(inputs, keep_partial: bool = False, aliases: dict | None = None,
                  prefer: dict | None = None, min_area: float = 1000.0) -> ReorgResult:
    """
    inputs: list of (file name, path-or-file-like) pairs, or a list of paths.
    """
    if aliases is not None:
        ALIASES.clear(); ALIASES.update(aliases)
    pairs = [(Path(x).name, x) if not isinstance(x, tuple) else x for x in inputs]
    per_mode, log = [], []
    for name, src in pairs:
        d = read_table(src, name)
        missing = [c for c in REQUIRED_COLS if c not in d.columns]
        if missing:
            raise ValueError(f"{name}: missing expected columns {missing[:6]}{'…' if len(missing) > 6 else ''}. "
                             "Export the Compound Discoverer SIL table with 'Relative Area 0…25' columns.")
        mode = mode_of(name)
        if mode in [m for m, *_ in per_mode]:
            mode = f"{mode}_{len(per_mode) + 1}"
        S_i = build_sample_table(d)
        d["Sample"] = d.Filename.map(dict(zip(S_i.Filename, S_i.Sample)))
        d["Mode"] = mode
        d["nC"] = d.Formula.map(carbon_count)
        per_mode.append((mode, name, d, S_i))
        log.append(f"{mode}: {name} — {len(d)} rows, {d.Name.nunique()} compound names, {len(S_i)} raw files")

    modes = [m for m, *_ in per_mode]
    S = pd.concat([S_i.assign(Mode=m) for m, _, _, S_i in per_mode])
    files = S.pivot(index="Sample", columns="Mode", values="Filename")
    S = S.drop_duplicates("Sample").drop(columns=["Mode", "Filename", "StudyFileID"])
    naming = "reference" if all(reference_recognises(d.Filename) for _, _, d, _ in per_mode) else "generic"
    if naming == "reference":   # exactly as sil_flux_pipeline.run()
        S = S.iloc[sorted(range(len(S)), key=lambda i: _ref_sort_key(S.iloc[i]))].reset_index(drop=True)
    else:
        S = order_samples(S.drop(columns=["Block"]))
        log.append("Note: raw-file names do not match FILE_PATTERN; groups and times were read by the generic "
                   "parser (check them in Step 2).")
    for m in modes:
        S[f"File {m}"] = S.Sample.map(files[m]) if m in files else np.nan
        absent = S.Sample[S[f"File {m}"].isna()].tolist()
        if absent:
            log.append(f"Note: {m} has no file for: {', '.join(absent)}")

    tables, dfs, dropped = [], [], []
    rule = ""
    for mode, _, d, _ in per_mode:
        d = assign_peaks(d)
        W_i, fixed, rule, rem, dp = build_wide_table(d, S, mode, require_all=not keep_partial, min_area=min_area)
        tables.append(W_i); dfs.append(d); dropped += [x for x in dp if "low-value rule removed" not in x]
        log += [x for x in dp if "low-value rule removed" in x]
        log.append(f"{mode}: {W_i[['Compound', 'Peak']].drop_duplicates().shape[0]} compound peaks, "
                   f"{len(W_i)} isotopologue rows kept ({rem} zero / low-value rows removed)")
    W = pd.concat(tables, ignore_index=True)
    df = pd.concat(dfs, ignore_index=True)
    W, W_removed, D, kept = resolve_unique(W, df, S, prefer)
    fixed = ["Method", "Selection"] + fixed[1:]
    # v2.2 criterion: metabolites whose kept peak has ONLY the M+0 isotopologue are excluded from the combined table
    W, W_removed, kept, m0_only = exclude_m0_only(W, W_removed, kept)
    if len(D):
        log.append(f"One peak per metabolite: {D.Metabolite.nunique()} metabolites had several candidates; "
                   f"{int((D.Kept == 'No').sum())} candidate peaks removed "
                   f"({int((D.CoElutesWith != '').sum())} flagged as possible in-source fragments)")
    if m0_only:
        log.append(f"{len(m0_only)} metabolite(s) excluded because only the M+0 isotopologue remained "
                   f"(no labeled isotopologue kept): " + ", ".join(m0_only[:40]) + ("…" if len(m0_only) > 40 else ""))
    log.append(f"Final: {W.Compound.nunique()} metabolites, {len(W)} isotopologue rows, {len(S)} samples "
               f"({int((S.SampleType == TYPE_LAB).sum())} labeled)")
    return ReorgResult(W=W, W_removed=W_removed, decisions=D, df=df, S=S, kept=kept, fixed=fixed, rule=rule,
                       dropped_peaks=dropped, sources=[(m, n) for m, n, *_ in per_mode], log=log,
                       m0_only=m0_only, naming=naming)
