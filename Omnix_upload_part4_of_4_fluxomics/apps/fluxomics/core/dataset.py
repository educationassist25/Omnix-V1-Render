"""
core.dataset — the single internal data model of the application.

Every analysis works on the CHANNEL TABLE: one row per metabolite × isotopologue (M+0 … M+n), one column per
sample, values = raw peak areas (NaN = not reported, 0 = reported as zero). Isotopologues are never summed for
intensity analysis.

Two input routes produce the same Dataset:
  * Compound Discoverer SIL exports  -> core.io_reorganize (one peak per metabolite) -> from_reorg()
  * Generic isotopologue table (+ optional sample-metadata file) -> read_generic()
      wide:  Metabolite | [Formula] | [Method] | Isotopologue | Sample1 | Sample2 | …
      long:  Metabolite | [Formula] | Isotopologue | Sample | Intensity
      metadata: Sample | Group | Time | [TimeUnit] | [Treatment] | [Replicate] | [Type] | [Include]
"""
from __future__ import annotations

import io
import re
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .io_reorganize import RA_COLS, MAX_ISO, TYPE_LAB, TYPE_REF, TYPE_BLANK, parse_filenames, order_samples
from .correction import parse_formula, TRACER_ELEMENT

ANNOT_COLS = {"kegg": "KEGG", "hmdb": "HMDB", "pubchem": "PubChemCID", "pubchemcid": "PubChemCID", "chebi": "ChEBI",
              "bigg": "BiGG"}
NON_SAMPLE = {"metabolite", "compound", "name", "feature", "isotopologue", "isotope", "label", "formula", "method",
              "mode", "rt", "rt (min)", "rt in min", "mz", "m/z", "mw", "molecular weight", "n", "labelatoms",
              "label atoms", "# c", "kegg", "hmdb", "pubchem", "pubchemcid", "chebi", "bigg", "pathway", "pathways",
              "class", "note", "notes", "peak", "selection", "duplicate decision", "sample", "intensity", "area",
              "value", "abundance", "channel", "i"}


@dataclass
class Dataset:
    CT: pd.DataFrame                    # Compound, Formula, Method, n, i, Channel, <samples…>
    S: pd.DataFrame                     # sample table (default design)
    source: str                         # "cd" or "generic"
    reorg: object = None                # ReorgResult for the CD route
    log: list = field(default_factory=list)
    annotations: pd.DataFrame | None = None   # IDs supplied in a generic table
    LT: pd.DataFrame | None = None      # v2.5.3 labeling table: ALL isotopologues M+0…M+n of each analysed metabolite
                                        # (incl. rows the zero / minimum-value rule removed from CT) — the
                                        # denominator of every isotopologue fraction (= CD "Exchange Rate in Percent")

    @property
    def labeling_table(self) -> pd.DataFrame:
        """Table used for isotopologue fractions, incorporation, mean label number and fractional enrichment."""
        return self.LT if self.LT is not None else self.CT

    @property
    def samples(self) -> list:
        return list(self.S.Sample)

    def channel_matrix(self, samples=None) -> pd.DataFrame:
        s = samples or self.samples
        X = self.CT.set_index("Channel")[[c for c in s if c in self.CT.columns]]
        return X.astype(float)

    def metabolites(self) -> pd.DataFrame:
        return self.CT[["Compound", "Formula", "Method"]].drop_duplicates("Compound").reset_index(drop=True)


def label_atoms(formula: str, tracer: str = "13C") -> int:
    el = TRACER_ELEMENT.get(tracer, "C")
    return int(parse_formula(formula).get(el, 0)) if isinstance(formula, str) and formula.strip() else 0


# ---------------------------------------------------------------------------
# Compound Discoverer route
# ---------------------------------------------------------------------------
def apply_zero_rule(CT: pd.DataFrame, labeled: list, min_area: float = 0.0) -> tuple[pd.DataFrame, list]:
    """
    Zero values are not valid data in labeled samples (rule of reorganize_sil_isotopologues.py):
      * if M+0 is 0 in any labeled sample, the compound is deleted with all its isotopologues;
      * an isotopologue row that is 0 in any labeled sample is deleted (M+0 is kept).
    Empty cells (compound not reported in a sample) do not count as 0.
    12C reference and blank samples are not considered.
    """
    lab = [c for c in labeled if c in CT.columns]
    log = []
    if not lab or CT.empty:
        return CT, log
    zero = (CT[lab] == 0).any(axis=1)
    if min_area and min_area > 0:                       # v2.3: values below min_area treated like zero
        zero = zero | (CT[lab] < min_area).any(axis=1)
    bad_m0 = set(CT.loc[zero & (CT.i == 0), "Compound"])
    out = CT[~CT.Compound.isin(bad_m0) & ~(zero & (CT.i > 0))].reset_index(drop=True)
    n_rows = int((zero & (CT.i > 0) & ~CT.Compound.isin(bad_m0)).sum())
    if bad_m0:
        log.append(f"{len(bad_m0)} compound(s) deleted with all isotopologues because M+0 = 0 (or below the minimum value) in a labeled sample: "
                   + ", ".join(sorted(bad_m0)[:30]) + ("…" if len(bad_m0) > 30 else ""))
    if n_rows:
        log.append(f"{n_rows} isotopologue row(s) deleted because the value is 0 or below the minimum value in a labeled sample "
                   "(zero is not valid data in labeled samples).")
    return out, log


def drop_m0_only(CT: pd.DataFrame) -> tuple[pd.DataFrame, list]:
    """Metabolites with only the M+0 isotopologue row carry no label information: excluded (v2.2)."""
    has_lab = CT[CT.i > 0].groupby("Compound").size()
    only = [c for c in dict.fromkeys(CT.Compound) if int(has_lab.get(c, 0)) == 0]
    if not only:
        return CT, []
    return (CT[~CT.Compound.isin(only)].reset_index(drop=True),
            [f"{len(only)} metabolite(s) excluded because only the M+0 isotopologue remained: "
             + ", ".join(only[:30]) + ("…" if len(only) > 30 else "")])


def from_reorg(r, tracer: str = "13C") -> Dataset:
    """
    Isotopologue rows of the reorganized (combined) table — one peak per metabolite, exactly as in the workbook of
    sil_flux_pipeline.py stage 1, minus metabolites with only M+0. Every isotopologue stays its own channel.
    """
    W = r.W
    samples = list(r.S.Sample)
    i = W.Isotopologue.str.replace("M+", "", regex=False).astype(int)
    n = [min(label_atoms(f, tracer) or int(c), MAX_ISO) for f, c in zip(W.Formula, W["# C"])]
    CT = pd.DataFrame({"Compound": W.Compound.values, "Formula": W.Formula.values, "Method": W.Method.values,
                       "n": n, "i": i.values})
    CT["Channel"] = CT.Compound + " | M+" + CT.i.astype(str)
    for smp in samples:
        CT[smp] = pd.to_numeric(W[smp], errors="coerce").values
    # v2.2: no second zero rule here — the table is exactly the reorganized (combined) table of Step 1,
    # built with the rules of sil_flux_pipeline.py (+ the M+0-only exclusion).
    S = r.S.copy()
    S["Treatment"] = ""
    LT = labeling_table_cd(r, CT, samples)
    log = list(r.log) + [f"Analysis table = combined table: {CT.Compound.nunique()} metabolites, "
                         f"{len(CT)} isotopologue rows.",
                         f"Labeling table: all {len(LT)} isotopologues M+0…M+n of these metabolites (Relative Area "
                         "0…n of the kept peak) — denominator of the isotopologue fractions, as in Compound "
                         "Discoverer's “Exchange Rate in Percent”."]
    return Dataset(CT=CT, S=S, source="cd", reorg=r, log=log, LT=LT)


def labeling_table_cd(r, CT: pd.DataFrame, samples: list) -> pd.DataFrame:
    """
    All isotopologues M+0…M+n (Relative Area 0…n) of the kept peak of every analysed metabolite, per sample.
    Fraction(M+i) = Relative Area i / Σ Relative Area 0…n — Compound Discoverer's "Exchange Rate in Percent".
    """
    from .io_reorganize import RA_COLS
    K = r.kept_rows()
    peak = r.W.drop_duplicates(["Method", "Compound"]).set_index(["Method", "Compound"]).Peak
    meta = CT.drop_duplicates("Compound").set_index("Compound")[["Formula", "Method", "n"]]
    rows = []
    for comp, m in meta.iterrows():
        pk = peak.get((m.Method, comp))
        g = K[(K.Mode == m.Method) & (K.Name == comp) & (K.PeakLabel == pk)].drop_duplicates("Sample").set_index("Sample")
        n = int(m.n)
        for i in range(n + 1):
            vals = {smp: (float(g.at[smp, RA_COLS[i]]) if smp in g.index and pd.notna(g.at[smp, RA_COLS[i]]) else np.nan)
                    for smp in samples}
            rows.append({"Compound": comp, "Formula": m.Formula, "Method": m.Method, "n": n, "i": i, **vals})
    LT = pd.DataFrame(rows)
    LT.insert(5, "Channel", LT.Compound + " | M+" + LT.i.astype(str))
    return LT


# ---------------------------------------------------------------------------
# Generic route
# ---------------------------------------------------------------------------
def _read_any(src, name: str) -> pd.DataFrame:
    name = name.lower()
    if name.endswith((".csv",)):
        return pd.read_csv(src)
    if name.endswith((".tsv", ".txt")):
        return pd.read_csv(src, sep="\t")
    return pd.read_excel(src)


def _find(cols, *cands):
    low = {c.lower().strip(): c for c in cols}
    for c in cands:
        if c in low:
            return low[c]
    return None


def parse_isotopologue(v) -> int | None:
    if pd.isna(v):
        return None
    s = str(v).strip().upper().replace(" ", "")
    m = re.fullmatch(r"(?:M\+?|M_|ISO|13C|C13)?(\d+)(?:\.0)?", s)
    if m:
        return int(m.group(1))
    m = re.search(r"M\+?(\d+)$", s)
    return int(m.group(1)) if m else None


def _type_of(v) -> str:
    s = str(v).strip().lower()
    if s in ("", "nan"):
        return TYPE_LAB
    if any(x in s for x in ("blank", "solvent", "blk")):
        return TYPE_BLANK
    if any(x in s for x in ("12c", "unlab", "reference", "natural", "ref", "control-unlabeled")):
        return TYPE_REF
    if s in ("qc", "pool", "pooled qc"):
        return "QC"
    return TYPE_LAB


def _time_parse(v, unit_default=""):
    if pd.isna(v) or str(v).strip() == "":
        return np.nan, unit_default
    s = str(v).strip().lower()
    m = re.fullmatch(r"(\d+(?:\.\d+)?)\s*(h|hr|hrs|hour|hours|min|mins|minute|minutes|m|d|day|days)?", s)
    if not m:
        return np.nan, unit_default
    u = m.group(2) or unit_default or "h"
    u = "h" if u.startswith("h") else ("d" if u.startswith("d") else "min")
    return float(m.group(1)), u


def read_metadata(src, name: str, samples: list) -> tuple[pd.DataFrame, list]:
    M = _read_any(src, name)
    log = []
    cs = _find(M.columns, "sample", "sample name", "sample id", "filename", "file")
    if cs is None:
        raise ValueError("The metadata file needs a 'Sample' column matching the sample columns of the data table.")
    cg = _find(M.columns, "group", "condition", "genotype", "cell line", "arm")
    ct = _find(M.columns, "time", "timepoint", "time point", "labeling time")
    cu = _find(M.columns, "timeunit", "time unit", "unit")
    ctr = _find(M.columns, "treatment", "drug", "intervention")
    cr = _find(M.columns, "replicate", "rep", "biological replicate")
    cty = _find(M.columns, "type", "sampletype", "sample type", "label status")
    ci = _find(M.columns, "include", "use")
    key = {str(s).strip().lower(): s for s in samples}
    used = {c for c in (cs, cg, ct, cu, ctr, cr, cty, ci) if c}
    extra = [c for c in M.columns if c not in used]
    rows, missing = [], []
    for _, r in M.iterrows():
        sname = key.get(str(r[cs]).strip().lower())
        if sname is None:
            missing.append(str(r[cs])); continue
        tv, tu = _time_parse(r[ct] if ct else np.nan, str(r[cu]).strip() if cu and pd.notna(r[cu]) else "")
        rows.append(dict(Sample=sname, SampleType=_type_of(r[cty]) if cty else TYPE_LAB,
                         Group=str(r[cg]).strip() if cg and pd.notna(r[cg]) else "All",
                         Treatment=str(r[ctr]).strip() if ctr and pd.notna(r[ctr]) else "",
                         TimeValue=tv, TimeUnit=tu, SampleNo=str(len(rows) + 1),
                         Replicate=int(r[cr]) if cr and pd.notna(r[cr]) else 0,
                         Include=bool(r[ci]) if ci and pd.notna(r[ci]) else True,
                         **{f"Meta: {c}": (str(r[c]) if pd.notna(r[c]) else "") for c in extra}))
    S = pd.DataFrame(rows)
    if missing:
        log.append(f"Metadata rows without a matching data column (ignored): {', '.join(missing[:10])}")
    absent = [s for s in samples if s not in set(S.Sample)]
    if absent:
        log.append(f"Data columns without metadata (excluded): {', '.join(absent[:10])}")
    S["Time"] = [f"{v:g} {u}" if pd.notna(v) else ("Unlabeled" if t == TYPE_REF else "—")
                 for v, u, t in zip(S.TimeValue, S.TimeUnit, S.SampleType)]
    return S, log


def read_generic(table_src, table_name: str, meta_src=None, meta_name: str | None = None,
                 tracer: str = "13C", min_area: float = 1000.0) -> Dataset:
    T = _read_any(table_src, table_name)
    T.columns = [str(c).strip() for c in T.columns]
    log = [f"Generic table: {table_name} — {len(T)} rows, {T.shape[1]} columns"]
    cm = _find(T.columns, "metabolite", "compound", "name", "feature")
    ciso = _find(T.columns, "isotopologue", "isotope", "label", "m+", "channel", "i")
    if cm is None or ciso is None:
        raise ValueError("The table needs a metabolite column (Metabolite / Compound / Name) and an isotopologue column "
                         "(Isotopologue: M+0, M+1, … or 0, 1, …).")
    cf = _find(T.columns, "formula")
    cmeth = _find(T.columns, "method", "mode")
    cn = _find(T.columns, "n", "labelatoms", "label atoms", "# c")
    csamp = _find(T.columns, "sample")
    cint = _find(T.columns, "intensity", "area", "value", "abundance")
    if csamp and cint:      # long format
        keys = [c for c in (cm, cf, cmeth, ciso, cn) if c]
        dup = T.duplicated(keys + [csamp]).sum()
        if dup:
            raise ValueError(f"{dup} duplicated rows for the same metabolite, isotopologue and sample. Each "
                             "isotopologue must appear once per sample (isotopologues are never summed).")
        T = T.pivot_table(index=keys, columns=csamp, values=cint, aggfunc="first").reset_index()
        T.columns = [str(c) for c in T.columns]
        log.append("Long format detected and pivoted to one column per sample.")
    samples = [c for c in T.columns if c.lower() not in NON_SAMPLE and c not in (cm, ciso, cf, cmeth, cn)
               and pd.api.types.is_numeric_dtype(pd.to_numeric(T[c], errors="coerce"))
               and pd.to_numeric(T[c], errors="coerce").notna().mean() > 0.3]
    if not samples:
        raise ValueError("No numeric sample columns were found.")
    T["_i"] = T[ciso].map(parse_isotopologue)
    bad = T["_i"].isna().sum()
    if bad:
        log.append(f"{bad} rows with an unreadable isotopologue label were skipped.")
    T = T[T["_i"].notna()].copy()
    rows, warn_n = [], []
    for comp, g in T.groupby(cm, sort=False):
        form = str(g[cf].iloc[0]) if cf and pd.notna(g[cf].iloc[0]) else ""
        meth = str(g[cmeth].iloc[0]) if cmeth and pd.notna(g[cmeth].iloc[0]) else "—"
        n = label_atoms(form, tracer)
        if not n and cn and pd.notna(g[cn].iloc[0]):
            n = int(g[cn].iloc[0])
        if not n:
            n = int(g["_i"].max()); warn_n.append(str(comp))
        if g["_i"].duplicated().any():
            raise ValueError(f"'{comp}' has the same isotopologue listed more than once. Give each isotopologue one "
                             "row per metabolite (isotopologues are never summed); use separate metabolite names for "
                             "separate peaks.")
        g = g.set_index("_i")[samples].apply(pd.to_numeric, errors="coerce")
        for i in range(0, n + 1):
            vals = g.loc[i] if i in g.index else pd.Series(np.nan, index=samples)
            rows.append({"Compound": str(comp), "Formula": form, "Method": meth, "n": n, "i": i,
                         **{s: float(vals[s]) if pd.notna(vals[s]) else np.nan for s in samples}})
    CT = pd.DataFrame(rows)
    CT.insert(5, "Channel", CT.Compound + " | M+" + CT.i.astype(str))
    if warn_n:
        log.append(f"{len(warn_n)} metabolites without a formula: number of labeled atoms taken as the highest "
                   "isotopologue present (atom enrichment is approximate for these).")
    ann = None
    acols = {c: ANNOT_COLS[c.lower()] for c in T.columns if c.lower() in ANNOT_COLS}
    if acols:
        ann = T[[cm] + list(acols)].drop_duplicates(cm).rename(columns={cm: "Compound", **acols})
    if meta_src is not None:
        S, mlog = read_metadata(meta_src, meta_name or "metadata.csv", samples)
        log += mlog
    else:
        S = parse_filenames(samples).rename(columns={"Filename": "_f"})
        S["Sample"] = samples
        S = S.drop(columns="_f")
        S["Treatment"] = ""
        log.append("No metadata file: groups and times were read from the sample names (edit them in Step 2).")
    keep = [c for c in ["Compound", "Formula", "Method", "n", "i", "Channel"]] + list(S.Sample)
    CT = CT[keep]
    full = CT.copy()                     # every isotopologue M+0…M+n as supplied (labeling denominator)
    # rows absent from the table = not reported: drop them instead of inventing zeros
    CT = CT[~CT[list(S.Sample)].isna().all(axis=1)].reset_index(drop=True)
    CT, zlog = apply_zero_rule(CT, list(S.loc[S.SampleType == TYPE_LAB, "Sample"]), min_area)
    log += zlog
    CT, mlog = drop_m0_only(CT)
    log += mlog
    S = order_samples(S.assign(TimeUnit=S.TimeUnit.fillna("")))
    LT = full[full.Compound.isin(set(CT.Compound))].reset_index(drop=True)
    log.append(f"Dataset: {CT.Compound.nunique()} metabolites, {len(CT)} isotopologue channels, {len(S)} samples.")
    log.append(f"Labeling table: all {len(LT)} isotopologues M+0…M+n of these metabolites — denominator of the "
               "isotopologue fractions (Area i / Σ Area 0…n, as Compound Discoverer's “Exchange Rate in Percent”).")
    return Dataset(CT=CT, S=S, source="generic", log=log, annotations=ann, LT=LT)


def to_generic(ds: Dataset) -> tuple[bytes, bytes]:
    """Export the channel table + metadata in the generic format (also used as templates)."""
    T = ds.CT.drop(columns=["Channel", "n"]).rename(columns={"Compound": "Metabolite", "i": "Isotopologue"})
    T["Isotopologue"] = "M+" + T.Isotopologue.astype(str)
    M = ds.S[["Sample", "SampleType", "Group", "Treatment", "TimeValue", "TimeUnit"]].rename(
        columns={"SampleType": "Type", "TimeValue": "Time"})
    M["Type"] = M.Type.map({TYPE_LAB: "Labeled", TYPE_REF: "Unlabeled reference", TYPE_BLANK: "Blank"}).fillna(M.Type)
    return T.to_csv(index=False).encode(), M.to_csv(index=False).encode()
