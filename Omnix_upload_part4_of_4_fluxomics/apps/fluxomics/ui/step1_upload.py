"""Step 1 — upload peak-intensity data (Compound Discoverer SIL exports or a generic isotopologue table) and
sample metadata; inspect the unified isotopologue-channel table."""
from __future__ import annotations

import io
import json
from pathlib import Path

import pandas as pd
import streamlit as st

from core import design as dz, dataset as dsm
from ui import common as C

EXAMPLE = Path(__file__).resolve().parent.parent / "example_data"
GENERIC = EXAMPLE / "generic"
TEMPLATES = EXAMPLE / "templates"
CD_TEMPLATES = [TEMPLATES / "cd_sil" / "CD_SIL_Template_RP_Pos.xlsx", TEMPLATES / "cd_sil" / "CD_SIL_Template_HILIC_Neg.xlsx"]
GEN_TEMPLATES = [("Table with 12C reference (CSV)", TEMPLATES / "generic" / "Generic_Isotopologue_Table_with_12C_Ref.csv"),
                 ("Table without 12C reference (CSV)",
                  TEMPLATES / "generic" / "Generic_Isotopologue_Table_without_12C_Ref.csv"),
                 ("Sample metadata (CSV)", TEMPLATES / "generic" / "Generic_Sample_Metadata.csv"),
                 ("Instructions (Excel)", TEMPLATES / "generic" / "Generic_Template_Instructions.xlsx")]
XLSX = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"


def _template_buttons(items, key):
    """One download button per template file (files missing from the installation are skipped)."""
    items = [(lab, f) for lab, f in items if f.exists()]
    if not items:
        return
    cols = st.columns(len(items))
    for c, (lab, f) in zip(cols, items):
        c.download_button(f"⬇ {lab}", f.read_bytes(), f.name, "text/csv" if f.suffix == ".csv" else XLSX,
                          key=f"{key}_{f.stem}")

TEMPLATE_TABLE = pd.DataFrame({
    "Metabolite": ["Citrate"] * 7 + ["Lactate"] * 4,
    "Formula": ["C6H8O7"] * 7 + ["C3H6O3"] * 4,
    "Method": ["HILIC_Neg"] * 11,
    "KEGG": ["C00158"] * 7 + ["C00186"] * 4,
    "Isotopologue": [f"M+{i}" for i in range(7)] + [f"M+{i}" for i in range(4)],
    "Ctrl_12h_1": [5.1e7, 2.0e6, 1.4e7, 3.0e6, 2.1e6, 4.0e5, 1.0e5, 8.0e8, 1.0e7, 2.0e7, 3.1e8],
    "KO_12h_1": [4.0e7, 1.9e6, 1.8e7, 3.9e6, 2.6e6, 5.0e5, 1.2e5, 7.0e8, 9.0e6, 1.7e7, 3.6e8],
})
TEMPLATE_META = pd.DataFrame({"Sample": ["Ctrl_12h_1", "KO_12h_1"], "Group": ["Ctrl", "KO"], "Treatment": ["", ""],
                              "Time": [12, 12], "TimeUnit": ["h", "h"], "Replicate": [1, 1],
                              "Type": ["Labeled", "Labeled"], "Include": [True, True]})


def _accept(ds, payload_meta):
    s = C.ss()
    s.ds = ds
    s.files_meta = payload_meta
    s.design = dz.default_design(ds.S)
    C.reset_downstream("upload")
    s.report_figs = []
    st.rerun()


def render():
    C.header(1, "Upload and inspect data",
             "Load isotopologue peak areas and sample information. Every isotopologue (M+0 … M+n) of every metabolite "
             "is kept as its own channel; isotopologues are never summed.")
    s = C.ss()
    route = st.radio("Input format", ["Compound Discoverer SIL export(s)", "Generic isotopologue table + metadata"],
                     horizontal=True, key="in_route")
    if route.startswith("Compound"):
        _cd_route(s)
    else:
        _generic_route(s)
    ds = s.ds
    if ds is None:
        _explain()
        return
    _inspect(ds)
    C.next_button("design", "Experimental design")


def _cd_route(s):
    files = st.file_uploader("Excel / CSV exports (columns Name, Formula, … Relative Area 0–25), one per method",
                             type=["xlsx", "xls", "csv", "tsv", "txt"], accept_multiple_files=True, key="cd_files")
    with st.expander("Format, templates and example", expanded=not files):
        st.markdown(
            "**One Compound Discoverer SIL export per LC-MS method** (Excel or CSV; the workbook name contains the "
            "method, e.g. `…_RP_Pos.xlsx`, `…_HILIC_Neg.xlsx`; data on the first sheet).  \n"
            "**Columns:** `Name` · `Formula` · `Molecular Weight` · `RT in min` · `Area` · `Status` · `Study File ID` · "
            "`Filename` · `Relative Area 0` … `Relative Area 25` (one row per compound × raw file).  \n"
            "**Raw-file names:** `<date>_<method>_<CellLine>_<Group>_<SampleNo>_<Time>.raw` (e.g. "
            "`…_CellLine_GeneKD_1001_12Hr.raw`); time `12C` = unlabeled reference; `Blank` = blank. The templates "
            "contain a sample list and full instructions on their extra sheets.")
        _template_buttons([("RP_Pos template (Excel)", CD_TEMPLATES[0]),
                           ("HILIC_Neg template (Excel)", CD_TEMPLATES[1])], "cd_tpl")
        cd_example = st.button("Load Compound Discoverer example (RP_Pos + HILIC_Neg templates; shControl vs GeneKD, "
                               "12 / 18 / 40 h)", key="cd_example", disabled=not all(f.exists() for f in CD_TEMPLATES),
                               help="Processes the two template workbooks exactly like uploaded exports.")
    s.min_area = float(st.number_input(
        "Minimum value in labeled samples", 0.0, 1e9, float(s.min_area), 100.0, key="min_area_in",
        help="An isotopologue with a value below this in ANY labeled sample is removed; if M+0 is below it, the "
             "whole peak is removed. 12C reference and blank are not checked; empty cells do not count. 0 = off."))
    payload = None
    if cd_example:
        payload = tuple((f.name, f.read_bytes()) for f in CD_TEMPLATES)
    elif files and st.button("Process files", type="primary", key="cd_go"):
        payload = tuple((f.name, f.getvalue()) for f in files)
    if payload:
        try:
            with st.spinner("Reading files, matching peaks across samples, selecting one peak per metabolite…"):
                r = C.run_reorg(payload, s.keep_partial, json.dumps(s.prefer), s.min_area)
                ds = dsm.from_reorg(r, s.tracer["element"])
        except Exception as e:
            st.error(f"Could not process the files: {e}")
            return
        _accept(ds, [(n, C.sha(b), None) for n, b in payload])


def _generic_route(s):
    c1, c2 = st.columns(2)
    tab = c1.file_uploader("Isotopologue intensity table (CSV / TSV / Excel)", type=["csv", "tsv", "txt", "xlsx", "xls"],
                           key="gen_table")
    meta = c2.file_uploader("Sample metadata (optional but recommended)", type=["csv", "tsv", "txt", "xlsx", "xls"],
                            key="gen_meta")
    with st.expander("Format and templates", expanded=not tab):
        st.markdown(
            "**Wide table:** `Metabolite` · `Formula` (recommended: sets n = number of labeled atoms) · `Method` · "
            "optional IDs `KEGG` `HMDB` `PubChemCID` `ChEBI` `BiGG` · `Isotopologue` (M+0, M+1 … or 0, 1 …) · "
            "one column per sample with peak areas.  \n"
            "**Long table** also works: `Metabolite`, `Isotopologue`, `Sample`, `Intensity`.  \n"
            "**Metadata:** `Sample` (matching the column names) · `Group` · `Time` · `TimeUnit` (min / h / d) · "
            "`Treatment` · `Replicate` · `Type` (Labeled / Unlabeled reference / Blank / QC) · `Include`.  \n"
            "Missing isotopologue rows count as not detected; empty cells are missing values (see Normalization).")
        if all(f.exists() for _, f in GEN_TEMPLATES):
            _template_buttons(GEN_TEMPLATES, "gen_tpl")
        else:   # fallback: minimal templates built in
            d1, d2 = st.columns(2)
            d1.download_button("⬇ Table template (CSV)", TEMPLATE_TABLE.to_csv(index=False).encode(),
                               "isotopologue_table_template.csv", "text/csv")
            d2.download_button("⬇ Metadata template (CSV)", TEMPLATE_META.to_csv(index=False).encode(),
                               "sample_metadata_template.csv", "text/csv")
    ex = GENERIC / "synthetic_isotopologue_table.csv"
    use_example = st.button("Load synthetic example data (shControl vs GeneKD, [U-13C6]-glucose, 2 / 8 / 24 h)",
                            disabled=not ex.exists(), key="gen_example",
                            help="Simulated data (no real experiment): 63 metabolites, 2 groups × 3 time points × 3 "
                                 "replicates, one unlabeled reference and one blank.")
    if use_example:
        tb, mb = ex.read_bytes(), (GENERIC / "synthetic_sample_metadata.csv").read_bytes()
        names = (ex.name, "synthetic_sample_metadata.csv")
    elif tab is not None and st.button("Read table", type="primary", key="gen_go"):
        tb, mb = tab.getvalue(), (meta.getvalue() if meta is not None else None)
        names = (tab.name, meta.name if meta is not None else None)
    else:
        return
    try:
        with st.spinner("Reading table…"):
            ds = dsm.read_generic(io.BytesIO(tb), names[0], io.BytesIO(mb) if mb else None, names[1],
                                  s.tracer["element"], min_area=s.min_area)
    except Exception as e:
        st.error(f"Could not read the table: {e}")
        return
    files = [(names[0], C.sha(tb), None)] + ([(names[1], C.sha(mb), None)] if mb else [])
    _accept(ds, files)


def _explain():
    st.markdown("#### What happens to your data")
    st.markdown(
        "- **Compound Discoverer exports:** each “Relative Area i” column becomes isotopologue channel M+i (up to "
        "the number of tracer-element atoms). Peaks are matched across samples by retention time; peaks whose M+0 is "
        "0 (or below the minimum value, default 1,000) in any labeled sample are removed; isotopologue rows "
        "below the minimum value in any labeled sample are removed; each metabolite is kept **once** (best peak / name / method); "
        "metabolites left with only M+0 are excluded.\n"
        "- **Generic table:** read as is; metadata assign groups, treatments and time points.\n"
        "- **One internal table:** metabolite × isotopologue channels × samples. Intensity analyses use each channel "
        "separately; fractions and incorporation are calculated within each metabolite.\n"
        "- **No values are changed** at this step.")


def _inspect(ds):
    s = C.ss()
    S = ds.S
    m = st.columns(5)
    m[0].metric("Source", "Compound Discoverer" if ds.source == "cd" else "Generic table")
    m[1].metric("Metabolites", ds.CT.Compound.nunique())
    m[2].metric("Isotopologue channels", len(ds.CT))
    m[3].metric("Samples", len(S))
    lab = (S.SampleType == "13C labeled").sum() if "SampleType" in S else len(S)
    m[4].metric("Labeled samples", int(lab))
    with st.expander("Processing log", expanded=False):
        for x in ds.log:
            st.markdown(f"- {x}")
        if ds.reorg is not None:
            for x in ds.reorg.dropped_peaks[:200]:
                st.markdown(f"- Peak removed (M+0 = 0 in a labeled sample): {x}")
    tabs = ["Isotopologue channels", "Samples", "Missing values"] + \
           (["Combined table (clean)", "Selection decisions", "Removed peaks"] if ds.reorg is not None else [])
    T = st.tabs(tabs)
    with T[0]:
        q = st.text_input("Filter metabolites", key="f1")
        CT = ds.CT if not q else ds.CT[ds.CT.Compound.str.contains(q, case=False, regex=False)]
        C.show_df(CT, height=420)
    with T[1]:
        C.show_df(S, height=420)
    with T[2]:
        X = ds.channel_matrix()
        miss = X.isna() | (X <= 0)
        per_s = (miss.mean() * 100).rename("Missing or zero (%)").to_frame()
        per_s["Type"] = S.set_index("Sample").SampleType.reindex(per_s.index)
        c1, c2 = st.columns([1, 2])
        with c1:
            C.show_df(per_s.round(1), height=380)
        with c2:
            import plotly.graph_objects as go
            by_i = miss.groupby(ds.CT.set_index("Channel").i.reindex(X.index).values).mean().mean(axis=1) * 100
            fig = go.Figure(go.Bar(x=[f"M+{int(i)}" for i in by_i.index], y=by_i.values, marker_color="#6f5fc7"))
            fig.update_layout(template="plotly_white", height=360, title="Missing / zero values by isotopologue",
                              yaxis_title="% of cells", margin=dict(l=50, r=10, t=50, b=40))
            C.show_plotly(fig, key="miss_bar")
        st.caption("Zeros in Compound Discoverer “0Fill” exports usually mean “not detected”; Normalization lets you treat "
                   "them as missing.")
    if ds.reorg is not None:
        r = ds.reorg
        with T[3]:
            st.caption(f"One peak per metabolite, built with the rules of sil_flux_pipeline.py "
                       f"({'FILE_PATTERN sample naming' if r.naming == 'reference' else 'generic sample naming'}). "
                       f"{len(r.m0_only)} metabolite(s) with only M+0 were excluded (listed in Removed peaks). "
                       "This is exactly the table used by every later step.")
            C.show_df(r.W, height=420)
            st.download_button("⬇ Combined workbook (.xlsx)", r.workbook_bytes(),
                               "Fluxomics_combined_one_peak_per_metabolite.xlsx",
                               "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
        with T[4]:
            C.show_df(r.decisions, height=420)
        with T[5]:
            C.show_df(r.W_removed, height=420)
    t, mtab = dsm.to_generic(ds)
    c1, c2, _ = st.columns([1, 1, 2])
    c1.download_button("⬇ Channel table (generic CSV)", t, "isotopologue_channel_table.csv", "text/csv")
    c2.download_button("⬇ Sample metadata (CSV)", mtab, "sample_metadata.csv", "text/csv")
