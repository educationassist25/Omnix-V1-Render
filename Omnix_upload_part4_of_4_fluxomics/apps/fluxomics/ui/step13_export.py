"""Step 10 — export processed data, statistics, figures and a reproducible report."""
from __future__ import annotations

import json

import pandas as pd
import streamlit as st

from core import design as dz, labeling as lb, report as rp
from ui import common as C


def _settings(s):
    return {"input": s.ds.source if s.ds is not None else None, "tracer": s.tracer, "correction": s.correction,
            "zero_rule": "isotopologue rows with 0 or < minimum value in any labeled sample removed; M+0 = 0 or "
                         "< minimum value removes the compound; metabolites with only M+0 excluded",
            "normalization": s.norm["settings"] if s.norm else None,
            "intensity_dataset": s.intensity_choice, "minimum_value_labeled": s.min_area,
            "comparison": {"A": s.group_a, "B": s.group_b},
            "statistics_runs": {k: v["settings"] for k, v in s.stats.items()}}


def render():
    C.header(10, "Export and report",
             "Download processed data and every statistics table in one workbook, all collected figures, and a "
             "self-contained HTML report with checksums, settings, methods and software versions.")
    if not C.require("upload", "design"):
        return
    s = C.ss()
    D = C.design_final()
    t1, t2, t3 = st.tabs(["Data & statistics", "Figures", "Report & settings"])
    with t1:
        data_tab(s, D)
    with t2:
        figs_tab(s)
    with t3:
        report_tab(s, D)


def data_tab(s, D):
    F, S = C.labeling()
    An = dz.analysed(D)
    samples = An.Sample.tolist()
    sheets = {"Design": D, "Channels_raw_areas": s.ds.CT}
    if s.ds.LT is not None:                      # every isotopologue M+0…M+n: denominator of the fractions
        sheets["Labeling_areas_all_isotopologues"] = s.ds.LT
    if s.norm:
        sheets["Intensity_log2"] = s.norm["L"].reset_index()
        sheets["Intensity_log2_IQR"] = s.norm["N"].reset_index()
        sheets["IQR_normalization"] = s.norm["info"].reset_index()
    sheets["Fractions_long"] = F
    sheets["Labeling_per_sample"] = S
    for m in ("incorporation", "enrichment", "mean_labels"):
        sheets[f"{m}_matrix"] = lb.metric_matrix(S, m, samples).reset_index()
    for k, (name, res) in enumerate(s.stats.items(), 1):
        sheets[f"Stats{k}_{name}"[:31]] = res["table"]
        if res.get("posthoc") is not None and len(res["posthoc"]):
            sheets[f"Stats{k}_posthoc"] = res["posthoc"]
    if s.annotation is not None:
        sheets["Annotation"] = s.annotation
    st.caption("Sheets: " + ", ".join(sheets))
    if s.stats:
        st.caption("Statistics runs included: " + "; ".join(f"Stats{k} = {n}" for k, n in enumerate(s.stats, 1)))
    if st.button("Prepare Excel workbook", type="primary", key="xl_prep"):
        with st.spinner("Writing workbook…"):
            st.session_state["xl_bytes"] = C.to_excel(sheets)
            s.results["exported"] = True
    if st.session_state.get("xl_bytes"):
        st.download_button("⬇ Processed data & statistics (.xlsx)", st.session_state["xl_bytes"],
                           "Fluxomics_processed_data_and_statistics.xlsx",
                           "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
    if st.button("Prepare CSV bundle (ZIP)", key="csv_prep"):
        files = {}
        for n, df in sheets.items():
            if df is not None and len(df):
                files[f"{n}.csv"] = df.to_csv(index=not isinstance(df.index, pd.RangeIndex)).encode()
        st.session_state["csv_zip"] = C.zip_bytes(files)
        s.results["exported"] = True
    if st.session_state.get("csv_zip"):
        st.download_button("⬇ All tables as CSV (.zip)", st.session_state["csv_zip"], "Fluxomics_tables_csv.zip",
                           "application/zip")
    if s.ds.reorg is not None:
        st.download_button("⬇ Reorganized source workbook (Step 1)", s.ds.reorg.workbook_bytes(),
                           "Fluxomics_reorganized_one_peak_per_metabolite.xlsx",
                           "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")


def figs_tab(s):
    figs = s.report_figs
    if not figs:
        st.info("Figures you render in the “⬇ Downloads” box of Steps 4–9 are collected here and in the "
                "report.")
        return
    st.caption(f"{len(figs)} figure(s): " + "; ".join(c for c, _ in figs))
    st.download_button("⬇ All collected figures (PNG, ZIP)",
                       C.zip_bytes({f"{i + 1:02d}_{c[:60].replace('/', '-')}.png": b for i, (c, b) in enumerate(figs)}),
                       "Fluxomics_figures.zip", "application/zip")
    for cap, png in figs[-6:]:
        st.caption(cap)
        C.show_image(png)
    if st.button("Clear collected figures", key="fig_clear"):
        s.report_figs = []
        st.rerun()


def report_tab(s, D):
    st.markdown("#### Reproducible analysis report (HTML)")
    st.caption("Input checksums, processing log, design, all settings, key result tables, collected figures, the "
               "software versions and a methods paragraph.")
    settings = _settings(s)
    if st.button("Build report", type="primary", key="rep_btn"):
        tables = {}
        for name, res in s.stats.items():
            T = res["table"]
            tables[f"Statistics — {name}"] = T.sort_values("p") if "p" in T else T
        if s.norm:
            tables["IQR normalization factors"] = s.norm["info"].reset_index()
        A = C.annotation()
        tables["Annotation (mapping confidence)"] = A[["Compound", "Formula", "Confidence", "StandardName", "KEGG",
                                                       "HMDB", "Pathways"]]
        html = rp.build_report(title="Stable-isotope tracing analysis report",
                               inputs=[(n, h, "") for n, h, _ in s.files_meta], design=D, settings=settings,
                               log=s.ds.log, tables=tables, figures=s.report_figs)
        st.session_state["report_html"] = html
        s.results["exported"] = True
    if st.session_state.get("report_html"):
        st.download_button("⬇ Download report (.html)", st.session_state["report_html"].encode(),
                           "Fluxomics_analysis_report.html", "text/html")
    st.markdown("#### Methods paragraph")
    st.text_area("Copy into your manuscript (edit as needed)", rp.methods_text(settings), height=220, key="methods_txt")
    st.markdown("#### Settings file")
    st.download_button("⬇ Settings (.json)", rp.settings_json({**settings, "design": D.to_dict("records")}).encode(),
                       "Fluxomics_settings.json", "application/json")
    up = st.file_uploader("Load settings (.json)", type=["json"], key="set_up")
    if up is not None and st.button("Apply settings", key="set_apply"):
        cfg = json.load(up)
        s.tracer.update(cfg.get("tracer") or {})
        s.correction.update(cfg.get("correction") or {})
        if cfg.get("normalization"):
            s.norm_settings.update({k: v for k, v in cfg["normalization"].items() if k in s.norm_settings})
        s.min_area = float(cfg.get("minimum_value_labeled", s.min_area))
        if cfg.get("intensity_dataset") in ("log2", "log2_iqr"):
            s.intensity_choice = cfg["intensity_dataset"]
        if cfg.get("design"):
            Dn = pd.DataFrame(cfg["design"])
            keep = [c for c in ["Sample", "SampleType", "Group", "Treatment", "TimeValue", "TimeUnit", "Include",
                                "Replicate"] if c in Dn and c in s.design]
            cur = s.design.set_index("Sample")
            cur.update(Dn[keep].set_index("Sample"))
            s.design = cur.reset_index()
        C.reset_downstream("design")
        st.success("Settings applied. Re-run Step 3 (normalization) and the analysis steps.")
