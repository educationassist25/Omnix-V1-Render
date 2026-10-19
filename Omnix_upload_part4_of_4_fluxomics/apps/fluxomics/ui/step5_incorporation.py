"""Step 5 — isotopologue distributions, fractional enrichment and incorporation."""
from __future__ import annotations

import numpy as np
import pandas as pd
import streamlit as st

from core import design as dz, plots as pl, correction as cr, labeling as lb
from ui import common as C

EXPLAIN = """
**Isotopologue fraction (distribution, MID)** — Fraction(M+i) = Area(M+i) / Σ Area(M+0…M+n) × 100, n = atoms of the
tracer element in the formula. This is Compound Discoverer's “Exchange Rate in Percent”: **every** isotopologue
M+0…M+n of the metabolite is used — also those the zero / minimum-value rule removed from the intensity analyses —
and an isotopologue not reported in a sample counts as 0. Fractions of one metabolite in one sample sum to 100 %.

**Incorporation (%)** — (1 − Fraction M+0) × 100: share of molecules carrying **at least one** label.

**Fractional enrichment, FE (%)** — Σ i·Fraction(M+i) / n × 100: share of the tracer-element atoms that are
labeled (atom percent enrichment). A pool that is 40 % M+2 and 60 % M+0 has 40 % incorporation but FE = 13 % for a
six-carbon metabolite.

**Mean label number** — Σ i·Fraction(M+i): average labels per molecule.

FE and mean label number reflect tracer-derived label only after natural-abundance correction.
"""


def render():
    C.header(4, "Isotope incorporation",
             "Where the label went: isotopologue distributions, fractional enrichment and incorporation for every "
             "metabolite, sample, arm and time point. These are ratios within a sample, so they do not depend on "
             "intensity normalization.")
    if not C.require("upload", "design"):
        return
    s = C.ss()
    D = C.design_final()
    F, S = C.labeling()
    An = dz.analysed(D)
    arms, times = dz.groups_of(D), dz.times_of(D)
    with st.expander("Definitions", expanded=False):
        st.markdown(EXPLAIN)
    bad = C.interference(S, D)
    chk = S[S.Sample.isin(An.Sample)].SumCheck
    c = st.columns(4)
    c[0].metric("Metabolite × sample profiles", int(S.Sample.isin(An.Sample).sum()))
    c[1].metric("Fractions sum to 1", f"{np.isclose(chk, 1).mean() * 100:.1f} %")
    c[2].metric("Labeled ≥ 5 % FE (any condition)",
                int((S[S.Sample.isin(An.Sample)].groupby("Compound").Enrichment.max() >= 5).sum()))
    c[3].metric("Possible interference", len(bad))
    mode = {"already": "as supplied (already corrected)", "apply": "corrected in this app", "none": "uncorrected"}
    meth = (f" — {cr.METHODS.get(s.correction.get('method', 'isocor_hr'), '')}"
            if s.correction["mode"] == "apply" else "")
    st.caption(f"Isotopologue data: {mode[s.correction['mode']]}{meth} · tracer {s.tracer['name']} "
               f"({s.tracer['element']}).")
    if bad:
        C.note(f"{len(bad)} metabolite(s) look labeled in the unlabeled reference (> 10 % incorporation), pointing to "
               "an isobaric interference rather than tracer label: " + ", ".join(bad[:12]) + ("…" if len(bad) > 12 else "")
               + ". They are excluded from labeling statistics, heatmaps and maps by default.")
    s.results["interference"] = bad
    t1, t2 = st.tabs(["Isotopologue distribution", "Isotope incorporation table"])
    with t1:
        comps = sorted(S.Compound.unique(), key=str.lower)
        default = "Citric acid" if "Citric acid" in comps else comps[0]
        comp = st.selectbox("Metabolite", comps, index=comps.index(default), key="mid_comp")
        hide0 = st.toggle("Hide M+0 (zoom into the labeled part)", key="hide_m0")
        C.show_plotly(pl.mid_bars_plotly(F, D, comp, arms, show_m0=not hide0), key="midbars")
        X = F[F.Compound == comp].merge(An[["Sample", "Condition"]], on="Sample")
        piv = X.pivot_table(index="i", columns="Condition", values="Fraction", aggfunc="mean")
        piv = piv.reindex(columns=[c for c in dz.conditions_of(D) if c in piv.columns])
        piv.index = [f"M+{i}" for i in piv.index]
        Sm = S[S.Compound == comp].merge(An[["Sample", "Condition"]], on="Sample").groupby("Condition")[
            ["Incorporation", "Enrichment", "MeanLabels"]].mean().reindex(piv.columns).T
        Sm.index = ["Incorporation (%)", "Fractional enrichment (%)", "Mean label number"]
        st.markdown("**Mean fraction (%) per condition, and summary metrics**")
        C.show_df(pd.concat([(piv * 100), Sm]).round(2))
        if s.correction["mode"] == "apply":
            cmp = X.groupby("i")[["Fraction_raw", "Fraction"]].mean() * 100
            cmp.index = [f"M+{i}" for i in cmp.index]
            st.markdown("**Corrected vs uncorrected (mean over included samples)**")
            C.show_df(cmp.rename(columns={"Fraction_raw": "Uncorrected (%)", "Fraction": "Corrected (%)"}).round(2))
        elif len(X) and X.Formula.iloc[0]:
            nat = cr.natural_mid(X.Formula.iloc[0], min(int(X.i.max()) + 1, 4), s.tracer["element"])
            st.caption("For reference, an unlabeled molecule with this formula shows "
                       + ", ".join(f"M+{i} {v * 100:.1f} %" for i, v in enumerate(nat))
                       + " before natural-abundance correction.")
        C.export_widget(lambda w, h: pl.mid_bars_mpl(F, D, comp, arms, not hide0, w, h), f"MID_{comp[:30]}",
                        "exp_mid", (6.5, 3.8), caption=f"Isotopologue distribution — {comp}",
                        tables={f"Mean fraction (%) and metrics — {comp}": pd.concat([(piv * 100), Sm]).round(4)})
    with t2:
        by_s, by_c = lb.incorporation_table(F[~F.Compound.isin(bad)], S[~S.Compound.isin(bad)], D)
        st.caption("Every metabolite with every isotopologue fraction (%), incorporation (%), fractional enrichment "
                   "(%) and mean label number — per sample, or as mean / SD / n per group × time point. "
                   + (f"{len(bad)} metabolite(s) with possible interference are excluded." if bad else ""))
        view = st.radio("Columns", ["Per sample", "Group × time (mean, SD, n)"], horizontal=True, key="inc_view")
        T = by_s if view == "Per sample" else by_c
        q = st.text_input("Filter metabolites", key="inc_q")
        TT = T[T.Compound.str.contains(q, case=False, regex=False)] if q else T
        C.show_df(TT.round(3), height=520)
        st.caption(f"{T.Compound.nunique()} metabolites · {len(T)} rows.")
        notes = pd.DataFrame({"Definition": [
            "Fraction M+i (%) = area(M+i) / Σ area(M+0…M+n, every isotopologue of the metabolite) × 100 "
            "(= Compound Discoverer “Exchange Rate in Percent”)",
            "Incorporation (%) = (1 − fraction M+0) × 100 = 100 − Fraction M+0 (%)",
            "Fractional enrichment (%) = Σ i × fraction(M+i) / n × 100 (n = carbon atoms)",
            "Mean label number = Σ i × fraction(M+i)",
            f"Natural-abundance correction: {s.correction['mode']}"
            + (f" ({cr.METHODS.get(s.correction.get('method', 'isocor_hr'), '')})" if s.correction["mode"] == "apply"
               else "") + f"; tracer {s.tracer['name']}"]})
        C.download_panel("Isotope_incorporation", "inc_dl",
                         tables={"Per sample": by_s, "Group × time (mean, SD, n)": by_c},
                         workbook={"Per_sample": by_s, "Group_x_time": by_c, "Definitions": notes})
    C.next_button("boxplots", "Dynamic bar plots")
