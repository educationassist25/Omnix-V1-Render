"""Step 5 — dynamic bar plots: each metabolite with all isotopologues, all groups and time points in one plot."""
from __future__ import annotations

import numpy as np
import pandas as pd
import streamlit as st

from core import design as dz, barplots as bp
from ui import common as C

SOURCES = {
    "incorporation": "Isotope incorporation (%)",
    "fraction": "Isotopologue fraction (%)",
    "int": "Intensity — active dataset (sidebar)",
}


def _long(src, comp, F, An, S=None):
    if src == "incorporation":                # (1 − fraction M+0) × 100, one value per sample
        V = S[(S.Compound == comp) & S.Sample.isin(An.Sample)].assign(
            Value=lambda d: d.Incorporation, Isotopologue="Incorporation")
        return V[["Sample", "Isotopologue", "Value"]]
    if src == "fraction":
        V = F[(F.Compound == comp) & F.Sample.isin(An.Sample)].assign(
            Value=lambda d: d.Fraction * 100, Isotopologue=lambda d: "M+" + d.i.astype(str))
        return V[["Sample", "Isotopologue", "Value"]]
    X = C.intensity()                     # active dataset chosen in the sidebar
    rows = [r for r in X.index if r.rsplit(" | ", 1)[0] == comp]
    if not rows:
        return pd.DataFrame(columns=["Sample", "Isotopologue", "Value"])
    V = X.loc[rows].stack().rename("Value").reset_index()
    V.columns = ["Channel", "Sample", "Value"]
    V["Isotopologue"] = V.Channel.str.rsplit(" | ", n=1).str[1]
    return V[["Sample", "Isotopologue", "Value"]]


def render():
    C.header(5, "Dynamic bar plots",
             "Each metabolite in one plot: isotope incorporation (%), every isotopologue fraction (%) or every "
             "isotopologue's intensity (x = M+0 … M+n); one bar per group × time point (group colour, lighter = earlier time). Bar = mean of the biological replicates, "
             "error bar = SD or SEM, dots = individual replicates.")
    if not C.require("upload", "design", "normalize"):
        return
    s = C.ss()
    D = C.design_final()
    F, S = C.labeling()
    An = dz.analysed(D)
    arms, times = dz.groups_of(D), dz.times_of(D)
    c = st.columns([3, 2, 2])
    src = c[1].radio("Values", list(SOURCES), format_func=SOURCES.get, key="box_src")
    comps = sorted(S.Compound.unique(), key=str.lower)
    default = [x for x in ["Citric acid", "L-Glutamic acid"] if x in comps] or comps[:1]
    sel = c[0].multiselect("Metabolites (one plot each)", comps, default=default, key="box_sel")
    with c[2]:
        gsel = st.multiselect("Groups", arms, default=arms, key="box_groups")
        tsel = st.multiselect("Time points", times, default=times, format_func=lambda h: f"{h:g} h", key="box_times")
    o = st.columns(3)
    hide0 = o[0].toggle("Hide M+0", False, key="box_m0")
    pts = o[1].toggle("Show replicates", True, key="box_pts")
    err = o[2].radio("Error bar", list(bp.ERRORS), format_func=bp.ERRORS.get, horizontal=True, key="bar_err")
    if not sel or not gsel or not tsel:
        st.info("Choose at least one metabolite, group and time point.")
        return
    ylab = SOURCES[src] if src in ("fraction", "incorporation") else "Intensity — " + C.INT_LABEL[s.intensity_choice]
    xlab = "" if src == "incorporation" else "Isotopologue"
    for k, comp in enumerate(sel[:8]):
        V = _long(src, comp, F, An, S)
        if V.empty:
            st.warning(f"{comp}: no values for this selection.")
            continue
        X = V.merge(An[["Sample", "Group", "Hours", "TimeLabel"]], on="Sample").rename(columns={"Group": "Arm"})
        X = X[X.Arm.isin(gsel) & X.Hours.isin(tsel)]
        isos = sorted(X.Isotopologue.unique(), key=lambda x: int(x[2:]) if x.startswith("M+") else -1)
        if hide0 and src != "incorporation":
            isos = [i for i in isos if i != "M+0"]
            X = X[X.Isotopologue.isin(isos)]
        if not isos:
            st.warning(f"{comp}: nothing left to plot."); continue
        title = f"{comp} — {ylab}"
        gsel_o = [g for g in arms if g in gsel]
        C.show_plotly(bp.iso_bar_plotly(X, gsel_o, sorted(tsel), isos, ylab, title, pts, err, xlab=xlab), key=f"box_fig_{k}")
        T = bp.bar_summary(X, gsel_o, sorted(tsel), isos)
        with st.expander(f"Values — {comp}"):
            C.show_df(T.drop(columns="Hours").round(4))
        C.export_widget(lambda w, h, X=X, isos=isos, title=title: bp.iso_bar_mpl(X, gsel_o, sorted(tsel), isos, ylab,
                                                                             title, pts, err, w, h, xlab=xlab),
                        f"barplot_{comp[:30]}", f"exp_box_{k}", bp.default_size(len(isos), len(gsel_o) * len(tsel)),
                        caption=f"Bar plot — {title}",
                        tables={"Mean ± SD / SEM per group × time": T.drop(columns="Hours"),
                                "Replicate values": X.rename(columns={"Arm": "Group"})[
                                    ["Sample", "Group", "TimeLabel", "Isotopologue", "Value"]]})
    if len(sel) > 8:
        st.caption("Showing the first 8 metabolites.")
    s.results["box_done"] = True
    C.next_button("stats", "Statistics")
