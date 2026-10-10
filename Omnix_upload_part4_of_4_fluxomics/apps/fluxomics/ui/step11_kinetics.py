"""Step 9 — label arrival: every group and time point side by side (isotopologues never summed)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import streamlit as st

from core import design as dz, kinetics as kn, pathways as pw, race as rc
from ui import common as C


def _values(F, D):
    """Replicate-level isotopologue fractions (%) with Key = 'Compound | M+i'."""
    An = dz.analysed(D)
    V = F.merge(An[["Sample", "Group", "Hours"]], on="Sample").assign(
        Channel=lambda d: "M+" + d.i.astype(str), Value=lambda d: d.Fraction * 100)
    V["Key"] = V.Compound + " | " + V.Channel
    return V.dropna(subset=["Hours"])


def render():
    C.header(9, "Label arrival",
             "Every isotopologue of every metabolite (M+1 … M+n, never summed) for every group and time point, side "
             "by side: one panel per group × time point. Labeling describes label appearance in a pool, not "
             "absolute flux.")
    if not C.require("upload", "design"):
        return
    s = C.ss()
    D = C.design_final()
    groups, times = dz.groups_of(D), dz.times_of(D)
    F, S = C.labeling()
    if len(times) < 2:
        st.warning("At least two time points are needed for label-arrival analysis."); return
    bad = set(C.interference(S, D))
    A = C.annotation()
    members = pw.pathway_members(A)
    V = _values(F[~F.Compound.isin(bad)], D)
    SK = kn.summarize(V.assign(Compound=V.Key))

    c = st.columns([3, 2, 2])
    default_pw = [x for x in ["Glycolysis", "TCA cycle", "Purine metabolism", "Pyrimidine metabolism"]
                  if x in members] or sorted(members)[:3]
    pw_sel = c[0].multiselect("Pathways", sorted(members), default=default_pw, key="race2_pw")
    gsel = c[1].multiselect("Groups", groups, default=groups, key="race_groups")
    tsel = c[2].multiselect("Time points", times, default=times, format_func=lambda h: f"{h:g} h", key="race_times")
    o = st.columns(4)
    unc = o[0].radio("Error bar", ["SEM", "95% CI"], horizontal=True, key="race_unc")
    thr = o[1].number_input("Show isotopologues reaching ≥ (%)", 0.0, 100.0, 0.0, 0.5, key="race_thr")
    maxn = o[2].slider("Max isotopologues per pathway", 1, 100, 12, key="race_max")
    m0 = o[3].toggle("Include M+0", False, key="race_m0")
    gsel = [g for g in groups if g in gsel]
    tsel = sorted(tsel)
    if not pw_sel or not gsel or not tsel:
        st.info("Choose at least one pathway, group and time point."); return
    mem = {k: members[k] for k in pw_sel}
    labels = {c_: C.label_of(c_) for cs in mem.values() for c_ in cs}
    R = rc.race_table(SK, mem, gsel, tsel, labels, thr, maxn, m0)
    if R.empty:
        st.info("No isotopologue in the selected pathways reaches the threshold."); return
    fig = rc.race_grid_plotly(R, gsel, tsel, unc, "Fraction (%)")
    C.show_plotly(fig, key="race_anim")
    s.results["race_done"] = True
    st.caption(f"{R.Label.nunique()} isotopologues · columns = {len(gsel)} groups × {len(tsel)} time points · bar "
               f"colour = pathway · error bar = {unc}.")
    W = R.assign(Column=R.Group + " · " + R.Hours.map(lambda h: f"{h:g} h")).pivot_table(
        index=["Pathway", "Metabolite", "Iso"], columns="Column", values="mean", aggfunc="first")
    W = W[[f"{g} · {h:g} h" for g in gsel for h in tsel if f"{g} · {h:g} h" in W.columns]].reset_index()
    W = W.rename(columns={"Iso": "Isotopologue"})
    with st.expander("Values (mean isotopologue fraction, %)"):
        C.show_df(W.round(3), height=360)
    n_rows = R.Label.nunique()
    C.export_widget(lambda w, h: rc.race_grid_mpl(R, gsel, tsel, unc, "Fraction (%)", width=w, height=h),
                    "label_arrival", "exp_race",
                    (max(8.0, 1.5 * len(gsel) * len(tsel) + 2.6), round(max(3.0, 0.15 * n_rows + 1.8), 1)),
                    caption="Label arrival — groups × time points side by side",
                    tables={"Mean fraction (%) per group × time": W,
                            "Mean, SEM, 95% CI (long)": R.drop(columns=["Key"])})
    st.caption("The exported figure is laid out in inches: the label column, title, group headers, time labels, "
               "axis label and legend each get their own band, and every row ≥ 0.15 in, so no text overlaps a bar "
               "or another label at any size (height grows automatically when needed).")
    C.next_button("export", "Export & report")
