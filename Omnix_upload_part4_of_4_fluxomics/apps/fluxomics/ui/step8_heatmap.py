"""Step 8 — interactive and publication heatmaps of isotopologue-channel intensities.

Columns: any subset of conditions (groups) × time points; the row-wise Z-score is computed across the samples shown.
Rows: the entire dataset or the metabolites of selected annotation pathways (one combined heatmap or one heatmap
per pathway), then filtered by significance (Step 6), variance or all. Rows can be grouped by metabolite: short
isotopologue labels (M+0 … M+n, '(n.s.)' when not significant) with a bracket and the metabolite name."""
from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import streamlit as st

from core import design as dz, heatmap as hm, pathways as pw, plots as pl
from ui import common as C


def _col_info(D, cols):
    """Annotation table for the heatmap columns (samples): group, time, treatment, replicate, metadata columns."""
    s = C.ss()
    I = D.set_index("Sample")
    extra = [c for c in s.ds.S.columns if c.startswith("Meta: ")]
    if extra:
        I = I.join(s.ds.S.set_index("Sample")[extra])
    I = I.loc[cols]
    return I.assign(Time=I.TimeLabel)


def _comp(r) -> str:
    return str(r).rsplit(" | ", 1)[0]


def _hours_txt(hs) -> str:
    return " & ".join(f"{h:g}" for h in hs) + " h"


def render():
    C.header(8, "Heatmaps",
             "Clustered heatmaps of individual isotopologue channels (active intensity dataset from the sidebar, row-wise "
             "Z-score, one column per sample). Choose the conditions and time points to compare, the whole dataset or "
             "selected pathways, and select rows by p value or FDR from a Step 6 result.")
    if not C.require("upload", "design", "normalize"):
        return
    # heatmaps always sit in the middle of the page, whatever their width
    st.markdown("<style>[class*='st-key-hm_plotly'], div[data-testid='stElementContainer']:has([data-testid="
                "'stPlotlyChart']), .element-container:has(.stPlotlyChart), [data-testid='stPlotlyChart']"
                "{align-self:center !important;margin-left:auto !important;margin-right:auto !important}</style>",
                unsafe_allow_html=True)
    s = C.ss()
    D = C.design_final()
    An = dz.analysed(D)
    arms = dz.groups_of(D)
    times = dz.times_of(D)
    X, info = C.layer("int")
    A = C.annotation()

    # ---- columns (conditions × time points) and pathway scope
    with st.expander("Conditions, time points and pathways", expanded=True):
        a = st.columns(2)
        gsel = a[0].multiselect("Conditions (groups)", arms, default=arms, key="hm_groups",
                                help="E.g. Control and KD. The Z-score is recomputed across the samples shown.")
        tsel = a[1].multiselect("Time points", times, default=times, format_func=lambda h: f"{h:g} h",
                                key="hm_times", help="E.g. 12 & 18 h, 12 & 40 h or 18 & 40 h.")
        present = {_comp(r) for r in X.index}
        members = {p: [c for c in dict.fromkeys(v) if c in present] for p, v in pw.pathway_members(A).items()}
        members = {p: v for p, v in members.items() if v}
        b = st.columns([1.3, 3.2, 1.5])
        scope = b[0].radio("Metabolites", ["all", "pw"], key="hm_scope",
                           format_func={"all": "Entire dataset", "pw": "Selected pathways"}.get,
                           help="Pathways come from the metabolite annotation (High, Medium and user-confirmed "
                                "assignments; edit them in Step 7 · Pathway mapping).")
        psel = b[1].multiselect("Pathways", sorted(members), key="hm_pws", disabled=scope == "all",
                                format_func=lambda p: f"{p}  ·  {len(members[p])} metabolites",
                                placeholder="Choose one or more pathways")
        layout = b[2].radio("Layout", ["combined", "each"], key="hm_pw_layout", disabled=scope == "all",
                            format_func={"combined": "One combined heatmap",
                                         "each": "One heatmap per pathway"}.get)
        c = st.columns([2, 1.4, 2.6])
        can_sep = not (scope == "pw" and layout == "each")
        sep_on = c[0].toggle("Separate rows by pathway", scope == "pw", key=f"hm_pw_sep_{scope}",
                             disabled=not can_sep,
                             help="Each pathway is its own block of rows (never mixed), labelled on the left with "
                                  "a vertical line and the pathway name. Clustering orders rows only within a "
                                  "pathway. A metabolite in two selected pathways appears in both blocks. "
                                  "Entire dataset: each metabolite's first assigned pathway, plus “Other / not "
                                  "mapped”.") and can_sep
        pw_gap = c[1].number_input("Gap between pathways (rows)", 0.0, 4.0, 0.5, 0.25, key="hm_pw_gap",
                                   disabled=not sep_on,
                                   help="Blank space between pathway blocks, in row heights (0 = no gap).")
    gsel = [g for g in arms if g in gsel]
    tsel = sorted(tsel)
    cols_sel = An[An.Group.isin(gsel) & An.Hours.isin(tsel)].Sample.tolist()
    cols_sel = [c for c in cols_sel if c in X.columns]
    if len(cols_sel) < 3:
        st.warning("Choose conditions and time points with at least three samples in total."); return
    if scope == "pw" and not psel:
        st.info("Choose at least one pathway, or switch to “Entire dataset”."); return
    M_all = X[cols_sel]

    # ---- rows: significance (Step 6), variance or all
    with st.expander("Rows to show", expanded=True):
        r = st.columns([2, 2, 1, 1])
        mode = r[0].radio("Select rows by", ["sig", "var", "all"], key="hm_rows",
                          format_func={"sig": "Significance (p value / FDR, Step 6)", "var": "Top N by variance",
                                       "all": "All channels"}.get)
        sig_set = tested = None
        thr_txt = rs = ""
        if s.stats:
            rs = r[1].selectbox("Statistics result", list(s.stats), index=len(s.stats) - 1, key="hm_res")
            T = s.stats[rs]["table"]
            comps = list(dict.fromkeys(T.Comparison))
            csel = r[1].multiselect("Comparisons (significant in any)", comps, default=comps, key="hm_cmps")
            crit = r[2].radio("Criterion", ["FDR", "P value", "P value and FDR"], key="hm_crit")
            thr_q = r[3].number_input("FDR threshold (≤ 1)", 0.0001, 1.0, 0.25, 0.01, format="%.4f",
                                      key="hm_thr_fdr", disabled=crit == "P value")
            thr_p = r[3].number_input("P value threshold (≤ 1)", 0.0001, 1.0, 0.05, 0.01, format="%.4f",
                                      key="hm_thr_p", disabled=crit == "FDR")
            TT = T[T.Comparison.isin(csel)]
            qc = ["FDR"] if "FDR" in TT else [x for x in TT.columns if x.startswith("FDR ")]
            pc_ = ["p"] if "p" in TT else [x for x in TT.columns if x.startswith("p ")]
            if crit == "FDR":
                hit = (TT[qc] <= thr_q).any(axis=1)
                thr_txt = f"FDR ≤ {thr_q:g}"
            elif crit == "P value":
                hit = (TT[pc_] <= thr_p).any(axis=1)
                thr_txt = f"p ≤ {thr_p:g}"
            else:
                hit = (TT[pc_] <= thr_p).any(axis=1) & (TT[qc] <= thr_q).any(axis=1)
                thr_txt = f"p ≤ {thr_p:g} and FDR ≤ {thr_q:g}"
            sig_set, tested = set(TT[hit].Feature), set(TT.Feature)
        elif mode == "sig":
            r[1].warning("Run Step 6 first; showing top 60 by variance.")
        top = r[1].number_input("Top N rows by variance", 2, 5000, 60, 10, key="hm_top") if mode == "var" else 60
        inc_m0 = r[1].toggle("Include M+0 channels", True, key="hm_m0")
        ns_show = r[2].radio("Non-significant rows (n.s.)", ["hide", "show"], key="hm_ns_rows",
                             disabled=sig_set is None or mode != "sig",
                             format_func={"hide": "Hide (significant only)",
                                          "show": "Show, marked “(n.s.)”"}.get,
                             help="Show: every channel is kept; channels that do not meet the criterion in any "
                                  "selected comparison are labelled e.g. “M+2 (n.s.)” (on screen and in the "
                                  "export).") == "show"
        ns_mark = ns_show and mode == "sig" and sig_set is not None
        pick_box = st.container()

    # ---- panels (one heatmap, or one per pathway)
    if scope == "pw":
        if layout == "each":
            panels = [(p, set(members[p])) for p in psel]
        else:
            panels = [(", ".join(psel) if len(psel) <= 3 else f"{len(psel)} pathways",
                       set().union(*[members[p] for p in psel]))]
    else:
        panels = [(None, None)]
    prepared = []
    for pname, cset in panels:
        M = M_all if cset is None else M_all.loc[[i for i in M_all.index if _comp(i) in cset]]
        if not inc_m0:
            M = M[~M.index.str.endswith("| M+0")]
        M = M.dropna(how="all")
        M = M.loc[M.std(axis=1, ddof=1).fillna(0) > 0]          # constant rows have no Z-score
        msg = ""
        if ns_mark:
            msg = (f"{len(M)} channels; {sum(1 for f in M.index if f in tested and f not in sig_set)} marked "
                   f"(n.s.) — not {thr_txt} in any selected comparison of “{rs}”.")
        elif mode == "sig" and sig_set is not None:
            Ms = M.loc[[f for f in M.index if f in sig_set]]
            if len(Ms) >= 2:
                M, msg = Ms, f"{len(Ms)} channels with {thr_txt} in at least one selected comparison of “{rs}”."
            else:
                msg = (f"Only {len(Ms)} channel(s) have {thr_txt} in “{rs}”"
                       + (f" for {pname}" if pname else "") + "; a heatmap needs ≥ 2 rows. Showing the 60 most "
                       "variable channels instead — try “P value”, a higher threshold or “All channels”.")
                M = M.loc[M.var(axis=1).sort_values(ascending=False).index[:60]]
        elif mode in ("var", "sig"):
            M = M.loc[M.var(axis=1).sort_values(ascending=False).index[:int(top)]]
        prepared.append((pname, M, msg))
    all_rows = sorted({i for _, M, _ in prepared for i in M.index})
    with pick_box:
        pick = st.multiselect("Highlight rows", all_rows, key="hm_pick")

    ci0 = _col_info(D, cols_sel)
    base_tracks = [t for t in ["Group", "Time", "Treatment", "Replicate"] if t in ci0.columns
                   and ci0[t].astype(str).str.strip().replace({"0": "", "nan": ""}).ne("").any()]
    meta_tracks = [t for t in ci0.columns if t.startswith("Meta: ")]
    with st.expander("Column annotations, clustering and colours", expanded=False):
        rl = st.columns([2, 1, 1])
        row_style = rl[0].radio("Row label style", ["full", "grouped"], horizontal=True, key="hm_row_style",
                                format_func={"full": "Metabolite | M+i", "grouped": "Grouped by metabolite"}.get,
                                help="Grouped by metabolite: the rows of each metabolite stay together (M+0, M+1 …), "
                                     "labelled M+i; the metabolite name sits in its own column, joined to its rows "
                                     "by a plain vertical line. Metabolites follow the clustering order (no row "
                                     "dendrogram, since clustering single rows would split a metabolite apart).")
        row_labels = rl[1].toggle("Show row labels", True, key="hm_rl")
        seps = rl[2].toggle("Lines between metabolites", True, key="hm_seps", disabled=row_style != "grouped")
        tracks = st.multiselect("Column annotations (sample group and metadata)", base_tracks + meta_tracks,
                                default=[t for t in ["Group", "Time"] if t in base_tracks], key="hm_tracks")
        z = st.columns(3)
        cr_ = z[0].toggle("Cluster rows", True, key="hm_cr",
                          help="Grouped labels: clusters metabolites by their mean channel profile.")
        cc_ = z[1].toggle("Cluster columns", False, key="hm_cc")
        k = st.columns(3)
        metric = k[0].selectbox("Distance", list(hm.DISTANCES), format_func=hm.DISTANCES.get, key="hm_dist")
        method = k[1].selectbox("Linkage", list(hm.LINKAGES), format_func=hm.LINKAGES.get, key="hm_link")
        cmap = k[2].selectbox("Colour scale", list(hm.CMAPS), key="hm_cmap")
        ap = st.columns(3)
        auto = ap[0].toggle("Automatic colour limits", True, key="hm_auto")
        lim = ap[1].number_input("Colour limit ±", 0.1, 100.0, 3.0, 0.5, key="hm_lim", disabled=auto)
    grouped = row_style == "grouped"

    n0 = max(len(prepared[0][1]), 2)
    bw, bh = hm.balanced_cells(n0, len(cols_sel))
    with st.expander("Figure size and fonts", expanded=False):
        sizing = st.radio("Proportions", ["balanced", "cells", "custom"], horizontal=True, key="hm_sizing",
                          format_func={"balanced": "Balanced (automatic, proportional)",
                                       "cells": "Fixed cell size", "custom": "Custom figure size"}.get,
                          help="Balanced: cell width and height follow the number of samples and rows so the "
                               "heatmap body keeps publication proportions (near-square cells for small heatmaps, "
                               "≤ 7 in wide for large ones). Fixed: type the cell size. Custom: type the figure "
                               "size; cells stretch to fill it.")
        f0 = st.columns(4)
        cell_w_in = f0[0].number_input("Cell width (in)", 0.08, 0.80, 0.30, 0.02, format="%.2f", key="hm_cell_w",
                                       disabled=sizing != "cells", help="Width of one sample column.")
        cell_h_in = f0[1].number_input("Cell height (in)", 0.06, 0.30, 0.10, 0.01, format="%.2f", key="hm_cell",
                                       disabled=sizing != "cells",
                                       help="Height of one row; row labels are sized to the row so they stay "
                                            "readable.")
        fw_in = f0[2].number_input("Figure width (in)", 3.0, 40.0, 8.0, 0.5, key="hm_fig_w",
                                   disabled=sizing != "custom",
                                   help="Used for the heatmap on screen (96 px per inch) and for the download box.")
        fh_in = f0[3].number_input("Figure height (in)", 3.0, 200.0, 8.0, 0.5, key="hm_fig_h",
                                   disabled=sizing != "custom")
        if sizing == "balanced":
            st.caption(f"Balanced cells for the {'first ' if len(prepared) > 1 else ''}heatmap: "
                       f"{bw:.2f} in wide × {bh:.2f} in high ({len(cols_sel)} samples × {n0} rows).")
        f1 = st.columns(4)
        sec_fs = f1[0].number_input("Pathway names (pt, 0 = auto)", 0.0, 30.0, hm.DEFAULT_FONTS["section"], 0.5,
                                    key="hm_fs_section", disabled=not sep_on,
                                    help="Auto = the metabolite-name size; long names wrap or shrink to fit their "
                                         "block.")
        fam = f1[3].selectbox("Font family", hm.FONT_FAMILIES, key="hm_font_family")
        f2 = st.columns(6)
        fonts = dict(
            family=fam,
            title=f2[0].number_input("Title (pt)", 4.0, 30.0, hm.DEFAULT_FONTS["title"], 0.5, key="hm_fs_title"),
            row=f2[1].number_input("Row labels (pt, 0 = auto-fit)", 0.0, 24.0, hm.DEFAULT_FONTS["row"], 0.5,
                                   key="hm_fs_row"),
            group=f2[2].number_input("Metabolite names (pt, 0 = auto)", 0.0, 30.0, hm.DEFAULT_FONTS["group"], 0.5,
                                     key="hm_fs_group", disabled=not grouped),
            col=f2[3].number_input("Column labels (pt)", 3.0, 24.0, hm.DEFAULT_FONTS["col"], 0.5, key="hm_fs_col"),
            ann=f2[4].number_input("Annotations & legend (pt)", 3.0, 24.0, hm.DEFAULT_FONTS["ann"], 0.5,
                                   key="hm_fs_ann"),
            scale=f2[5].number_input("Colour scale (pt)", 3.0, 24.0, hm.DEFAULT_FONTS["scale"], 0.5,
                                     key="hm_fs_scale"),
            section=sec_fs)
        st.caption("Sizes in points apply to the exported figure; the on-screen heatmap uses the same family and "
                   "proportional sizes. Label space grows with the font size, so text never overlaps.")

    vmin = vmax = None
    if not auto:
        vmin, vmax = -lim, lim
    order = An.sort_values(["Group", "Hours"], key=lambda x: x.map({g: i for i, g in enumerate(arms)})
                           if x.name == "Group" else x).Sample
    col_order = [x for x in order if x in cols_sel]
    sel_txt = ("" if (len(gsel) == len(arms) and len(tsel) == len(times))
               else f" — {' vs '.join(str(g) for g in gsel)}, {_hours_txt(tsel)}")
    names = {}
    shown = 0
    for n, (pname, M, msg) in enumerate(prepared):
        sfx = "" if len(prepared) == 1 else f"_{n}"
        if len(prepared) > 1:
            st.markdown(f"#### {pname}")
        if msg:
            if msg.startswith("Only "):
                C.note(msg)
            else:
                st.caption(msg)
        if len(M) < 2:
            st.warning("Fewer than two rows match the selection" + (f" for {pname}." if pname else ".")); continue
        for c_ in {_comp(i) for i in M.index} - set(names):
            names[c_] = C.label_of(c_)
        Z = hm.zscore(M, "row", None).dropna(how="all")
        if not cc_:
            Z = Z[[x for x in col_order if x in Z.columns]]
        ns = (tested - sig_set) if (ns_mark and sig_set is not None) else set()
        row_sections = None
        if sep_on:
            # pathway blocks: rows of one pathway together, ordered (clustered) within the pathway only
            if scope == "pw":
                sec_def = [(p_, set(members[p_])) for p_ in psel]
            else:
                first = C.pathway_of()
                secn = sorted({first.get(_comp(i), "Other / not mapped") for i in Z.index},
                              key=lambda x: (x == "Other / not mapped", x.lower()))
                sec_def = [(p_, {_comp(i) for i in Z.index
                                 if first.get(_comp(i), "Other / not mapped") == p_}) for p_ in secn]
            pairs = []
            for p_, cs_ in sec_def:
                ch_ = [i for i in Z.index if _comp(i) in cs_]
                if not ch_:
                    continue
                Zs = Z.loc[ch_]
                if grouped:
                    o_ = hm.grouped_order(Zs, cr_, metric, method, names)
                elif cr_ and len(Zs) > 2:
                    o_ = [ch_[j] for j in hm.cluster_order(Zs.fillna(0), metric, method)[0]]
                else:
                    o_ = ch_
                pairs += [(p_, c_) for c_ in o_]
            rid = [f"{p_} ‖ {c_}" for p_, c_ in pairs]
            chan_of = dict(zip(rid, [c_ for _, c_ in pairs]))
            row_sections = dict(zip(rid, [p_ for p_, _ in pairs]))
            Z = Z.loc[[c_ for _, c_ in pairs]].set_axis(rid, axis=0)
        else:
            if grouped:
                Z = Z.loc[hm.grouped_order(Z, cr_, metric, method, names)]
            chan_of = {i: i for i in Z.index}
        if grouped:
            ct = hm.channel_text(list(dict.fromkeys(chan_of.values())), ns)
            row_text = {i: ct[chan_of[i]] for i in Z.index}
            row_groups = {i: names.get(_comp(chan_of[i]), _comp(chan_of[i])) for i in Z.index}
        else:
            row_text = {i: chan_of[i] + (" (n.s.)" if chan_of[i] in ns else "") for i in Z.index}
            row_groups = None
        hl = [i for i in Z.index if chan_of[i] in pick]
        ci = ci0.loc[Z.columns]
        cr_eff = cr_ and not grouped and not sep_on
        Zf = Z.fillna(0)
        ro = hm.cluster_order(Zf, metric, method)[0] if cr_eff and len(Z) > 2 else list(range(len(Z)))
        co = hm.cluster_order(Zf.T, metric, method)[0] if cc_ and Z.shape[1] > 2 else list(range(Z.shape[1]))
        ZZ = Z.iloc[ro, co]
        title = (f"Isotopologue channels — {info['unit']}, row Z-score" + (f" — {pname}" if pname else "")
                 + sel_txt)
        # size: balanced (from rows × samples), fixed cells or custom figure
        if sizing == "balanced":
            cw, ch = hm.balanced_cells(len(ZZ), ZZ.shape[1])
        else:
            cw, ch = cell_w_in, cell_h_in
        lab_len = max((len(str(x)) for x in (row_text or {i: i for i in ZZ.index}).values()), default=10)
        grp_len = max((len(v) for v in (row_groups or {}).values()), default=0)
        n_sec = len(set(row_sections.values())) if row_sections else 0
        auto_w = round(0.9 + cw * ZZ.shape[1] + 0.052 * lab_len + (0.07 * grp_len + 0.3 if grp_len else 0)
                       + (0.25 if n_sec else 0) + 1.65, 1)
        auto_h = round(max(3.0, ch * (len(ZZ) + pw_gap * max(n_sec - 1, 0)) + 3.0), 1)
        custom = sizing == "custom"
        fig_w, fig_h = (fw_in, fh_in) if custom else (auto_w, auto_h)
        # keep the download box in step with the size chosen here (only when this size changes)
        if s.get(f"_hm_size_prev{sfx}") != (fig_w, fig_h):
            st.session_state[f"exp_hm{sfx}_w"], st.session_state[f"exp_hm{sfx}_h"] = float(fig_w), float(fig_h)
            st.session_state[f"_hm_size_prev{sfx}"] = (fig_w, fig_h)
        _show_sized(_plotly(ZZ, ci.iloc[co], tracks, True, vmin, vmax, "Row Z-score", cmap, row_labels, hl,
                            title, arms, fonts, fig_w, fig_h, cell_px=None if custom else cw * PX_PER_IN,
                            cell_h_px=None if custom else ch * PX_PER_IN, row_text=row_text,
                            row_groups=row_groups, separators=seps, row_sections=row_sections,
                            section_gap=pw_gap), key="hm_plotly" + sfx)
        shown += 1
        st.caption(f"{len(ZZ)} rows × {ZZ.shape[1]} samples ({' · '.join(str(g) for g in gsel)}; "
                   f"{_hours_txt(tsel)}). Row-wise Z-score of {info['unit']} across the samples shown"
                   + (f"; metabolites ordered by clustering their mean channel profile ({hm.DISTANCES[metric]}, "
                      f"{hm.LINKAGES[method]})" if grouped and cr_ else
                      f"; rows clustered with {hm.DISTANCES[metric]} distance and {hm.LINKAGES[method]} linkage"
                      if cr_eff else "")
                   + (f"; {n_sec} pathway block(s), rows ordered within each pathway only" if n_sec else "")
                   + (f"; (n.s.) = not {thr_txt} in any selected comparison of “{rs}”" if ns else "") + ".")

        def make(w=11.0, h=None, Z=Z, ci=ci, title=title, cw=cw, ch=ch, fig_w=fig_w, fig_h=fig_h, custom=custom,
                 row_text=row_text, row_groups=row_groups, cr_eff=cr_eff, row_sections=row_sections, hl=hl):
            # size unchanged from the automatic one → exact cell size; otherwise stretch to the typed size
            sized = custom or abs(w - fig_w) > 0.05 or abs((h or fig_h) - fig_h) > 0.05
            kw = dict(width=w, height=h) if sized else dict(width=w, height=None, cell_w=cw, cell_h_in=ch)
            f, _, _ = hm.draw_heatmap(Z.fillna(0) if (cr_eff or cc_) else Z, ci, title=title,
                                      value_label="Row Z-score", diverging=True, vmin=vmin, vmax=vmax,
                                      cluster_rows=cr_eff, cluster_cols=cc_, metric=metric, method=method,
                                      row_labels=row_labels, **kw, group_order=arms, highlight=hl,
                                      annotations=tracks, cmap=cmap, fonts=fonts, row_text=row_text,
                                      row_groups=row_groups, separators=seps, row_sections=row_sections,
                                      section_gap=pw_gap)
            return f
        stem = "heatmap_isotopologue_channels" + (f"_{pname}" if pname and len(prepared) > 1 else "")
        chans = [chan_of[i] for i in ZZ.index]
        Tz = ZZ.set_axis(chans, axis=0).reset_index().rename(columns={"index": "Channel"})
        Tz.insert(1, "Metabolite", [names.get(_comp(c_), _comp(c_)) for c_ in chans])
        Tv = M.loc[chans].reset_index().rename(columns={"index": "Channel"})
        if row_sections:
            Tz.insert(0, "Pathway", [row_sections[i] for i in ZZ.index])
            Tv.insert(0, "Pathway", [row_sections[i] for i in ZZ.index])
        C.export_widget(lambda w, h, make=make: make(w, h), stem, "exp_hm" + sfx,
                        (float(fig_w), float(fig_h)), caption=f"Heatmap — {title}",
                        tables={"Z-score matrix (as plotted)": Tz, f"{info['unit']} values of these rows": Tv})
    if shown:
        s.results["heatmap_done"] = True
    C.next_button("kinetics", "Label arrival")


PX_PER_IN = 96
PT = 1.4          # on-screen px per exported pt (keeps the screen and the export in proportion)


def _show_sized(fig, key="hm_plotly"):
    """
    Show the heatmap at its natural pixel width, shrinking to the page when the window is narrower (never cut off):
    the figure is autosized and the element width capped at the natural width, so a narrow window gives narrower
    cells while the label column, legend and colour scale (pixel-anchored) keep their size. Exports are unaffected.
    """
    cfg = {"displaylogo": False, "toImageButtonOptions": {"format": "png", "scale": 5}, "scrollZoom": True,
           "responsive": True}
    w, h = int(fig.layout.width), int(fig.layout.height)
    fig.update_layout(plot_bgcolor="white", paper_bgcolor="white", width=None, autosize=True)
    for kw in (dict(width=w, height=h), dict(use_container_width=True)):
        try:
            return st.plotly_chart(fig, key=key, config=cfg, theme=None, **kw)
        except Exception:
            continue
    return st.plotly_chart(fig, key=key, config=cfg, theme=None)


def _plotly(ZZ, ci, tracks, diverging, vmin, vmax, label, cmap, row_labels, pick, title, arms, fonts=None,
            width_in=11.0, height_in=None, cell_px=None, cell_h_px=None, row_text=None, row_groups=None,
            separators=True, row_sections=None, section_gap=0.5):
    fo = {**hm.DEFAULT_FONTS, **(fonts or {})}
    fam = f"{fo['family']}, Liberation Sans, Helvetica, sans-serif"
    nt = len(tracks)
    nrows = len(ZZ)
    if cell_px:   # fixed cell size; a fixed row-label font keeps its size and makes the rows taller if needed
        ch = cell_h_px or cell_px
        row_px = (fo["row"] * PT) if fo["row"] else max(7.5, min(10.0, ch * 0.85))
        cell_h_px = max(ch, row_px * 1.12)
    else:
        row_px = (fo["row"] * PT) if fo["row"] else 10
    col_px, ann_px, sc_px, ti_px = fo["col"] * PT, fo["ann"] * PT, fo["scale"] * PT, fo["title"] * PT
    disp = [str((row_text or {}).get(r, r)) for r in ZZ.index]
    secs = [str(row_sections.get(r, "")) for r in ZZ.index] if row_sections else None
    pos = hm.row_positions(secs, section_gap) if secs else list(range(nrows))
    span = (pos[-1] + 1) if nrows else 1
    runs = hm.keyed_runs([str((row_groups or {}).get(r, "")) for r in ZZ.index], secs) if row_groups else []
    grp_px = (float(fo.get("group") or 0) * PT) or max(row_px * 1.25, ann_px)
    b_off = 5                                      # small gap between the widest row label and the line
    grp_w = int(hm.text_width([n_ for n_, _, _ in runs], grp_px, hm._families(fo["family"]), bold=True) * 72
                + b_off + 16) if runs else 0
    h_tr = max(16, int(ann_px * 1.5))
    xlab = int(max(len(str(c)) for c in ZZ.columns) * col_px * 0.55 + 12)
    leg_w = int(max([len("Time" if t == "TimeLabel" else t) for t in tracks] +
                    [len(str(v)) for t in tracks for v in ci[t].unique()] + [4]) * ann_px * 0.52 + 30) if tracks else 10
    if row_labels:   # measured label width (font metrics, px) + plotly's tick-label standoff
        rlab = int(hm.text_width(disp, row_px, hm._families(fo["family"]), bold=bool(pick)) * 72 + 3) if runs \
            else int(max(len(d) for d in disp) * row_px * 0.5 + 8)
    else:
        rlab = 6
    right = rlab + grp_w + leg_w + 10
    cb_title_px = int(sc_px * 1.6 + 6)           # colour-scale title (its own line, above the bar)
    cb_px = int(cb_title_px + 14 + sc_px * 1.6 + 16)   # title + bar + tick labels underneath the heatmap
    t_px = int(ti_px * 2.2 + 20)
    if cell_px:
        h_main = int(round(cell_h_px * span))
    elif height_in:
        total = int(height_in * PX_PER_IN)
        h_main = max(total - h_tr * nt - xlab - cb_px - t_px - 20, 120)
    else:
        h_main = max(380, int(max(15, row_px * 1.5) * nrows) + 120)
    # pathway sections (left): vertical line + vertical name, wrapped / reduced to fit the block height
    famL = hm._families(fo["family"])
    sec_px0 = (float(fo.get("section") or 0) * PT) or grp_px
    row_hp = h_main / max(span, 1)
    sec_lab = []
    for name, a0, a1 in (hm.row_runs(secs) if secs else []):
        block = ((pos[a1] - pos[a0] + 1) * row_hp - 4) / 72          # px as "pt" units → inches for the measure
        lines, fs_ = hm.wrap_to_length(name, sec_px0, famL, block)
        sec_lab.append((lines, fs_, a0, a1))
    S_LINE, S_TXT = 7, 6
    sec_w = max((S_LINE + S_TXT + int(len(l_) * f_ * 1.25) + 4 for l_, f_, _, _ in sec_lab), default=0)
    trk_w = int(max([len("Time" if t == "TimeLabel" else t) for t in tracks] + [0]) * ann_px * 0.6 + 12)
    base_l = max(60, sec_w + 4, trk_w)
    width_px = int(width_in * PX_PER_IN)
    if cell_px:   # plot area = exactly n_cols cells; figure grows to fit the title and the colour scale
        plot_w = int(round(cell_px * len(ZZ.columns)))
        width_px = max(base_l + plot_w + right, int(len(title) * ti_px * 0.55) + 40, base_l + 140 + right)
    fig = make_subplots(rows=nt + 1, cols=1, shared_xaxes=True, vertical_spacing=0.004,
                        row_heights=[h_tr] * nt + [h_main])
    cols = list(ZZ.columns)
    acol = hm.annotation_colors(ci, tracks, arms)
    for k, t in enumerate(tracks):
        vals = list(ci[t])
        cats = list(acol[t])
        code = [cats.index(v) if v in cats else np.nan for v in vals]
        n = max(len(cats), 1)
        scale = []
        for i, cv in enumerate(cats):
            scale += [[i / n, acol[t][cv]], [(i + 1) / n, acol[t][cv]]]
        fig.add_trace(go.Heatmap(z=[code], x=cols, y=["Time" if t == "TimeLabel" else t], zmin=-0.5, zmax=n - 0.5,
                                 colorscale=scale or [[0, "#eee"], [1, "#eee"]], showscale=False,
                                 customdata=[[str(v) for v in vals]],
                                 hovertemplate="%{x}<br>" + t + ": %{customdata}<extra></extra>"), row=k + 1, col=1)
    if cmap in (None, "Omnix diverging"):
        cs = pl.plotly_div_scale() if diverging else pl.plotly_seq_scale()
    else:
        cs = hm.PLOTLY_SCALE.get(cmap, cmap)
    lim_ = float(np.nanmax(np.abs(ZZ.values))) if vmax is None else vmax
    dom = fig.layout[f"yaxis{nt + 1}" if nt else "yaxis"].domain or (0, 1)
    dom_mid, dom_len = (dom[0] + dom[1]) / 2, dom[1] - dom[0]
    # numeric row positions (gaps between pathway blocks); one trace per block so cells never span a gap
    blocks = hm.row_runs(secs) if secs else [("", 0, nrows - 1)]
    hov = np.array([[str(r)] * len(cols) for r in ZZ.index], dtype=object)
    for b_, (_, a0, a1) in enumerate(blocks):
        fig.add_trace(go.Heatmap(z=ZZ.values[a0:a1 + 1], x=cols, y=pos[a0:a1 + 1], colorscale=cs,
                                 zmid=0 if diverging else None, showscale=b_ == 0,
                                 zmin=-lim_ if diverging else vmin, zmax=lim_ if diverging else vmax,
                                 customdata=hov[a0:a1 + 1],
                                 colorbar=dict(orientation="h", thickness=14, len=120, lenmode="pixels",
                                               x=1, xanchor="right", xpad=0, yanchor="top",
                                               y=-(xlab + cb_title_px) / (h_main + h_tr * nt + 20), ypad=0,
                                               tickfont=dict(size=max(sc_px - 1, 6), color="black"),
                                               nticks=5, outlinecolor="black", outlinewidth=0.5),
                                 hovertemplate="%{customdata}<br>%{x}<br>%{z:.3g}<extra></extra>"), row=nt + 1, col=1)
    fig.update_yaxes(range=[span - 0.5, -0.5], showticklabels=row_labels, tickfont=dict(size=row_px, color="black"),
                     showgrid=False, zeroline=False, tickmode="array", tickvals=pos,
                     ticktext=[f"<b><span style='color:#a8402c'>{d}</span></b>" if (runs and r_ in pick) else d
                               for r_, d in zip(ZZ.index, disp)],
                     side="right", row=nt + 1, col=1)
    for k in range(nt):
        fig.update_yaxes(tickfont=dict(size=ann_px, color="black"), row=k + 1, col=1)
    fig.update_xaxes(showgrid=False, zeroline=False)
    fig.update_xaxes(tickangle=-90, tickfont=dict(size=col_px, color="black"), tickmode="array", tickvals=cols,
                     ticktext=[str(c) for c in cols], row=nt + 1, col=1)
    yref = f"y{nt + 1}" if nt else "y"
    for rname, y_ in zip(ZZ.index, pos):
        if rname in pick and not runs:      # grouped labels: highlighted rows are bold red (as in the export)
            fig.add_annotation(x=1.0, xref="paper", y=y_, yref=yref, text="◀",
                               showarrow=False, xanchor="left", xshift=rlab, font=dict(color="#a8402c"))
    if cell_px:   # extra width (long title): split evenly left / right so the heatmap block stays centred
        extra = max(0, width_px - (base_l + plot_w + right))
        l_px, right = base_l + extra // 2, right + extra - extra // 2
    else:
        l_px = base_l
    # grouped rows: white separators, a bracket per metabolite and its name right of the row labels
    shapes, g_ann = [], []
    if runs:
        xb = rlab + b_off                          # px right of the plot edge (stays put when the plot shrinks)
        for k, (name, a0, a1) in enumerate(runs):
            if separators and k > 0 and (not secs or secs[a0] == secs[a0 - 1]):
                shapes.append(dict(type="line", xref="paper", yref=yref, x0=0, x1=1, y0=pos[a0] - 0.5,
                                   y1=pos[a0] - 0.5, line=dict(color="white", width=1.5), layer="above"))
            shapes.append(dict(type="line", xref="paper", xsizemode="pixel", xanchor=1, yref=yref, x0=xb, x1=xb,
                               y0=pos[a0] - 0.38, y1=pos[a1] + 0.38, line=dict(color="black", width=1.6)))
            g_ann.append(dict(text=f"<b>{name}</b>", x=1, xref="paper", xanchor="left", xshift=xb + 6,
                              y=(pos[a0] + pos[a1]) / 2, yref=yref, yanchor="middle", showarrow=False,
                              font=dict(size=grp_px, color="black")))
    # pathway sections: vertical line + rotated name on the left (pixel-anchored, stay put when the plot shrinks)
    for lines, fs_, a0, a1 in sec_lab:
        shapes.append(dict(type="line", xref="paper", xsizemode="pixel", xanchor=0, yref=yref, x0=-S_LINE,
                           x1=-S_LINE, y0=pos[a0] - 0.38, y1=pos[a1] + 0.38, line=dict(color="black", width=1.6)))
        g_ann.append(dict(text="<br>".join(f"<b>{l_}</b>" for l_ in lines), x=0, xref="paper", xanchor="center",
                          xshift=-(S_LINE + S_TXT + len(lines) * fs_ * 1.25 / 2), y=(pos[a0] + pos[a1]) / 2,
                          yref=yref, yanchor="middle", textangle=-90, showarrow=False, align="center",
                          font=dict(size=fs_, color="black")))
    # annotation legend: right side, vertically centred on the heatmap
    leg = "<br><br>".join(f"<b>{'Time' if t == 'TimeLabel' else t}</b><br>" + "<br>".join(
        f"<span style='color:{c}'>■</span> {v:g}" if isinstance(v, float) else f"<span style='color:{c}'>■</span> {v}"
        for v, c in acol[t].items()) for t in tracks)
    cb_title = dict(text=label, x=1, xref="paper", xanchor="right", y=0, yref="paper", yanchor="top",
                    yshift=-(xlab + 2), showarrow=False, font=dict(size=sc_px, color="black"))
    n_leg = sum(1 + len(acol[t]) for t in tracks) + max(len(tracks) - 1, 0)
    extra_b = max(0, int(n_leg * ann_px * 1.3) - (h_main + xlab + cb_px)) if tracks else 0
    fig.update_layout(template="plotly_white", width=width_px,
                      height=h_main + h_tr * nt + xlab + cb_px + t_px + 20 + extra_b,
                      font=dict(family=fam, color="black"),
                      title=dict(text=title, x=0.5, xanchor="center",
                                 font=dict(size=ti_px, family=fam, color="black")),
                      margin=dict(l=l_px, r=right, t=t_px, b=xlab + cb_px + extra_b),
                      shapes=shapes,
                      annotations=list(fig.layout.annotations) + [cb_title] + g_ann
                      + ([dict(text=leg, x=1.0, xref="paper", xanchor="left", xshift=rlab + grp_w + 6,
                               y=dom_mid, yref="paper", yanchor="middle", showarrow=False, align="left",
                               font=dict(size=ann_px, color="black"))] if tracks else []))
    return fig
