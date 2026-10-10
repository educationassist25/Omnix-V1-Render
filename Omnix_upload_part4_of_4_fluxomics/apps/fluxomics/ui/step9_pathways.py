"""Step 7 — pathway mapping: one hand-drawn pathway at a time; every detected metabolite carries a Z-score heatmap
(Z-score from the same log2 / log2 + IQR data as the Heatmap step). Metabolite IDs of the shown pathway are verified
online (KEGG / PubChem) automatically."""
from __future__ import annotations

import numpy as np
import pandas as pd
import streamlit as st

from core import design as dz, pathways as pw, pathmap as pmap, pathway_templates as ptpl
from ui import common as C

CONF_HELP = ("**High** exact name + formula · **Medium** same compound ignoring D/L (not resolved by LC-MS) · "
             "**Low** isomer / probable mis-annotation (not assigned unless you accept it) · **Class rule** compound "
             "class only · **Unmapped** not in the library · **User-confirmed** edited, accepted or supplied by you.")


def render():
    C.header(7, "Pathway mapping",
             "One pathway at a time, drawn as a curated pathway map. Every detected metabolite is a card with a "
             "Z-score heatmap of its isotopologue channels (row-wise Z-score from the intensity dataset in the "
             "sidebar, the same data as the Heatmap step). Rows = M+0 … M+n, columns = groups × time points. "
             "Metabolites in italics were not detected. Identifiers of the shown pathway are verified online "
             "automatically.")
    if not C.require("upload", "design", "normalize"):
        return
    s = C.ss()
    A = C.annotation()
    sel = selection(s, A)
    if sel is not None:
        A = auto_verify(s, A, sel)
    if sel is not None:
        map_tab(s, A, sel)
    with st.expander("Metabolite annotation (identifiers, pathways, confidence) — edit if needed", expanded=False):
        annotation_tab(s, A)
    C.next_button("heatmap", "Heatmaps")


# ---------------------------------------------------------------------------
def selection(s, A):
    """One pathway, groups, time points, panels and colour settings."""
    D = C.design_final()
    groups, times = dz.groups_of(D), dz.times_of(D)
    F, S = C.labeling()
    bad = set(C.interference(S, D))
    detected = set(S.Compound) - bad
    opts = ptpl.options()
    n_det = {k: len(set(pmap.bind(ptpl.T[k], A, detected).values())) for k, _ in opts}
    keys = [k for k, _ in opts]
    with st.container(border=True) if _has_border() else st.container():
        c = st.columns([3, 2, 2])
        key = c[0].selectbox("Pathway (one at a time)", keys, index=keys.index("glycolysis"), key="pm_tpl",
                             format_func=lambda k: f"{ptpl.T[k].title}  ·  {n_det[k]} detected",
                             help="Each pathway is shown on its own, as a curated map with its reactions.")
        gsel = c[1].multiselect("Groups", groups, default=groups, key="pm_groups")
        tsel = c[2].multiselect("Time points", times, default=times, format_func=lambda h: f"{h:g} h",
                                key="pm_times")
        d = st.columns(5)
        m0 = d[0].toggle("Include M+0 row", True, key="pm_m0")
        ext = d[1].toggle("Other members below", True, key="pm_extra",
                          help="Detected metabolites of the matching annotation pathway(s) that are not on this "
                               "drawing.")
        vic = d[2].toggle("Print values in cells", False, key="pm_vals")
        z_cmap = d[3].selectbox("Z-score colours", list(pmap.DIV_CMAPS), key="pm_z_cmap")
        zlim = d[4].number_input("Z-score colour limit ±", 0.5, 10.0, 2.0, 0.5, key="pm_zlim")
    gsel = [g for g in groups if g in gsel]
    cols = [(g, h) for g in gsel for h in sorted(tsel)]
    return dict(key=key, groups=gsel, times=sorted(tsel), cols=cols, panels=["z"], m0=m0, extra=ext, vic=vic,
                z_cmap=z_cmap, zlim=zlim, bad=bad, detected=detected)


def auto_verify(s, A, sel):
    """Verify the identifiers of the selected pathway's metabolites online (KEGG / PubChem), once per compound and
    session. Runs automatically; quietly skipped when offline (retry button provided)."""
    tpl = ptpl.T[sel["key"]]
    binding = pmap.bind(tpl, A, sel["detected"])
    comps = set(binding.values()) | set(pmap.extras(tpl, A, sel["detected"], binding))
    tried = st.session_state.setdefault("_pm_verify_tried", {}).setdefault(id(s.ds), set())   # per dataset
    todo = sorted(c for c in comps if c not in tried and c in set(A.loc[A.LibraryID.astype(str) != "", "Compound"]))
    box = st.empty()
    if st.session_state.get("_pm_offline"):
        with box.container():
            c = st.columns([4, 1])
            c[0].caption("Online ID verification skipped — KEGG / PubChem not reachable from this computer.")
            if c[1].button("Retry verification", key="pm_ver_retry"):
                st.session_state["_pm_offline"] = False
                st.rerun()
        return A
    if not todo:
        v = A[A.Compound.isin(comps)]
        if "Verification" in v and len(v):
            ok = int(v.Verification.astype(str).eq("verified").sum())
            flag = int((v.Verification.astype(str).str.len() > 0).sum()) - ok
            box.caption(f"Online verification (KEGG / PubChem): {ok} verified"
                        + (f", {flag} with notes — see the Verification column in the annotation table below"
                           if flag else "") + ".")
        return A
    if not pw.online_available():
        st.session_state["_pm_offline"] = True
        box.caption("Online ID verification skipped — KEGG / PubChem not reachable from this computer.")
        return A
    bar = box.progress(0.0, f"Verifying {len(todo)} metabolite(s) of this pathway online…")
    try:
        A2 = pw.verify_online(A, compounds=todo,
                              progress=lambda f, name: bar.progress(f, f"Verifying online: {name}"))
        tried |= set(todo)
        s.annotation = A2
        box.empty()
        return A2
    except Exception as e:
        box.caption(f"Online verification failed: {e}")
        tried |= set(todo)
        return A


def _has_border():
    import inspect
    try:
        return "border" in inspect.signature(st.container).parameters
    except Exception:
        return False


def map_tab(s, A, sel):
    D = C.design_final()
    F, S = C.labeling()
    if not sel["cols"]:
        st.info("Choose at least one group and one time point."); return
    tpl = ptpl.T[sel["key"]]
    binding = pmap.bind(tpl, A, sel["detected"])
    extra = pmap.extras(tpl, A, sel["detected"], binding) if sel["extra"] else []
    comps = set(binding.values()) | set(extra)
    # Z-score: row-wise Z of each channel over all analysed samples, from the active intensity dataset
    X, info = (C.layer("int") if s.norm is not None else (None, None))
    z_b = pmap.zscore_blocks(X, D, sel["cols"], comps) if (X is not None and comps) else {}
    z_label = f"Z-score ({info['unit']})" if info else "Z-score"
    if X is None:
        C.note("Z-scores need the normalized data — run <b>Step 3 · Normalization</b> first.")
        return
    blocks = {"z": z_b}
    rows = {}
    for c_ in comps:
        idx = set()
        for p_ in sel["panels"]:
            if c_ in blocks[p_]:
                idx |= set(int(i) for i in blocks[p_][c_].index)
        if not sel["m0"]:
            idx.discard(0)
        if idx:
            rows[c_] = sorted(idx)
    labels = {c_: C.label_of(c_) for c_ in comps}
    pm_ = pmap.layout(tpl, binding, rows, sel["cols"], sel["panels"], labels=labels, extra=extra,
                      values_in_cells=sel["vic"])
    norm_txt = C.INT_LABEL[s.intensity_choice] if s.norm is not None else "—"
    subtitle = f"Z-score = row-wise Z of each isotopologue channel ({norm_txt}), mean per group × time"
    vh = st.select_slider("Viewer height (px)", [500, 650, 800, 950, 1100, 1300], 800, key="pm_vh",
                          help="Drag to move, mouse wheel or + / − to zoom, Fit = whole width. The map is a vector "
                               "drawing: names, values and labels scale with the zoom. Hover a heatmap cell for its "
                               "value.")
    svg = pmap.render_svg(pm_, blocks, rows, z_lim=sel["zlim"], z_cmap=sel["z_cmap"], subtitle=subtitle, values_in_cells=sel["vic"], z_label=z_label)
    _show_html(pmap.viewer_html(svg, vh), vh + 50)
    s.results["maps_done"] = True
    on_map = int((pm_.nodes.kind == "card").sum()) - len([c_ for c_ in extra if c_ in rows])
    n_met = sum(1 for n in tpl.nodes.values() if n["kind"] == "met")
    st.caption(f"{tpl.title}: {on_map} of {n_met} metabolites on the map were detected"
               + (f"; {len([c_ for c_ in extra if c_ in rows])} other detected member(s) listed below the map"
                  if extra else "")
               + f". Columns: {pmap.column_key(sel['cols'])}. Intensity dataset for Z-scores: {norm_txt} (sidebar). "
               "Metabolites are matched by KEGG ID (name as fallback); italic = not detected; dashed arrows = "
               "several steps or links to other pathways.")
    T = pmap.values_table(blocks, rows, sel["cols"], labels)
    tables = {}
    if len(T):
        tables[z_label] = T[["Metabolite", "Compound", "Isotopologue"]
                                + [c_ for c_ in T.columns if c_.endswith("Z-score")]]
    stem = "pathway_" + tpl.key
    x0, x1, y0, y1 = pm_.extent
    n_cb = len(sel["panels"])
    C.export_widget(lambda w, h: pmap.render_mpl(pm_, blocks, rows, z_lim=sel["zlim"],
                                                 z_cmap=sel["z_cmap"], subtitle=subtitle,
                                                 width=w, values_in_cells=sel["vic"], z_label=z_label),
                    stem, "exp_pm",
                    (round((x1 - x0) * 6.5 / 72 + 0.45 + 1.15 * n_cb, 1),
                     round((y1 - y0) * 6.5 / 72 + 0.8, 1)),
                    caption=f"Pathway map — {tpl.title}", tables=tables or None,
                    workbook=None)
    st.caption("Figure height follows the width automatically (fonts scale with the map), so the layout stays free "
               "of overlaps at any size.")
    node_tab = pd.DataFrame([dict(Node=n["label"], KEGG=n["kegg"], Compound=binding.get(nid, ""))
                             for nid, n in tpl.nodes.items() if n["kind"] == "met"])
    with st.popover("Map node ↔ metabolite table"):
        C.show_df(node_tab, height=320)


def _show_html(page: str, height: int):
    """Embed the interactive viewer (works on old and new Streamlit)."""
    if hasattr(st, "iframe"):                       # Streamlit ≥ 1.5x
        try:
            st.iframe(page, height=height)
            return
        except Exception:
            pass
    import streamlit.components.v1 as components
    components.html(page, height=height, scrolling=False)


def annotation_tab(s, A):
    st.markdown(CONF_HELP)
    cnt = A.Confidence.value_counts()
    m = st.columns(6)
    for i, k in enumerate(["High", "Medium", "User-confirmed", "Low", "Class rule", "Unmapped"]):
        m[i].metric(k, int(cnt.get(k, 0)))
    flt = st.multiselect("Show confidence", ["High", "Medium", "User-confirmed", "Low", "Class rule", "Unmapped"],
                         default=["High", "Medium", "User-confirmed", "Low", "Class rule", "Unmapped"], key="ann_flt")
    cols = [c for c in ["Compound", "Formula", "Method", "Confidence", "MappingStatus", "StandardName", "Candidate",
                        "KEGG", "HMDB", "PubChemCID", "ChEBI", "BiGG", "Pathways", "CandidatePathways", "Note",
                        "IDSource", "Verification"] if c in A.columns]
    view = A[A.Confidence.isin(flt)][cols]
    st.caption("Edit identifiers or pathways directly (use “;” between pathways). Edited rows become User-confirmed.")
    ed = st.data_editor(view, key="ann_editor", hide_index=True, height=420,
                        disabled=[c for c in ["Compound", "Formula", "Method", "Confidence", "MappingStatus",
                                              "Candidate", "CandidatePathways", "IDSource", "Verification"] if c in cols])
    c1, c2, c3 = st.columns(3)
    if c1.button("Save edits", key="ann_save"):
        base = A.set_index("Compound")
        chg = [r for _, r in ed.iterrows() if any(str(r[c]) != str(base.loc[r.Compound, c]) for c in
                                                   ["StandardName", "KEGG", "HMDB", "PubChemCID", "ChEBI", "Pathways"])]
        if chg:
            s.annotation = pw.apply_overrides(A, pd.DataFrame(chg))
            st.success(f"{len(chg)} row(s) saved as User-confirmed.")
            st.rerun()
        else:
            st.info("No changes.")
    cand = A[A.Confidence.isin(["Low", "Class rule"]) & (A.CandidatePathways != "")].Compound.tolist()
    with c2.popover("Accept suggestions…"):
        acc = st.multiselect("Low-confidence / class-rule rows to accept", cand, key="ann_acc")
        if st.button("Accept selected", key="ann_acc_btn") and acc:
            s.annotation = pw.accept_suggestions(A, acc)
            st.rerun()
    with c3.popover("Verify / fill IDs online"):
        st.caption("Queries KEGG REST and PubChem from this computer (internet needed). Fills blank IDs and flags "
                   "mismatches; results are cached in data/id_cache.json.")
        if st.button("Run online verification", key="ann_ver"):
            bar = st.progress(0.0, "Starting…")
            try:
                s.annotation = pw.verify_online(A, progress=lambda f, name: bar.progress(f, f"Checking {name}"))
                st.success("Verification finished — see the Verification column.")
                st.rerun()
            except Exception as e:
                st.error(f"Online verification failed: {e}")
    C.table_downloads(A, "metabolite_annotation", "ann_tab")
    up = st.file_uploader("Restore an edited annotation table (CSV)", type=["csv"], key="ann_up")
    if up is not None and st.button("Apply uploaded annotations", key="ann_apply"):
        O = pd.read_csv(up, dtype=str).fillna("")
        s.annotation = pw.apply_overrides(A, O[O.Confidence.eq("User-confirmed")] if "Confidence" in O else O)
        st.rerun()


