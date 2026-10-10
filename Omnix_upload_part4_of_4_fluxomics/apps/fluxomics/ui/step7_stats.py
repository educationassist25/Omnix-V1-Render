"""Step 6 — statistics on isotopologue-channel intensities (log2 or log2 + IQR): two-group t-tests and ANOVA."""
from __future__ import annotations

import numpy as np
import pandas as pd
import streamlit as st

from core import design as dz, stats as stt
from ui import common as C

KINDS = {"ttest": "Two-Group (T-Test)", "anova": "ANOVA"}


def render():
    C.header(6, "Statistical analysis",
             "Tests on every isotopologue channel separately (isotopologues are never summed), on the log2 or "
             "log2 + IQR-normalized intensities (active dataset, sidebar). Every result table lists the p value and the Benjamini–Hochberg FDR "
             "for each channel, sorted by p value.")
    if not C.require("upload", "design", "normalize"):
        return
    s = C.ss()
    D = C.design_final()
    groups, times = dz.groups_of(D), dz.times_of(D)
    lay = "int"                         # active intensity dataset (sidebar switch)
    X, info = C.layer(lay)
    crit, alpha = "P value", None       # v2.3: no threshold — p and FDR are reported for every channel
    st.caption(f"{len(X)} isotopologue channels · {X.shape[1]} samples · effects are log2 fold changes (B − A).")
    kind = st.radio("Analysis", list(KINDS), format_func=KINDS.get, horizontal=True, key="st_kind")
    an = dz.analysed(D)
    settings = dict(layer=lay, kind=kind, n=len(X))

    if kind == "ttest":
        q = st.columns(4)
        method = q[0].selectbox("Test", list(stt.TTESTS), format_func=stt.TTESTS.get, key="st_tt")
        if len(groups) < 2:
            st.warning("Two-group tests need at least two groups."); return
        ga = q[1].selectbox("Group A (reference)", groups, index=groups.index(s.group_a) if s.group_a in groups else 0,
                            key="st_a")
        rest = [g for g in groups if g != ga]
        gb = q[2].selectbox("Group B", rest, index=rest.index(s.group_b) if s.group_b in rest else 0, key="st_b")
        tsel = q[3].multiselect("Time points (each tested separately)", times, default=times,
                                format_func=lambda h: f"{h:g} h", key="st_times")
        C.direction_banner(ga, gb, True)
        if method == "paired":
            st.caption("Paired: samples of A and B are matched by the biological-replicate number (Step 2) within "
                       "each time point.")
        settings.update(method=method, group_a=ga, group_b=gb, times=list(tsel))
    else:
        q = st.columns(3)
        method = q[0].selectbox("ANOVA method", list(stt.ANOVAS), format_func=stt.ANOVAS.get, key="st_am")
        gsel = q[1].multiselect("Groups", groups, default=groups, key="st_groups")
        tsel = q[2].multiselect("Time points", times, default=times, format_func=lambda h: f"{h:g} h",
                                key="st_atimes")
        if method != "twoway":
            r = st.columns(2)
            factor = r[0].radio("Compare", ["groups", "times"], horizontal=True, key="st_factor",
                                format_func={"groups": "Groups, at each time point",
                                             "times": "Time points, within each group"}.get)
            ph = r[1].selectbox("Post hoc test", stt.POSTHOC_BY[method], format_func=stt.POSTHOC.get, key="st_ph")
            st.caption("Post hoc comparisons are run for every tested channel.")
            settings.update(method=method, groups=gsel, times=list(tsel), factor=factor, posthoc=ph)
        else:
            st.caption("Model per channel: y ~ group + time + group:time (Type II sums of squares). The interaction "
                       "asks whether the groups differ in how the channel changes over time.")
            settings.update(method=method, groups=gsel, times=list(tsel))

    if st.button("Run statistics", type="primary", key="st_go"):
        try:
            with st.spinner(f"Testing {len(X)} isotopologue channels…"):
                res = _run(X, D, an, kind, settings, crit, alpha)
        except Exception as e:
            st.error(f"Statistics failed: {e}")
            return
        name = (f"{C.INT_LABEL[s.intensity_choice].split(' intensity')[0]} · " +
                (f"{stt.TTESTS[settings['method']].split(' (')[0]} · {settings['group_b']} vs {settings['group_a']}"
                 if kind == "ttest" else f"{stt.ANOVAS[settings['method']].split(' (')[0]}"
                 + (f" · {settings['factor']}" if settings["method"] != "twoway" else "")))
        res.update(settings=settings, info=info, label=name)
        s.stats[name] = res
        s.stats_current = name
    if not s.stats:
        return
    keys = list(s.stats)
    cur = st.selectbox("Result set", keys, index=keys.index(s.get("stats_current")) if s.get("stats_current") in keys
                       else len(keys) - 1, key="st_cur")
    _show(s.stats[cur], cur, crit, alpha)
    C.next_button("pathways", "Pathway mapping")


def _pairs(an, ga, gb, h):
    a = an[(an.Group == ga) & np.isclose(an.Hours, h)].set_index("Replicate").Sample
    b = an[(an.Group == gb) & np.isclose(an.Hours, h)].set_index("Replicate").Sample
    return [(a[r], b[r]) for r in a.index if r in b.index]


def _annotate(T):
    parts = C.split_channel(T.Feature)
    T.insert(T.columns.get_loc("Feature") + 1, "Compound", parts.Compound.values)
    T.insert(T.columns.get_loc("Compound") + 1, "Isotopologue", parts.Isotopologue.values)
    return T


def _run(X, D, an, kind, p, crit, alpha):
    out, ph = [], []
    if kind == "ttest":
        for h in p["times"]:
            sa, sb = dz.samples_of(D, p["group_a"], h), dz.samples_of(D, p["group_b"], h)
            if not sa or not sb:
                continue
            T = stt.two_group(X, sa, sb, p["method"], True, _pairs(an, p["group_a"], p["group_b"], h))
            T = T.rename(columns={"n A": f"n {p['group_a']}", "n B": f"n {p['group_b']}",
                                  "Mean A": f"Mean {p['group_a']}", "Mean B": f"Mean {p['group_b']}",
                                  "SD A": f"SD {p['group_a']}", "SD B": f"SD {p['group_b']}"})
            T.insert(0, "Comparison", f"{p['group_b']} vs {p['group_a']} · {h:g} h")
            out.append(T)
    elif p["method"] == "twoway":
        sub = an[an.Group.isin(p["groups"]) & an.Hours.isin(p["times"])][["Sample", "Group", "TimeLabel"]]
        T = stt.two_way_anova(X[sub.Sample], sub, "Group", "TimeLabel", "bh")
        T = T.rename(columns={c: c.replace("q ", "FDR ") for c in T.columns if c.startswith("q ")})
        T.insert(0, "Comparison", "Group × Time")
        out.append(T)
    else:
        if p["factor"] == "groups":
            blocks = [(f"Groups at {h:g} h", {g: dz.samples_of(D, g, h) for g in p["groups"]}) for h in p["times"]]
        else:
            blocks = [(f"Time points within {g}", {f"{h:g} h": dz.samples_of(D, g, h) for h in p["times"]})
                      for g in p["groups"]]
        for name, grp in blocks:
            grp = {k: v for k, v in grp.items() if v}
            if len(grp) < 2:
                continue
            T, P = stt.anova_method(X, grp, p["method"], p["posthoc"], True, alpha, crit)
            T.insert(0, "Comparison", name)
            out.append(T)
            if len(P):
                P.insert(0, "Comparison", name)
                ph.append(_annotate(P))
    if not out:
        raise ValueError("No comparison had samples in both / all groups.")
    T = _annotate(pd.concat(out, ignore_index=True))
    return dict(table=T, posthoc=pd.concat(ph, ignore_index=True) if ph else pd.DataFrame())


def _show(res, name, crit, alpha):
    T = res["table"].copy()
    st.markdown(f"#### {name}")
    summ = T.groupby("Comparison", sort=False).apply(
        lambda d: pd.Series({"Channels tested": int(d.filter(regex=r"^p").notna().any(axis=1).sum()),
                             "Metabolites": int(d.Compound.nunique())}))
    C.show_df(summ)
    comp = st.selectbox("Comparison", list(summ.index), key="st_cmp")
    d = T[T.Comparison == comp]
    sortc = "p" if "p" in d else [c for c in d.columns if c.startswith("p ")][-1]
    C.show_df(d.sort_values(sortc), height=420)
    P = res.get("posthoc")
    if P is not None and len(P):
        with st.expander(f"Post hoc comparisons ({len(P)} rows)"):
            C.show_df(P[P.Comparison == comp] if "Comparison" in P else P, height=320)
    # downloads
    sett = pd.DataFrame([(k, ", ".join(map(str, v)) if isinstance(v, (list, tuple)) else v)
                         for k, v in res["settings"].items()], columns=["Setting", "Value"])
    sheets = {"Summary": summ.reset_index(), "Results": T, "Post_hoc": P if P is not None else pd.DataFrame(),
              "Settings": sett, "Design": C.design_final()}
    stem = "Statistics_" + "".join(ch if ch.isalnum() else "_" for ch in name)[:60]
    C.download_panel(stem, f"st_dl_{abs(hash(name)) % 10**8}",
                     tables={"Results (all comparisons)": T, "Post hoc": P if P is not None else pd.DataFrame()},
                     workbook=sheets)
