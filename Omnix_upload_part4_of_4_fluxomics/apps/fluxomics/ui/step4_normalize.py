"""Step 4 — log2 transformation followed by robust IQR normalization; both datasets stay available."""
from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from core import design as dz, normalization as nm
from core.style import group_colors
from ui import common as C


def render():
    C.header(3, "Normalization",
             "log2 transformation, then robust IQR normalization of each labeled sample. Both datasets are kept; the "
             "sidebar switch “Intensity dataset for downstream steps” chooses which one every later intensity "
             "analysis uses (default: log2 + IQR).")
    if not C.require("upload", "design"):
        return
    s = C.ss()
    D = C.design_final()
    An = dz.analysed(D)
    with st.expander("How the normalization works", expanded=False):
        st.markdown(
            "1. **log2:** L = log2(area). Zero / low values in labeled samples were already removed in Step 1; values "
            "are not imputed, so an empty cell stays missing.\n"
            "2. **Robust IQR:** for each sample s, N = (L − median_s) / IQR_s × median(IQR) + median(median_s), with "
            "median and IQR over all isotopologue channels of that sample. Each sample gets the same median and "
            "interquartile range; values stay in log2 units, so group differences remain log2 fold changes.\n"
            "3. **Every later intensity analysis** (bar plots, statistics, heatmaps, pathway abundance) uses the "
            "dataset selected in the sidebar box — log2 + IQR-normalized (default) or log2.")
    excl = D[~D.Sample.isin(An.Sample)]
    ref_bl = excl[excl.SampleType.isin(["12C reference", "Blank"])]
    other = excl[~excl.SampleType.isin(["12C reference", "Blank"])]
    C.note(f"<b>Samples normalized:</b> {len(An)} labeled samples. <b>Excluded from normalization:</b> "
           + (", ".join(f"{r.Sample} ({r.SampleType})" for r in ref_bl.itertuples()) or "none")
           + (f"; switched off in Step 2: {', '.join(other.Sample)}" if len(other) else "")
           + ". 12C reference and blank samples never enter the median / IQR and are not part of the normalized "
             "file used downstream.", ok=True)
    if st.button("Re-run log2 + IQR normalization", key="norm_go") or s.norm is None:
        R = C.raw_matrix()                       # labeled, included samples only (no 12C reference, no blank)
        R = R[[c_ for c_ in An.Sample if c_ in R.columns]]
        L, bad = nm.log2_transform(R)
        N, info = nm.robust_iqr(L, None)         # median / IQR over all isotopologue channels
        s.norm = dict(L=L, N=N, info=info, bad=bad, settings=dict(reference="all channels", factors="none"))
        s.stats, s.results = {}, {}
    Nn = s.norm
    if Nn["bad"]:
        C.note(f"{Nn['bad']} non-positive values could not be log-transformed and are missing.")
    chk = nm.check(Nn["L"], Nn["N"])
    m = st.columns(4)
    m[0].metric("SD of sample medians — log2", f"{chk['SD of sample medians (log2)']:.3f}")
    m[1].metric("… after IQR normalization", f"{chk['SD of sample medians after']:.3f}")
    m[2].metric("SD of sample IQRs — log2", f"{chk['SD of sample IQRs (log2)']:.3f}")
    m[3].metric("… after IQR normalization", f"{chk['SD of sample IQRs after']:.3f}")
    arms = dz.groups_of(D)
    gc = group_colors(arms)
    order = An.sort_values(["Group", "Hours"], key=lambda x: x.map({g: i for i, g in enumerate(arms)})
                           if x.name == "Group" else x).Sample.tolist()
    t1, t2, t3 = st.tabs(["Sample distributions", "Scale factors", "Normalized data"])
    with t1:
        fig = go.Figure()
        for k, (lab, M) in enumerate((("log2", Nn["L"]), ("log2 + IQR", Nn["N"]))):
            for smp in order:
                g = An.set_index("Sample").Group[smp]
                fig.add_trace(go.Box(y=M[smp].dropna(), name=smp, legendgroup=g, marker_color=gc[g], boxpoints=False,
                                     showlegend=False, xaxis=f"x{k + 1}" if k else "x", yaxis=f"y{k + 1}" if k else "y",
                                     line=dict(width=1)))
        fig.update_layout(grid=dict(rows=2, columns=1, pattern="independent"), template="plotly_white", height=640,
                          margin=dict(l=60, r=10, t=40, b=120), title="Per-sample distributions (top: log2, bottom: "
                          "log2 + IQR; colour = group)")
        fig.update_layout(yaxis=dict(title="log2"), yaxis2=dict(title="log2 + IQR"))
        fig.update_xaxes(tickangle=-60, tickfont=dict(size=9))
        C.show_plotly(fig, key="norm_box")
        C.export_widget(lambda w, h: _dist_mpl(Nn, order, An, gc, w, h), "normalization_sample_distributions",
                        "exp_norm", (max(8.0, 0.32 * len(order) + 1.5), 6.5),
                        caption="Per-sample distributions before / after IQR normalization",
                        title="⬇ Downloads — both panels (log2 and log2 + IQR)")
        C.export_widget(lambda w, h: _dist_one_mpl(Nn["L"], "log2 intensity", order, An, gc, w, h),
                        "sample_distributions_log2", "exp_norm_l2", (max(8.0, 0.32 * len(order) + 1.5), 4.0),
                        caption="Per-sample distributions — log2", title="⬇ Downloads — log2 plot")
        C.export_widget(lambda w, h: _dist_one_mpl(Nn["N"], "log2 + IQR-normalized intensity", order, An, gc, w, h),
                        "sample_distributions_log2_IQR", "exp_norm_iqr", (max(8.0, 0.32 * len(order) + 1.5), 4.0),
                        caption="Per-sample distributions — log2 + IQR", title="⬇ Downloads — log2 + IQR plot")
    with t2:
        C.show_df(Nn["info"].round(3))
        C.table_downloads(Nn["info"].reset_index(), "IQR_normalization_scale_factors", "norm_sf",
                          "Scale factors")
    with t3:
        st.caption(f"{Nn['N'].shape[0]} isotopologue channels × {Nn['N'].shape[1]} labeled samples — this "
                   "log2 + IQR-normalized table is the default dataset for later intensity analyses (switch in the sidebar).")
        C.show_df(Nn["N"].round(3), height=380)
        C.download_panel("Fluxomics_normalized", "norm_dl",
                         tables={"log2 + IQR normalized (used downstream)": _annotated(Nn["N"], s.ds),
                                 "log2 (before IQR)": _annotated(Nn["L"], s.ds)},
                         workbook={"log2_IQR_normalized": _annotated(Nn["N"], s.ds), "log2": _annotated(Nn["L"], s.ds),
                                   "Normalized_samples": An[["Sample", "SampleType", "Group", "TimeLabel",
                                                              "Replicate"]],
                                   "Excluded_samples": excl[["Sample", "SampleType", "Group"]].assign(
                                       Reason=np.where(excl.SampleType.isin(["12C reference", "Blank"]),
                                                       "12C reference / blank — not normalized",
                                                       "switched off in Step 2")),
                                   "Scale_factors": Nn["info"].reset_index()})
    C.next_button("incorporation", "Isotope incorporation")


def _dist_mpl(Nn, order, An, gc, w, h):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 1, figsize=(w, h), sharex=True)
    grp = An.set_index("Sample").Group
    for ax, (lab, M) in zip(axes, (("log2", Nn["L"]), ("log2 + IQR normalized", Nn["N"]))):
        data = [M[smp].dropna().to_numpy() for smp in order]
        bp = ax.boxplot(data, patch_artist=True, showfliers=False, widths=0.6,
                        medianprops=dict(color="#1d1d22", lw=1))
        for patch, smp in zip(bp["boxes"], order):
            patch.set(facecolor=gc[grp[smp]], alpha=0.55, edgecolor=gc[grp[smp]], lw=0.8)
        ax.set_ylabel(lab, fontsize=9)
        ax.grid(axis="y", color="#e7e6ec", lw=0.5)
    axes[1].set_xticks(range(1, len(order) + 1))
    axes[1].set_xticklabels(order, rotation=60, ha="right", fontsize=7)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c, alpha=0.6) for c in gc.values()]
    axes[0].legend(handles, list(gc), frameon=False, fontsize=8, ncol=len(gc), loc="lower left",
                   bbox_to_anchor=(0, 1.02))
    fig.tight_layout()
    return fig


def _dist_one_mpl(M, label, order, An, gc, w, h):
    """One panel: per-sample box plots of one dataset (log2 or log2 + IQR)."""
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(w, h))
    grp = An.set_index("Sample").Group
    data = [M[smp].dropna().to_numpy() for smp in order]
    bp_ = ax.boxplot(data, patch_artist=True, showfliers=False, widths=0.6, medianprops=dict(color="#000000", lw=1))
    for patch, smp in zip(bp_["boxes"], order):
        patch.set(facecolor=gc[grp[smp]], alpha=0.55, edgecolor=gc[grp[smp]], lw=0.8)
    ax.set_ylabel(label, fontsize=9)
    ax.grid(axis="y", color="#e7e6ec", lw=0.5)
    ax.set_xticks(range(1, len(order) + 1))
    ax.set_xticklabels(order, rotation=60, ha="right", fontsize=7)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c, alpha=0.6) for c in gc.values()]
    ax.legend(handles, list(gc), frameon=False, fontsize=8, ncol=len(gc), loc="lower left", bbox_to_anchor=(0, 1.02))
    ax.set_title(f"Per-sample distributions — {label}", loc="left", fontsize=10, pad=18)
    fig.tight_layout()
    return fig


def _annotated(M, ds):
    """Channel matrix with Compound / Formula / Method / Isotopologue columns in front."""
    meta = ds.CT.set_index("Channel")[["Compound", "Formula", "Method", "i"]].reindex(M.index)
    out = meta.assign(Isotopologue="M+" + meta.i.astype("Int64").astype(str)).drop(columns="i")
    return pd.concat([out, M], axis=1).reset_index(names="Channel")
