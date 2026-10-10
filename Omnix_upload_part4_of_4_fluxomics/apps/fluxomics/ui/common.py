"""Shared UI helpers: session state, derived data layers, cached computations, export widgets."""
from __future__ import annotations

import hashlib
import io
import zipfile

import numpy as np
import pandas as pd
import streamlit as st

from core import design as dz, labeling as lb, pathways as pw
from core.io_reorganize import process_files, TYPE_LAB
from core.style import EXPORT_FORMATS, MIME, fig_bytes

STEPS = [
    ("upload", "1 · Upload & inspect"),
    ("design", "2 · Experimental design"),
    ("normalize", "3 · Normalization"),
    ("incorporation", "4 · Isotope incorporation"),
    ("boxplots", "5 · Dynamic bar plots"),
    ("stats", "6 · Statistics"),
    ("pathways", "7 · Pathway mapping"),
    ("heatmap", "8 · Heatmaps"),
    ("kinetics", "9 · Label arrival"),
    ("export", "10 · Export & report"),
]
N_STEPS = len(STEPS)
INT_LABEL = {"log2": "log2 intensity", "log2_iqr": "log2 + IQR-normalized intensity"}
INT_LAYERS = {"int": "Intensity — active dataset (sidebar)"}
LAYERS = {
    "int": "Intensity — active dataset (switch in the sidebar)",
    "int_log2": "Intensity — log2",
    "int_iqr": "Intensity — log2 + IQR normalized",
    "fraction": "Isotopologue fraction (%)",
    "incorporation": "Incorporation (%)",
    "enrichment": "Fractional enrichment (%)",
    "mean_labels": "Mean label number",
}
METRIC_COL = {"incorporation": "Incorporation", "enrichment": "Enrichment", "mean_labels": "MeanLabels"}
METRIC_NAME = {"incorporation": "Incorporation (%)", "enrichment": "Fractional enrichment (%)",
               "mean_labels": "Mean label number"}


def ss():
    return st.session_state


def init_state():
    defaults = dict(step="upload", ds=None, files_meta=[], design=None, design_confirmed=False,
                    tracer=dz.Tracer().to_dict(), correction=dz.CorrectionSettings().to_dict(),
                    norm_settings=dict(reference="all", factors={}), norm=None, intensity_choice="log2_iqr",
                    annotation=None, labeling=None, labeling_sig=None, results={}, stats={}, fits=None,
                    fits_metric=None, interpretation="", ai_text="", report_figs=[], prefer={}, keep_partial=False,
                    group_a=None, group_b=None, flux=None, min_area=1000.0)
    for k, v in defaults.items():
        if k not in st.session_state:
            st.session_state[k] = v


def reset_downstream(level: str = "design"):
    """Clear results that depend on an earlier step."""
    s = ss()
    if level in ("upload",):
        s.design_confirmed = False
        s.annotation = None
    if level in ("upload", "design"):
        s.norm = None
    s.labeling = None
    s.results, s.stats, s.fits, s.interpretation, s.ai_text = {}, {}, None, "", ""


def done(step: str) -> bool:
    s = ss()
    return {"upload": s.ds is not None, "design": s.design_confirmed,
            "normalize": s.norm is not None, "incorporation": s.labeling is not None,
            "boxplots": bool(s.results.get("box_done")), "stats": bool(s.stats),
            "heatmap": bool(s.results.get("heatmap_done")), "pathways": bool(s.results.get("maps_done")),
            "kinetics": bool(s.results.get("race_done")),
            "export": bool(s.results.get("exported"))}.get(step, False)


def require(*steps) -> bool:
    names = dict(STEPS)
    for k in steps:
        if not done(k):
            st.info(f"Complete **{names[k]}** first.", icon="ℹ️")
            if st.button(f"Go to {names[k]}", key=f"goto_{k}_{st.session_state.step}"):
                st.session_state.step = k
                st.rerun()
            return False
    return True


# ---------------------------------------------------------------------------
# cached heavy computations
# ---------------------------------------------------------------------------
@st.cache_data(show_spinner=False, max_entries=4)
def run_reorg(payload: tuple, keep_partial: bool, prefer_json: str, min_area: float = 1000.0):
    import json
    files = [(name, io.BytesIO(b)) for name, b in payload]
    return process_files(files, keep_partial=keep_partial, prefer=json.loads(prefer_json) or None,
                         min_area=min_area)


def design_final() -> pd.DataFrame:
    return dz.finalize_design(ss().design)


def analysed_samples() -> list:
    return dz.analysed(design_final()).Sample.tolist()


def labeling():
    """(fractions long, summary long) for analysed + reference samples; recomputed when settings change."""
    s = ss()
    D = design_final()
    samples = D[D.Include.astype(bool) | (D.SampleType != TYPE_LAB)].Sample.tolist()
    c = {**dz.CorrectionSettings().to_dict(), **s.correction}
    sig = (tuple(samples), c["mode"], c["resolution"], c["method"], c["instrument"], float(c["resolving_power"]),
           float(c["ref_mz"]), s.tracer["element"], round(float(s.tracer["purity"]), 4), id(s.ds))
    if s.labeling is None or s.labeling_sig != sig:
        with st.spinner("Calculating isotopologue fractions…"):
            F, S = lb.compute_from_ct(s.ds.labeling_table, samples, s.tracer["element"], c["mode"],
                                      float(s.tracer["purity"]), c["resolution"], c["method"],
                                      float(c["resolving_power"]), float(c["ref_mz"]), c["instrument"])
        s.labeling, s.labeling_sig = (F, S), sig
    return s.labeling


def interference(S: pd.DataFrame | None = None, D: pd.DataFrame | None = None, threshold: float = 10.0) -> list:
    """Metabolites that look labeled in the unlabeled reference (likely isobaric interference)."""
    D = D if D is not None else design_final()
    if S is None:
        S = labeling()[1]
    refs = dz.reference_samples(D)
    if not refs:
        return []
    r = S[S.Sample.isin(refs)].groupby("Compound").Incorporation.max()
    return sorted(r[r > threshold].index)


def annotation():
    s = ss()
    if s.annotation is None and s.ds is not None:
        A = pw.map_metabolites(s.ds.metabolites())
        if s.ds.annotations is not None:
            O = s.ds.annotations.copy()
            for c in ("KEGG", "HMDB", "PubChemCID", "ChEBI"):
                if c in O:
                    O[c] = O[c].fillna("").astype(str)
            O = O[O.drop(columns="Compound").apply(lambda r: any(str(v).strip() for v in r), axis=1)]
            if len(O):
                A = pw.apply_overrides(A, O)
                if "BiGG" in O:
                    A = A.merge(O[["Compound", "BiGG"]], on="Compound", how="left")
        s.annotation = A
    return s.annotation


def label_of(comp: str) -> str:
    A = annotation()
    if A is None:
        return comp
    r = A.set_index("Compound")
    if comp in r.index and r.loc[comp, "Confidence"] in ("High", "Medium", "User-confirmed") and r.loc[comp, "Label"]:
        return r.loc[comp, "Label"]
    return comp


def pathway_of() -> dict:
    """First assigned pathway per metabolite (for grouping); unmapped -> 'Other / not mapped'."""
    A = annotation()
    out = {}
    for r in A.itertuples():
        p = str(r.Pathways).split(";")[0].strip() if isinstance(r.Pathways, str) and r.Pathways.strip() else ""
        out[r.Compound] = p if p and r.Confidence in ("High", "Medium", "User-confirmed") else "Other / not mapped"
    return out


def cond_means(metric: str = "incorporation") -> pd.DataFrame:
    F, S = labeling()
    D = design_final()
    X = S.merge(dz.analysed(D)[["Sample", "Condition"]], on="Sample")
    M = X.groupby(["Compound", "Condition"])[METRIC_COL[metric]].mean().unstack()
    return M.reindex(columns=dz.conditions_of(D))


def raw_matrix() -> pd.DataFrame:
    """Measured isotopologue areas (rows = 'Compound | M+i'), analysed labeled samples; nothing imputed or summed."""
    s = ss()
    return s.ds.channel_matrix(analysed_samples())


def intensity(choice: str | None = None) -> pd.DataFrame | None:
    """Channel-level intensity matrix (log2 or log2+IQR) of analysed samples."""
    s = ss()
    if s.norm is None:
        return None
    c = choice or s.intensity_choice                             # sidebar switch (default log2 + IQR)
    return s.norm["L"] if c == "log2" else s.norm["N"]


def layer(key: str):
    """
    Returns (X, info) with X = features × analysed samples and info = dict(log, unit, kind, title).
    Intensity layers are per isotopologue channel ('Compound | M+i'); fraction is per channel; incorporation,
    enrichment and mean labels are per metabolite.
    """
    s = ss()
    D = design_final()
    samples = dz.analysed(D).Sample.tolist()
    if key.startswith("int"):
        ch = {"int": s.intensity_choice, "int_log2": "log2", "int_iqr": "log2_iqr"}[key]
        X = intensity(ch)
        return X[[c for c in samples if c in X.columns]], dict(log=True, unit=INT_LABEL[ch], kind="channel",
                                                                title=INT_LABEL[ch])
    F, S = labeling()
    bad = set(interference(S, D))
    if key == "fraction":
        FM = lb.fraction_matrix(F, samples) * 100
        FM = FM[~FM.index.str.split(" \\| ").str[0].isin(bad)]
        return FM, dict(log=False, unit="fraction (%)", kind="channel", title="Isotopologue fraction (%)")
    M = lb.metric_matrix(S, {"incorporation": "incorporation", "enrichment": "enrichment",
                             "mean_labels": "mean_labels"}[key], samples)
    M = M.drop(index=[b for b in bad if b in M.index])
    return M, dict(log=False, unit=METRIC_NAME[key], kind="metabolite", title=METRIC_NAME[key])


def split_channel(idx) -> pd.DataFrame:
    s = pd.Series(list(idx), dtype=str)
    parts = s.str.rsplit(" | ", n=1, expand=True)
    if parts.shape[1] == 1:
        parts[1] = ""
    return pd.DataFrame({"Compound": parts[0], "Isotopologue": parts[1].fillna("")})


# ---------------------------------------------------------------------------
# downloads and display helpers
# ---------------------------------------------------------------------------
def to_excel(sheets: dict) -> bytes:
    buf = io.BytesIO()
    with pd.ExcelWriter(buf, engine="openpyxl") as xw:
        for name, df in sheets.items():
            if df is None or (hasattr(df, "empty") and df.empty):
                continue
            sheet = "".join("-" if ch in '[]:*?/\\' else ch for ch in str(name))[:31] or "Sheet"
            base, k = sheet, 2
            while sheet in xw.sheets:
                sheet = f"{base[:28]}_{k}"; k += 1
            df.to_excel(xw, sheet_name=sheet, index=not isinstance(df.index, pd.RangeIndex))
    return buf.getvalue()


XLSX_MIME = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
MAX_PIXELS = 60e6          # raster export budget (e.g. 10 × 16.7 in at 600 DPI)
DPI_OPTS = [300, 600, 1200]


def _box():
    import inspect
    try:
        if "border" in inspect.signature(st.container).parameters:
            return st.container(border=True)
    except Exception:
        pass
    return st.container()


def _safe(stem: str) -> str:
    return "".join(ch if (ch.isalnum() or ch in "-_.") else "_" for ch in str(stem))[:80].strip("_") or "table"


def download_panel(stem: str, key: str, make_fig=None, default_size=(8.0, 6.0), tables: dict | None = None,
                   workbook: dict | None = None, caption: str = "", add_to_report: bool = True,
                   title: str = "⬇ Downloads"):
    """
    One consistent download box used under every figure / table of the app.
      Figure row : format · resolution (DPI, default 600) · width · height (inches) → Render → Download
      Table rows : one line per table with CSV and Excel buttons (direct download)
      Workbook   : optional multi-sheet Excel (all tables of the step in one file)
    """
    stem = _safe(stem)
    with _box():
        st.markdown(f"**{title}**")
        if make_fig is not None:
            c = st.columns([1.1, 1.1, 0.9, 0.9, 1.3, 1.9])
            fmt = c[0].selectbox("Figure format", EXPORT_FORMATS, key=f"{key}_fmt",
                                 help="PNG / TIFF / JPEG = high-resolution raster; PDF / SVG = vector "
                                      "(infinitely sharp, editable text).")
            dpi = c[1].selectbox("Resolution (DPI)", DPI_OPTS, index=1, key=f"{key}_dpi",
                                 disabled=fmt in ("PDF", "SVG"))
            # a step may preset the size through Session State (e.g. the heatmap size controls); then no default
            vw = {} if f"{key}_w" in st.session_state else dict(value=float(min(max(default_size[0], 2.0), 60.0)))
            vh = {} if f"{key}_h" in st.session_state else dict(value=float(min(max(default_size[1], 2.0), 200.0)))
            w = c[2].number_input("Width (in)", 2.0, 60.0, step=0.5, key=f"{key}_w", **vw)
            h = c[3].number_input("Height (in)", 2.0, 200.0, step=0.5, key=f"{key}_h", **vh)
            c[4].markdown("<div style='height:1.75rem'></div>", unsafe_allow_html=True)
            if c[4].button("Render figure", key=f"{key}_prep", type="primary"):
                with st.spinner("Rendering high-resolution figure…"):
                    import matplotlib.pyplot as plt
                    try:
                        fig = make_fig(w, h)
                    except Exception as e:
                        st.error(f"Could not draw the figure: {e}")
                        fig = None
                    if fig is not None:
                        wi, hi = fig.get_size_inches()
                        # very tall / wide figures: keep raster exports within a safe pixel budget
                        # (memory on small servers); vector formats are unaffected
                        eff = dpi
                        if fmt not in ("PDF", "SVG") and wi * hi * dpi ** 2 > MAX_PIXELS:
                            eff = max(72, int((MAX_PIXELS / (wi * hi)) ** 0.5))
                        data = fig_bytes(fig, fmt, eff if fmt not in ("PDF", "SVG") else 300)
                        pv = min(150, max(40, int((8e6 / (wi * hi)) ** 0.5)))
                        preview = fig_bytes(fig, "PNG", pv)
                        if add_to_report:
                            rp = min(200, max(40, int((12e6 / (wi * hi)) ** 0.5)))
                            ss().report_figs = [x for x in ss().report_figs if x[0] != (caption or stem)] + \
                                               [(caption or stem, fig_bytes(fig, "PNG", rp))]
                        plt.close(fig)
                        st.session_state[f"{key}_data"] = (data, fmt, preview, eff, wi, hi)
                        if eff < dpi:
                            st.session_state[f"{key}_note"] = (f"Figure is {wi:.1f} × {hi:.1f} in: resolution lowered to "
                                                               f"{eff} DPI to stay within {MAX_PIXELS / 1e6:.0f} "
                                                               "million pixels. Use PDF or SVG for full vector "
                                                               "quality, or a smaller height.")
                        else:
                            st.session_state.pop(f"{key}_note", None)
            if st.session_state.get(f"{key}_data"):
                data, f, preview, d_, wi, hi = st.session_state[f"{key}_data"]
                ext = f.lower().replace("jpeg", "jpg")
                c[5].markdown("<div style='height:1.75rem'></div>", unsafe_allow_html=True)
                c[5].download_button(f"⬇ {stem}.{ext}", data, f"{stem}.{ext}", MIME[f], key=f"{key}_dl")
                size = ("vector" if f in ("PDF", "SVG") else f"{int(wi * d_):,} × {int(hi * d_):,} px at {d_} DPI")
                st.caption(f"Ready: {f} · {wi:.1f} × {hi:.1f} in · {size} · {len(data) / 1e6:.1f} MB"
                           + (" · added to the report (Step 10)" if add_to_report else ""))
                if st.session_state.get(f"{key}_note"):
                    st.caption("ℹ " + st.session_state[f"{key}_note"])
                with st.expander("Preview", expanded=False):
                    show_image(preview)
        for name, df in (tables or {}).items():
            if df is None or (hasattr(df, "empty") and df.empty):
                continue
            k = _safe(name)
            c = st.columns([3.2, 1, 1, 1.8])
            c[0].markdown(f"<div style='padding-top:0.45rem'>Table · <b>{name}</b> "
                          f"<span style='color:#000'>({len(df):,} rows)</span></div>", unsafe_allow_html=True)
            c[1].download_button("⬇ CSV", df.to_csv(index=not isinstance(df.index, pd.RangeIndex)).encode(),
                                 f"{stem}_{k}.csv" if len(tables) > 1 else f"{stem}.csv", "text/csv",
                                 key=f"{key}_{k}_csv")
            c[2].download_button("⬇ Excel", to_excel({name: df}),
                                 f"{stem}_{k}.xlsx" if len(tables) > 1 else f"{stem}.xlsx", XLSX_MIME,
                                 key=f"{key}_{k}_xlsx")
        if workbook:
            c = st.columns([3.2, 2, 1.8])
            c[0].markdown("<div style='padding-top:0.45rem'>Workbook · <b>all tables of this view</b> "
                          f"<span style='color:#000'>({len(workbook)} sheets)</span></div>", unsafe_allow_html=True)
            c[1].download_button("⬇ Excel workbook", to_excel(workbook), f"{stem}_all.xlsx", XLSX_MIME,
                                 key=f"{key}_wb")


def table_downloads(df: pd.DataFrame, stem: str, key: str, name: str | None = None):
    download_panel(stem, key, tables={name or "Values": df})


def export_widget(make_fig, stem: str, key: str, default_size=(8.0, 6.0), add_to_report=True, caption="",
                  tables: dict | None = None, workbook: dict | None = None, title: str = "⬇ Downloads"):
    download_panel(stem, key, make_fig=make_fig, default_size=default_size, tables=tables, workbook=workbook,
                   caption=caption, add_to_report=add_to_report, title=title)


def show_image(data: bytes):
    for kw in (dict(width="stretch"), dict(use_container_width=True), dict(use_column_width=True)):
        try:
            st.image(data, **kw)
            return
        except Exception:
            continue
    st.image(data)


def show_df(df, height=None, **kw):
    try:
        st.dataframe(df, width="stretch", height=height or "auto", **kw)
    except Exception:
        st.dataframe(df, use_container_width=True, height=height, **kw)


def show_plotly(fig, key=None, select=False):
    cfg = {"displaylogo": False, "toImageButtonOptions": {"format": "png", "scale": 5}, "scrollZoom": True}
    sel = dict(on_select="rerun", selection_mode=("points",)) if select else {}
    # white plot area and solid black text everywhere (the app's own style, not the Streamlit chart theme)
    fig.update_layout(template="plotly_white", plot_bgcolor=fig.layout.plot_bgcolor or "white",
                      paper_bgcolor=fig.layout.paper_bgcolor or "white", font_color="black")
    h = getattr(fig.layout, "height", None)
    sizes = ([dict(width="stretch", height=int(h))] if h else []) + [dict(width="stretch"), dict(use_container_width=True)]
    for size in sizes:
        try:
            return st.plotly_chart(fig, key=key, config=cfg, theme=None, **size, **sel)
        except Exception:
            continue
    return st.plotly_chart(fig, key=key, config=cfg, theme=None, use_container_width=True)


def show_mpl(fig):
    import matplotlib.pyplot as plt
    show_image(fig_bytes(fig, "PNG", 200))
    plt.close(fig)


def direction_banner(group_a, group_b, log_scale=True):
    st.markdown(
        f"<div class='direction'>Comparison: <b>{group_b}</b> vs <b>{group_a}</b> &nbsp;·&nbsp; "
        f"{'log2FC = mean(' + group_b + ') − mean(' + group_a + ')' if log_scale else 'difference = ' + group_b + ' − ' + group_a}"
        f" &nbsp;·&nbsp; <span class='pos'>positive = higher in {group_b}</span></div>", unsafe_allow_html=True)


def sig_stars(q):
    if q is None or not np.isfinite(q):
        return ""
    return "***" if q < 0.001 else "**" if q < 0.01 else "*" if q < 0.05 else ""


def zip_bytes(files: dict) -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for name, data in files.items():
            z.writestr(name, data)
    return buf.getvalue()


def sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def header(n: int, title: str, desc: str):
    st.markdown(f"<div class='step-head'><div class='step-kicker'>Step {n} of {N_STEPS}</div>"
                f"<div class='step-title'>{title}</div><p class='step-desc'>{desc}</p></div>", unsafe_allow_html=True)


def note(text: str, ok: bool = False):
    st.markdown(f"<div class='{'okbox' if ok else 'note'}'>{text}</div>", unsafe_allow_html=True)


def next_button(to: str, label: str):
    st.write("")
    if st.button(f"Continue to {label} →", type="primary", key=f"next_{to}"):
        st.session_state.step = to
        st.rerun()


def arms_times():
    D = design_final()
    return dz.groups_of(D), dz.times_of(D)
