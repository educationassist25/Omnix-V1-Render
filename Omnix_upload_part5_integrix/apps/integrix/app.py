"""
Omnix Integrix - tri-omics integration (transcriptomics, proteomics, metabolomics) (standalone Streamlit app, v1).

Run:  streamlit run app.py
"""
from __future__ import annotations

import io
import json
import os
import sys
import time

import numpy as np
import pandas as pd
import streamlit as st

# make the app's own folder importable however it is launched (streamlit run from another folder, an Omnix
# multipage portal, st.navigation, Render/Cloud working directories)
_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from trace_core.compat import patch_streamlit  # noqa: E402

ST_HAS = patch_streamlit(st)      # older Streamlit: unsupported options are dropped instead of crashing
from trace_core import __version__
from trace_core import figures as F
from trace_core import report as R
from trace_core.cascade import PATTERN_HELP
from trace_core.data import LAYER_NAMES, build_design, cramers_v, make_layer, read_table
from trace_core.engine import run_trace
from trace_core.kb import load_kb
from trace_core import ui as UI

HERE = os.path.dirname(os.path.abspath(__file__))
DEMO = os.path.join(HERE, "demo_data")
LOGO = os.path.join(HERE, "assets", "omnix_logo_web.png")

# Inside the Omnix portal (omnix_portal.portal.run_app sets OMNIX = True) the portal owns the page title, icon and
# header, decides who may analyse their own data, and the modules can hand their normalized data to Integrix.
OMNIX = globals().get("OMNIX", False)
if not OMNIX:
    st.set_page_config(page_title="Omnix™ Integrix", page_icon=LOGO if os.path.exists(LOGO) else None,
                       layout="wide")
INBOX_KEY = "omnix_integrix_inbox"     # portal-level session key (survives switching between Omnix apps)
DEMO_ONLY = bool(OMNIX and st.session_state.get("omnix_demo_only", False))
_PORTAL_STOP = globals().get("OMNIX_STOP")   # inside the portal: end this page without stopping the portal's side panel


def _stop():
    if _PORTAL_STOP is not None:
        raise _PORTAL_STOP()
    st.stop()
try:
    _theme = st.context.theme.type
except Exception:  # older Streamlit
    _theme = None
# the viewer's light / dark choice (menu > Settings, or the system setting); remembered in case a rerun
# reports no theme
if _theme:
    st.session_state["_mode"] = _theme
MODE = "dark" if st.session_state.get("_mode", "light") == "dark" else "light"
st.markdown(UI.css(MODE), unsafe_allow_html=True)
if OMNIX:      # the same white page as the other Omnix platforms
    st.markdown("<style>.stApp{background:#FFFFFF !important}</style>", unsafe_allow_html=True)


@st.cache_data(show_spinner=False)
def _demo_csv(name: str) -> pd.DataFrame:
    return pd.read_csv(os.path.join(DEMO, name))


@st.cache_data(show_spinner=False)
def _demo_bytes(name: str) -> bytes:
    return open(os.path.join(DEMO, name), "rb").read()


@st.cache_resource(show_spinner="Loading knowledge base")
def get_kb():
    return load_kb()


_SHOWN: dict = {}       # charts drawn in this run, per name: every chart gets its own key


def show(fig, name="trace_figure", notitle=False, width_px: int = 1100, pub=None):
    """Draw a Plotly chart and, under it, the download panel for the same figure.

    width_px: about how wide the chart is on screen (full width ~1100 px, half width ~620 px); the download uses it
    as its default width, so the file has the screen's proportions. pub: (builder, w, h, fs) of the
    publication-style figure, used only when the screen figure cannot be exported on this server."""
    cfg = F.config(st.session_state.get("fig_fmt", "png"), st.session_state.get("fig_scale", 6), name)
    # An explicit key: two charts with the same figure on one page (e.g. a pathway chart shown in two sections)
    # would otherwise get the same automatic ID, which Streamlit 1.52 rejects (StreamlitDuplicateElementId).
    _SHOWN[name] = _SHOWN.get(name, 0) + 1
    key = f"plt_{name}_{_SHOWN[name]}"
    base = fig
    if notitle:   # the surrounding card already carries the title
        base = go_figure(fig)
        base.update_layout(title_text="", margin_t=min(40, fig.layout.margin.t or 40))
    fig = UI.theme_figure(base, MODE)
    try:
        st.plotly_chart(fig, width="stretch", config=cfg, theme=None, key=key)
    except TypeError:  # older Streamlit
        st.plotly_chart(fig, use_container_width=True, config=cfg, theme=None, key=key)
    screen_download_panel(f"dl_{name}_{_SHOWN[name]}", base, name, width_px, pub)


def go_figure(fig):
    import plotly.graph_objects as go
    return go.Figure(fig)


def _blacken():
    """The portal's black-text rule (omnix_portal.figure_style), so files match the screen inside Omnix."""
    try:
        from omnix_portal.figure_style import blacken_plotly
        return blacken_plotly
    except Exception:
        return None


def screen_download_panel(key: str, fig, name: str, width_px: int = 1100, pub=None):
    """Publication download of the chart on screen: the same chart (data, marks, colours) typeset for a journal
    (trace_core.printfig): column width, real point sizes, height from the content, tight margins, no title.
    PNG · TIFF · JPEG · SVG · PDF · EPS or all six, 300-1200 DPI. Built only when asked for."""
    from trace_core import plotly_export as PX
    from trace_core import printfig as PRT
    with st.expander("⬇ Download publication figure (PNG · TIFF · JPEG · SVG · PDF · EPS)"):
        if not PX.available():
            if pub is not None:
                st.caption("The typeset version of this chart cannot be drawn on this server right now; this is "
                           "the matplotlib version of the same chart.")
                download_panel(key + "_pub", *pub[:4], name=name, expander=False)
            else:
                st.caption("Use the camera icon at the top right of the chart to save it as PNG or SVG.")
            return
        widths = PRT.widths_for(name)
        sizes = list(widths)
        d_w = PRT.DEFAULT_WIDTH.get(name, 5.5)
        d_i = next((i for i, k in enumerate(sizes) if widths[k] == d_w), len(sizes) - 1)
        c = st.columns(3) + st.columns(3)
        size = c[0].selectbox("Width", sizes, index=d_i, key=f"{key}_size",
                              help="Journal column widths. The chart is laid out again at this width; text keeps "
                                   "its point size, so it prints at exactly that size.")
        font = c[1].selectbox("Text size (pt)", [7, 8, 9, 10, 12], index=1, key=f"{key}_pt",
                              help="Size of labels and axis titles in the printed figure; tick labels and legend "
                                   "are 1 pt smaller (never below 6 pt).")
        hh = c[2].number_input("Height (in, 0 = auto)", 0.0, 30.0, 0.0, 0.25, key=f"{key}_h")
        fmt = c[3].selectbox("Format", list(PX.FORMATS) + ["All formats (ZIP)"], key=f"{key}_fmt")
        dpi = c[4].selectbox("Resolution (DPI)", [300, 600, 1200], index=1, key=f"{key}_dpi")
        transp = c[5].checkbox("Transparent background", key=f"{key}_tr",
                               help="PNG, SVG and PDF only; TIFF, JPEG and EPS are always on white.")
        if st.button("Prepare file", key=f"{key}_go", type="secondary"):
            with st.spinner("Preparing the file"):
                sf = PX.screen_figure(fig, UI.theme_figure, _blacken(), transp)
                out, ww, hh_ = PRT.typeset(sf, name, widths[size], font, hh or None)
                if transp:
                    out.update_layout(paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)")
                if fmt.startswith("All"):
                    data, fn = PX.all_formats_zip(out, name, dpi, ww, hh_, transp), f"{name}_all_formats.zip"
                else:
                    ext = PX.FORMATS[fmt]
                    data, fn = PX.figure_bytes(out, ext, dpi, ww, hh_, transp), f"{name}.{ext}"
                eff = PX.effective_dpi(ww, hh_, dpi)
                note = f"{ww:g} × {hh_:.2f} in · text {font} pt"
                if fmt not in ("SVG", "PDF", "EPS"):
                    note += f" · {int(ww * eff):,} × {int(hh_ * eff):,} px at {eff} DPI"
                    if eff < dpi:
                        note += " (lowered for this large figure; SVG or PDF keep full quality)"
                st.session_state[f"_scr_{key}"] = (fn, data, note)
        if f"_scr_{key}" in st.session_state:
            fn, data, note = st.session_state[f"_scr_{key}"]
            st.download_button(f"Download {fn} ({len(data) / 1e6:.1f} MB)", data, fn,
                               mime=PX.MIME[fn.rsplit(".", 1)[1]], key=f"{key}_btn", type="primary",
                               icon=":material/download:")
            st.caption(note)


def show_mpl(fig, max_scale: float = 1.15, cache_key=None):
    """Show a matplotlib figure as a crisp vector image (SVG, text converted to outlines so it renders exactly as
    laid out, at any zoom or screen density). It scales down with the window and never beyond 115 % of its
    natural size. `fig` may be a figure or a zero-argument builder; with `cache_key` the rendered SVG is kept for
    this session, so a rerun (any click) does not redraw it. Returns the natural size in inches."""
    import base64
    import matplotlib.pyplot as plt
    cache = st.session_state.setdefault("_svg_cache", {})
    key = (cache_key, MODE) if cache_key is not None else None
    if key is not None and key in cache:
        svg, w_in, h_in = cache[key]
    else:
        fig = fig() if callable(fig) else fig
        w_in, h_in = fig.get_size_inches()
        buf = io.BytesIO()
        with plt.rc_context({"svg.fonttype": "path"}):
            fig.savefig(buf, format="svg", bbox_inches="tight", pad_inches=0.05)
        plt.close(fig)
        svg = UI.fluid_svg(buf.getvalue())
        if MODE == "dark":
            svg = UI.darken_svg(svg)
        if key is not None:
            if len(cache) > 24:
                cache.pop(next(iter(cache)))
            cache[key] = (svg, float(w_in), float(h_in))
    b64 = base64.b64encode(svg).decode()
    st.markdown(f'<div style="width:100%;text-align:center"><img src="data:image/svg+xml;base64,{b64}" '
                f'style="width:100%;max-width:{int(w_in * 96 * max_scale)}px;height:auto;display:inline-block"/></div>',
                unsafe_allow_html=True)
    return float(w_in), float(h_in)


def _is_p(c) -> bool:
    c = str(c)                    # case matters: "P" is the product balance, "p" a p-value
    return c in ("p", "q", "fdr", "FDR", "padj") or c.startswith(("p_", "q_")) or c.endswith(("_p", "_q"))


def _fmt_p(v):
    if pd.isna(v):
        return "—"
    v = float(v)
    return "0" if v == 0 else (f"{v:.1e}" if v < 1e-3 else f"{v:.4f}")


def _fmt_num(v):
    if pd.isna(v):
        return "—"
    v = float(v)
    if v == int(v) and abs(v) < 1e15:
        return f"{int(v)}"
    return f"{v:.4f}".rstrip("0").rstrip(".") if abs(v) >= 1e-4 else f"{v:.2e}"


def table(df, height=380, key=None):
    """Tables: p- and q-values in scientific notation when small (never rounded to 0), other numbers to 4
    decimals, missing values as a dash (not 'None'), no row-number column, and no empty rows below short tables."""
    df = df.reset_index(drop=True)
    fmt = {c: (_fmt_p if _is_p(c) else _fmt_num) for c in df.columns if pd.api.types.is_numeric_dtype(df[c])
           and not pd.api.types.is_bool_dtype(df[c])}
    # Streamlit shows any empty cell as "None" whatever the formatting: columns with gaps are shown as text
    # (a dash for "not measured"); columns without gaps stay numeric and sort as numbers
    df = df.copy()
    for c in list(fmt):
        if df[c].isna().any():
            df[c] = df[c].map(fmt.pop(c))
    for c in df.columns:
        if df[c].dtype == object and df[c].isna().any():
            df[c] = df[c].where(df[c].notna(), "—")
    h = int(min(height, 38 + 35 * max(len(df), 1) + 2))
    if key is None:            # same reason as show(): identical tables on one page need their own IDs
        _SHOWN["__table"] = _SHOWN.get("__table", 0) + 1
        key = f"tbl_auto_{_SHOWN['__table']}"
    try:
        st.dataframe(df.style.format(fmt, na_rep="—"), height=h, hide_index=True, key=key)
    except Exception:          # Styler not available for this table: plain display
        st.dataframe(df, height=h, hide_index=True, key=key)


def dl_card(title, file_stem, desc="", key="", buttons=()):
    """One download card, the same style in every Omnix platform: title, a short description, the file name and
    the buttons (label, data, file name, mime, key). The files are unchanged."""
    import html as _h
    with st.container(key=f"omxdl_{key}", border=True):
        st.markdown(f"<div class='omx-dl-t'>{_h.escape(title)}</div>"
                    + (f"<div class='omx-dl-d'>{_h.escape(desc)}</div>" if desc else "")
                    + f"<div class='omx-dl-f'>{_h.escape(file_stem)}</div>", unsafe_allow_html=True)
        with st.container(horizontal=True, gap="small"):
            for lab, data, fn, mime, bk in buttons:
                st.download_button(lab, data, fn, mime=mime, key=bk, icon=":material/download:", help=fn)


def _xlsx_of(df, sheet):
    buf = io.BytesIO()
    with pd.ExcelWriter(buf, engine="openpyxl") as w:
        df.to_excel(w, sheet_name=str(sheet)[:31] or "Table", index=False)
    return buf.getvalue()


def csv_button(df, name, label="Download table (CSV)", key=None):
    """A table as a download card (CSV and Excel). The CSV is the same file as before."""
    import functools
    import re
    stem = name.rsplit(".", 1)[0]
    title = label if label != "Download table (CSV)" else stem.replace("_", " ").capitalize()
    title = re.sub(r"^Download\s+", "", title)
    bk = key or "dl_" + re.sub(r"[^A-Za-z0-9]+", "_", stem)
    dl_card(title[:1].upper() + title[1:], stem, f"{len(df):,} rows", key=bk, buttons=[
        ("CSV", df.to_csv(index=False).encode(), name, "text/csv", bk),
        ("Excel", functools.partial(_xlsx_of, df, stem), f"{stem}.xlsx", "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", f"{bk}_xlsx")])


def _fonts():
    import matplotlib.font_manager as fm
    have = {f.name for f in fm.fontManager.ttflist}
    pref = ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans", "FreeSans", "Inter", "Carlito",
            "Times New Roman", "Liberation Serif", "DejaVu Serif", "FreeSerif"]
    return [f for f in pref if f in have] or ["DejaVu Sans"]


def download_panel(key: str, builder, w: float, h: float | None, fs: float = 7.0, name: str = "figure",
                   expander: bool = True):
    """Publication export for one figure: any format, 300-1200 DPI, custom size / font, or all formats
    at once. `builder(w, h, fs)` returns a matplotlib figure (rendered only when requested)."""
    import contextlib
    from trace_core import pubfig as PF
    box = st.expander("⬇ Download publication-quality figure (PNG · TIFF · JPEG · SVG · PDF · EPS)") if expander \
        else contextlib.nullcontext()
    with box:
        c = st.columns(3) + st.columns(3)          # two rows of three: readable at half-screen width
        fmt = c[0].selectbox("Format", list(PF.FORMATS) + ["All formats (ZIP)"], key=f"{key}_fmt")
        dpi = c[1].selectbox("Resolution (DPI)", [300, 600, 1200], index=1, key=f"{key}_dpi")
        font = c[2].selectbox("Font family", _fonts(), key=f"{key}_font")
        ww = c[3].number_input("Width (in)", 2.0, 30.0, float(w), 0.5, key=f"{key}_w")
        hh = c[4].number_input("Height (in, 0 = auto)", 0.0, 30.0, float(h or 0), 0.5, key=f"{key}_h")
        ff = c[5].number_input("Font (pt)", 4.0, 20.0, float(fs), 0.5, key=f"{key}_fs")
        t1, t2 = st.columns([1, 1])
        transp = t1.checkbox("Transparent background", key=f"{key}_tr")
        if t2.button("Prepare file", key=f"{key}_go"):
            with st.spinner("Rendering"):
                PF.style(font, ff)
                fmts = list(PF.FORMATS.values()) if fmt.startswith("All") else [PF.FORMATS[fmt]]
                out, low = {}, None
                for f in fmts:
                    fig = builder(ww, hh or None, ff)
                    wi, hi = fig.get_size_inches()
                    eff = dpi
                    if wi * hi * dpi ** 2 > 60e6:          # raster budget (memory), as in every Omnix platform
                        eff = max(72, int((60e6 / (wi * hi)) ** 0.5))
                        low = eff
                    out[f"{name}.{f}"] = PF.export(fig, f, eff, transp)
                st.session_state[f"_pubnote_{key}"] = (
                    f"Raster files lowered to {low} DPI for this large figure (60 million pixels at most); "
                    "PDF, SVG and EPS are vector files." if low else None)
                if len(out) == 1:
                    (fn, data), = out.items()
                    st.session_state[f"_pub_{key}"] = (fn, data, PF.MIME[fn.rsplit('.', 1)[1]])
                else:
                    st.session_state[f"_pub_{key}"] = (f"{name}_all_formats.zip", R.bundle(out), "application/zip")
        if f"_pub_{key}" in st.session_state:
            fn, data, mime = st.session_state[f"_pub_{key}"]
            st.download_button(f"Download {fn} ({len(data) / 1e6:.1f} MB)", data, fn, mime=mime, key=f"{key}_btn",
                               type="primary")
            if st.session_state.get(f"_pubnote_{key}"):
                st.caption(st.session_state[f"_pubnote_{key}"])


kb = get_kb()

# ===================================================================================== setup (sidebar)
NORMDEMO = "Demo: Omnix normalized"
FACTORIAL = "Demo: Raw intensities"
OMNIX_UP = "Upload Omnix exports"          # the Omnix modules' "Download for Integrix" zips (Omnix platform build)
SESSION = "From this Omnix session"        # data sent with "Send to Integrix" in the Omnix modules
UPLOAD = "Upload my data"
SRC_OPTS = [NORMDEMO, FACTORIAL, UPLOAD]
SRC_CAPS = ["Condition × Treatment × Sex · 40 samples, log2", "the same 40 samples before normalization",
            "your own tables: raw or normalized (log2)"]
if OMNIX:
    _inbox_n = len(st.session_state.get(INBOX_KEY) or {})
    SRC_OPTS = [SESSION, NORMDEMO, FACTORIAL, OMNIX_UP, UPLOAD]
    SRC_CAPS = [(f"{_inbox_n} dataset{'s' if _inbox_n != 1 else ''} sent from the Omnix modules" if _inbox_n else
                 "nothing sent yet: use Send to Integrix in a module") + (" · subscription" if DEMO_ONLY else ""),
                *SRC_CAPS[:2], "the zips from Download for Integrix in each module"
                + (" · subscription" if DEMO_ONLY else ""), SRC_CAPS[2] + (" · subscription" if DEMO_ONLY else "")]
UP_KINDS = ["Raw intensities", "Normalized (log2)"]
CASE_PREF = ["Tumor", "Cancer", "Treatment", "Female", "Cancer Treatment"]
REF_PREF = ["Control", "Normal", "Male", "Cancer Control", "Normal Control"]


def _pick_case(groups: list) -> str:
    return next((g for g in CASE_PREF if g in groups), groups[-1] if groups else "")


def _pick_ref(groups: list, case: str) -> str:
    return next((g for g in REF_PREF if g in groups and g != case), next((g for g in groups if g != case), ""))


def _reset_groups():
    """The comparison column changed: choose Case / Reference again from its own levels."""
    for k_ in ("des_case", "des_ref"):
        st.session_state.pop(k_, None)
    _keep_open("step_design")


def _keep_open(step: str):
    """A setting changed the step's title (its summary): keep that step open on older Streamlit, where a new
    title would otherwise close it."""
    st.session_state["_open_step"] = step


README_DEMO = {
    NORMDEMO: """Omnix Integrix - example data: Omnix normalized, Condition x Treatment x Sex (SYNTHETIC)
======================================================================================
Synthetic data made for learning and testing the app. Not patient data.
This is the format the Omnix modules export for Integrix ("Omnix export contract" 1.0).

Files
  transcriptomics.csv + transcriptomics_omnix.json   genes x 40 samples, log2, from BulkRNAAI Pro
  proteomics.csv      + proteomics_omnix.json        proteins x 40 samples, log2, from ProteoAI Pro
  metabolomics.csv    + metabolomics_omnix.json      metabolites x 40 samples, PQN + log2, from MetaboAI Pro
  metadata.csv        Sample, Group, Condition, Treatment, Sex, Batch, Age (the same samples in every module)
  ground_truth.json   the biology planted in the data (the answers)

Matrices: one row per feature (identifier columns first), one column per sample, values already normalized on
log2 scale. Missing values are left empty (no imputation). No batch correction was applied (metabolomics batch C
has a shift: try Study design > Batch-correct Metabolomics).

Each *_omnix.json says how its matrix was made: module, scale, normalization, imputation, batch correction.
Integrix reads it and skips those steps; it starts at the per-layer statistics and the integration.
To use your own normalized tables: Data > Upload my data > Normalized (log2).

Comparisons to try (Study design > Group column, then Case over Reference)
  Condition  Cancer over Normal     glycolysis up, ASS1 bottleneck (citrulline up), tryptophan -> kynurenine,
                                    fatty-acid-oxidation proteins up without mRNA, cholesterol mRNA up only
  Treatment  Treatment over Control  glutaminase blocked: glutamine up, glutamate down, enzyme unchanged
  Group      Cancer Treatment over Cancer Control  the drug effect plus glycolysis partly reversed
  Sex        Female over Male       aromatase (CYP19A1): testosterone down, estrogens up; otherwise quiet
""",
    FACTORIAL: """Omnix Integrix - example data: Condition x Treatment x Sex (SYNTHETIC)
=====================================================================
Synthetic data made for learning and testing the app. Not patient data.

Files
  transcriptomics.csv  genes x 40 samples (column 'Gene' = gene symbol)
  proteomics.csv       proteins x 40 samples ('Protein' = UniProt, 'Gene' = symbol)
  metabolomics.csv     metabolites x 40 samples ('Metabolite' = name, 'HMDB' = ID where known)
  metadata.csv         one row per sample: Sample, Group, Condition, Treatment, Sex, Batch, Age
  ground_truth.json    the biology planted in the data (the answers)

Design: 2 x 2 x 2, five samples per combination. Values are raw intensities (the app log2-transforms them).
About 8% of proteomics values are missing at random; low values are censored below a detection limit in every
layer. Metabolomics batch C has a shift that is not tied to any group (try batch correction).

Comparisons to try (Study design > Group column, then Case over Reference)
  Condition  Cancer over Normal     glycolysis up, ASS1 bottleneck (citrulline up), tryptophan -> kynurenine,
                                    fatty-acid-oxidation proteins up without mRNA, cholesterol mRNA up only
  Treatment  Treatment over Control  glutaminase blocked: glutamine up, glutamate down, enzyme unchanged
  Group      Cancer Treatment over Cancer Control  the drug effect plus glycolysis partly reversed
  Sex        Female over Male       aromatase (CYP19A1): testosterone down, estrogens up; otherwise quiet

To analyze your own data, upload tables in the same layout (Data > Upload my data > Raw intensities).
""",
}


@st.cache_data(show_spinner=False)
def _example_zip(source: str) -> bytes:
    """The selected demo as a ZIP: the three omics tables, the sample metadata, a README and the ground truth."""
    import zipfile
    raw_, meta_ = _demo_tables(source)
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for k_, nm_ in (("rna", "transcriptomics.csv"), ("prot", "proteomics.csv"), ("met", "metabolomics.csv")):
            z.writestr(nm_, raw_[k_].to_csv(index=False))
        z.writestr("metadata.csv", meta_.to_csv(index=False))
        z.writestr("README.txt", README_DEMO.get(source, ""))
        for k_, info_ in _demo_infos(source).items():
            z.writestr(f"{info_['matrix'].rsplit('.', 1)[0]}_omnix.json", json.dumps(info_, indent=2))
        _gtd = {NORMDEMO: "normalized", FACTORIAL: "factorial"}.get(source, "")
        gt = os.path.join(DEMO, _gtd, "ground_truth.json" if _gtd else "demo_ground_truth.json")
        if os.path.exists(gt):
            z.write(gt, "ground_truth.json")
    return buf.getvalue()


def _demo_sub(source: str) -> str:
    return {NORMDEMO: "normalized/", FACTORIAL: "factorial/"}.get(source, "demo_")


@st.cache_data(show_spinner=False)
def _demo_infos(source: str) -> dict:
    """The Omnix sidecars of a demo (only the normalized demo has them)."""
    out = {}
    if source == NORMDEMO:
        for k_, nm_ in (("rna", "transcriptomics"), ("prot", "proteomics"), ("met", "metabolomics")):
            pth = os.path.join(DEMO, "normalized", f"{nm_}_omnix.json")
            if os.path.exists(pth):
                out[k_] = json.load(open(pth))
    return out


def _demo_tables(source: str):
    sub = _demo_sub(source)
    names = (("rna", "transcriptomics.csv"), ("prot", "proteomics.csv"), ("met", "metabolomics.csv"))
    raw_ = {k_: _demo_csv(sub + f_).copy() for k_, f_ in names}
    meta_ = _demo_csv(sub + "metadata.csv").copy()
    return raw_, meta_
SECTION_NAMES = ["Results", "Networks & targets", "Advanced", "Data & export"]
SECTION_ICONS = [":material/insights:", ":material/hub:", ":material/science:", ":material/folder_open:"]
STEP_KEYS = ["step_data", "step_design", "step_stats", "step_engines"]


def _accordion(opened: str):
    """Only one setup step open at a time."""
    if st.session_state.get(opened):
        for k in STEP_KEYS:
            if k != opened:
                st.session_state[k] = False


if OMNIX:
    _ib = st.session_state.get(INBOX_KEY) or {}
    _sig = tuple(sorted((k_, v_.get("sent", "")) for k_, v_ in _ib.items()))
    if "src" not in st.session_state:            # first visit: the data the modules sent, if any, else the demo
        st.session_state["src"] = SESSION if (len(_ib) >= 2 and not DEMO_ONLY) else NORMDEMO
        if not DEMO_ONLY:                         # subscribers see the data choices (own data) straight away
            st.session_state["_open_step"] = "step_data"
    elif _sig != st.session_state.get("_ib_seen") and _ib and not DEMO_ONLY:
        # new data sent from a module since the last visit: show it, instead of the demo chosen earlier
        if len(_ib) >= 2:
            st.session_state["src"] = SESSION
        st.session_state["_open_step"] = "step_data"
        st.session_state.pop("step_data", None)   # the Data step is drawn again, open
    st.session_state["_ib_seen"] = _sig
source = st.session_state.get("src", SRC_OPTS[0])
if source not in SRC_OPTS:
    source = SRC_OPTS[0]
if st.session_state.get("_src_prev") != source:      # another dataset: its own columns and groups
    for _k in ("des_group", "des_case", "des_ref", "des_sample", "des_batch"):
        st.session_state.pop(_k, None)
    st.session_state["_src_prev"] = source
up = {k: st.session_state.get(f"up_{k}") for k in ("rna", "prot", "met", "meta", "map")} \
    if source in ("Upload my data", OMNIX_UP) else {}

# ---------------------------------------------------------------- load raw tables
raw, meta, inputs = {}, None, {}
norm, omnix_msgs = {}, []          # norm: Omnix sidecar (or a user declaration) per layer that is already normalized
if source.startswith("Demo"):
    raw, meta = _demo_tables(source)
    _sub = _demo_sub(source)
    for f in ("transcriptomics.csv", "proteomics.csv", "metabolomics.csv"):
        inputs[_sub + f] = _demo_bytes(_sub + f)
    norm = _demo_infos(source)
elif source == SESSION:
    from trace_core.omnix_io import read_exports
    _ib = st.session_state.get(INBOX_KEY) or {}
    _use = [k_ for k_ in _ib if st.session_state.get(f"ib_use_{k_}", _ib[k_].get("default_use", True))]
    if _use and not DEMO_ONLY:
        _pairs = [(_ib[k_]["zip_name"], _ib[k_]["zip"]) for k_ in _use]
        raw, meta, norm, omnix_msgs = read_exports(_pairs)
        inputs.update({n_: b_ for n_, b_ in _pairs})
elif source == OMNIX_UP:
    from trace_core.omnix_io import read_exports
    _files = st.session_state.get("up_omnix") or []
    if _files:
        _pairs = [(f_.name, f_.getvalue()) for f_ in _files]
        raw, meta, norm, omnix_msgs = read_exports(_pairs)
        inputs.update({n_: b_ for n_, b_ in _pairs})
    if up.get("meta") is not None:            # separate metadata file wins over one inside the zips
        up["meta"].seek(0)
        meta = read_table(up["meta"])
else:
    for k in ("rna", "prot", "met"):
        if up.get(k) is not None:
            up[k].seek(0)
            raw[k] = read_table(up[k])
            up[k].seek(0)
            inputs[up[k].name] = up[k].read()
    if up.get("meta") is not None:
        up["meta"].seek(0)
        meta = read_table(up["meta"])
    if up.get("map") is not None and "met" in raw:
        up["map"].seek(0)
        mp = read_table(up["map"])
        key = mp.columns[0]
        raw["met"] = raw["met"].merge(mp.rename(columns={key: raw["met"].columns[0]}), on=raw["met"].columns[0],
                                      how="left", suffixes=("", "_map"))

if source == UPLOAD and st.session_state.get("up_kind") == "Normalized (log2)":
    norm = {k: {"declared": True, "scale": "log2"} for k in raw}
if meta is not None:              # stray spaces in sample names or group labels must not split a group
    meta = meta.apply(lambda c: c.astype(str).str.strip().where(c.notna(), c) if c.dtype == object else c)
# which layers to integrate: all three, transcriptomics + metabolomics, or proteomics + metabolomics
LAYER_MODES = {"Transcriptomics + Proteomics + Metabolomics": ("rna", "prot", "met"), "Transcriptomics + Metabolomics": ("rna", "met"),
               "Proteomics + Metabolomics": ("prot", "met")}
if DEMO_ONLY and source in (SESSION, OMNIX_UP, UPLOAD):     # Omnix demo access: own data needs a subscription
    raw, meta, norm = {}, None, {}
if source == UPLOAD:        # own data: choose the layers first, then upload only what they need
    mode_opts = list(LAYER_MODES)
else:                       # demos: the combinations their layers allow
    mode_opts = [m_ for m_, ks_ in LAYER_MODES.items() if all(k_ in raw for k_ in ks_)]
if st.session_state.get("layers_mode") not in mode_opts:
    st.session_state.pop("layers_mode", None)
layers_mode = st.session_state.get("layers_mode", mode_opts[0] if mode_opts else None)
need = LAYER_MODES[layers_mode] if layers_mode else ()
missing_layers = [k for k in need if k not in raw]
if layers_mode:             # files of layers not chosen are set aside (kept in their uploader), not used
    raw = {k: v for k, v in raw.items() if k in need}
    norm = {k: v for k, v in norm.items() if k in raw}
data_ok = "met" in raw and len(raw) >= 2 and meta is not None and not missing_layers
n_samples = len(meta) if meta is not None else 0
layer_lbl = {"rna": "mRNA", "prot": "Protein", "met": "Metabolite"}
data_summary = (f"{('Demo · ' + source.split(': ')[-1]) if source.startswith('Demo') else ('Sent from the modules' if source == SESSION else 'Your data')} · "
                f"{len(raw)} layers · {n_samples} samples") if data_ok else (
    ("Needs " + " and ".join(f"a {LAYER_NAMES[k].lower()} table" for k in missing_layers))
    if missing_layers and (raw or meta is not None) else
    "Needs sample metadata" if raw and meta is None else "Load the data tables")

with st.sidebar:
    _setup_head = st.empty()

    # ---- step 1: data
    _open = st.session_state.get("_open_step")         # sticky: a changing "expanded" value would reset the panel
    with st.expander(f"**Data** · {data_summary}", key="step_data", expanded=(not data_ok) or _open == "step_data",
                     icon=":material/check_circle:" if data_ok else ":material/counter_1:",
                     on_change=_accordion, args=("step_data",)):
        source = st.radio("Dataset", SRC_OPTS, key="src", captions=SRC_CAPS, label_visibility="collapsed",
                          on_change=_keep_open, args=("step_data",))
        id_cols = {}
        if source in (SESSION, OMNIX_UP, UPLOAD) and DEMO_ONLY:
            st.caption(":material/lock: Analyzing your own data in Integrix needs an Omnix plan that includes "
                       "Integrix (Multi-Omics, Lab or Enterprise). The two demo datasets are open to everyone.")
        if source == SESSION:
            _ib = st.session_state.get(INBOX_KEY) or {}
            if not _ib:
                st.caption("Nothing has been sent yet. In Transcriptomics, Proteomics, Metabolomics or Lipidomics, "
                           "finish **Normalization**, then press **Send to Integrix** in the side panel. Every "
                           "module's normalized matrix, metadata and processing notes arrive here.")
            else:
                st.caption("Normalized data sent from the Omnix modules in this session. Integrix uses each matrix "
                           "as the module made it (no second log transform or normalization).")
                for k_, v_ in _ib.items():
                    st.checkbox(f"**{v_['platform_name']}** · {v_['n_features']:,} features × {v_['n_samples']} "
                                f"samples · sent {v_['sent']}", key=f"ib_use_{k_}",
                                value=v_.get("default_use", True), on_change=_keep_open, args=("step_data",),
                                disabled=DEMO_ONLY)
                if sum(1 for k_ in _ib if k_ in ("metabolomics", "lipidomics")) == 2:
                    st.caption(":material/info: Metabolomics and Lipidomics both fill the metabolite layer; "
                               "tick the one to integrate (if both are ticked, the first is used).")
                if st.button("Clear the sent data", key="ib_clear", icon=":material/delete:"):
                    st.session_state[INBOX_KEY] = {}
                    st.rerun()
            for m_ in omnix_msgs:
                st.caption(f":material/info: {m_}")
            for k, df in raw.items():
                idc_ = norm.get(k, {}).get("id_columns") or []
                cols_ = list(df.columns)
                id_cols[k] = idc_[0] if idc_ and idc_[0] in cols_ else st.selectbox(
                    f"{LAYER_NAMES[k]} identifier column", cols_, key=f"id_{k}")
        if source == OMNIX_UP:
            st.caption("Upload the normalized export of each Omnix module (Transcriptomics, Proteomics, "
                       "Metabolomics / Lipidomics): one zip per module or one combined zip. Integrix reads each "
                       "export's sidecar and skips the steps the module already did.")
            st.file_uploader("Omnix exports (.zip)", type=["zip", "csv", "tsv", "json"], key="up_omnix",
                             accept_multiple_files=True)
            st.file_uploader("Sample metadata, if not inside the zips", type=["csv", "tsv", "xlsx"], key="up_meta")
            st.download_button("Example of an Omnix export", _example_zip(NORMDEMO),
                               "Integrix_example_Omnix_normalized.zip", mime="application/zip",
                               icon=":material/description:", width="stretch", key="dl_contract",
                               help="The normalized demo in the Omnix export format (tables, metadata and one "
                                    "*_omnix.json sidecar per module), with a README describing the format.")
            for m_ in omnix_msgs:
                st.caption(f":material/info: {m_}")
            for k, df in raw.items():
                idc_ = norm.get(k, {}).get("id_columns") or []
                cols_ = list(df.columns)
                id_cols[k] = idc_[0] if idc_ and idc_[0] in cols_ else st.selectbox(
                    f"{LAYER_NAMES[k]} identifier column", cols_, key=f"id_{k}")
        if source == "Upload my data":
            st.radio("Data type", UP_KINDS, key="up_kind", horizontal=True, on_change=_keep_open, args=("step_data",),
                     help="Raw intensities: Integrix log2-transforms them. Normalized (log2): the tables are used as "
                          "provided (for example the normalized matrices from the Omnix modules); no transform.")
            st.radio("Layers to integrate", mode_opts, key="layers_mode", on_change=_keep_open, args=("step_data",),
                     help="Choose the layers first; only their uploads are shown. Files of a layer you leave out are "
                          "kept, so switching back uses them again.")
            st.caption("Wide tables: one row per feature, one column per sample. The first column is the "
                       "identifier unless you choose another.")
            for k_, lbl_ in (("rna", "Transcriptomics (genes)"), ("prot", "Proteomics (proteins)"),
                             ("met", "Metabolomics / Lipidomics (required)")):
                # an uploader that is not drawn loses its file: unused layers stay drawn but hidden
                with st.container(key=f"upbox_{k_}" + ("" if k_ in need else "_off")):
                    st.file_uploader(lbl_, type=["csv", "tsv", "txt", "xlsx"], key=f"up_{k_}")
            st.file_uploader("Sample metadata (Sample, Group, optional Batch)", type=["csv", "tsv", "xlsx"],
                             key="up_meta")
            st.file_uploader("Optional: metabolite ID table (feature, HMDB / KEGG / ChEBI)",
                             type=["csv", "tsv", "xlsx"], key="up_map")
            for k, df in raw.items():
                id_cols[k] = st.selectbox(f"{LAYER_NAMES[k]} identifier column", list(df.columns), key=f"id_{k}")
        if source == UPLOAD:
            if missing_layers and (raw or meta is not None):
                st.caption(":material/info: Still needed: " + ", ".join(LAYER_NAMES[k] for k in missing_layers))
        elif len(mode_opts) > 1:
            st.radio("Layers to integrate", mode_opts, key="layers_mode", on_change=_keep_open, args=("step_data",),
                     help="Integrate all three layers, or metabolomics with only one enzyme layer (for example to "
                          "see what a two-layer study would find).")
        elif mode_opts:
            st.caption(f"Layers to integrate: **{mode_opts[0]}**")
        if raw:
            st.markdown(" ".join(f"`{layer_lbl[k]}`" for k in raw) + (f" `{n_samples} samples`" if n_samples else ""))
        if norm:
            from trace_core.omnix_io import describe as _desc
            st.markdown(UI.norm_status({LAYER_NAMES[k]: _desc(v) for k, v in norm.items() if k in raw},
                                       [LAYER_NAMES[k] for k in raw if k not in norm]), unsafe_allow_html=True)
        if source.startswith("Demo"):
            st.download_button("Download this example data", _example_zip(source),
                               "Integrix_example_" + "_".join(__import__("re").findall(
                                   r"[A-Za-z0-9]+", source.split(": ")[-1].replace("×", "x"))) + ".zip",
                               mime="application/zip", icon=":material/download:", width="stretch", key="dl_example",
                               help="The three omics tables, the sample metadata, a README with comparisons to try, "
                                    "and the planted ground truth: use it as a template for your own data.")

    # ---- step 2: study design
    sample_col, group_col, batch_col, case, ref = "Sample", "Group", None, None, None
    if meta is not None:
        mc = list(meta.columns)
        _gc0 = "Condition" if "Condition" in mc else ("Group" if "Group" in mc else mc[min(1, len(mc) - 1)])
        _gc = st.session_state.get("des_group", _gc0)
        _groups = sorted(meta[_gc].astype(str).unique()) if _gc in mc else []
        _case = st.session_state.get("des_case", _pick_case(_groups))
        _ref = st.session_state.get("des_ref", _pick_ref(_groups, _case))
        design_summary = f"{_case} over {_ref}" + (f" ({_gc})" if _gc not in ("Group",) else "")
    else:
        design_summary = "Needs sample metadata"
    with st.expander(f"**Study design** · {design_summary}", key="step_design", expanded=_open == "step_design",
                     icon=":material/check_circle:" if meta is not None else ":material/counter_2:",
                     on_change=_accordion, args=("step_design",)):
        if meta is not None:
            mc = list(meta.columns)
            sample_col = st.selectbox("Sample column", mc, index=mc.index("Sample") if "Sample" in mc else 0,
                                      key="des_sample")
            group_col = st.selectbox("Group column", mc, index=mc.index(_gc0),
                                     key="des_group", on_change=_reset_groups,
                                     help="The metadata column that defines the comparison (e.g. Condition, "
                                          "Treatment, Sex or the combined Group).")
            groups = sorted(meta[group_col].astype(str).unique())
            if st.session_state.get("des_case") not in groups:
                st.session_state.pop("des_case", None)
            _long = max((len(g) for g in groups), default=0) > 10     # long labels: one per row, fully readable
            cg1, cg2 = (st, st) if _long else st.columns(2)
            _c0 = _pick_case(groups)
            case = cg1.selectbox("Case", groups, key="des_case", index=groups.index(_c0) if _c0 in groups else 0,
                                 on_change=_keep_open, args=("step_design",))
            _refs = [g for g in groups if g != case]
            if st.session_state.get("des_ref") not in _refs:
                st.session_state.pop("des_ref", None)
            _r0 = _pick_ref(_refs, case)
            ref = cg2.selectbox("Control (reference)", _refs, key="des_ref", index=_refs.index(_r0) if _r0 in _refs else 0,
                                on_change=_keep_open, args=("step_design",))
            if len(groups) > 2:
                st.caption(f"{len(groups)} groups in '{group_col}': only the two chosen are compared.")
            bopts = ["(none)"] + [c for c in mc if c not in (sample_col, group_col)]
            batch_col = st.selectbox("Batch column", bopts, key="des_batch",
                                     index=bopts.index("Batch") if "Batch" in bopts else 0)
            batch_col = None if batch_col == "(none)" else batch_col
            bc = {}
            for k in ("rna", "prot", "met"):
                _done = bool(norm.get(k, {}).get("batch_corrected"))
                bc[k] = st.checkbox(f"Batch-correct {LAYER_NAMES[k]}" + (" (done by the module)" if _done else ""),
                                    value=False, key=f"bc_{k}",
                                    disabled=batch_col is None or _done)
                if _done:
                    bc[k] = False
        else:
            st.caption("Upload a sample metadata table in step 1.")
            bc = {k: False for k in ("rna", "prot", "met")}

    # ---- step 3: statistics
    _B = st.session_state.get("st_B", 1000)
    _a = st.session_state.get("st_alpha", 0.05)
    _pt = st.session_state.get("st_pthr", 1.5)
    with st.expander(f"**Statistics** · {_B:,} perms · q {_a:g}", key="step_stats", expanded=_open == "step_stats",
                     icon=":material/check_circle:", on_change=_accordion, args=("step_stats",)):
        B = st.select_slider("Label permutations", [200, 500, 1000, 2000, 5000], value=1000, key="st_B",
                             on_change=_keep_open, args=("step_stats",))
        c1_, c2_ = st.columns(2)
        alpha = c1_.selectbox("FDR (q)", [0.01, 0.05, 0.1], index=1, key="st_alpha",
                              on_change=_keep_open, args=("step_stats",))
        pthr = c2_.selectbox("Pattern |z|", [1.0, 1.5, 2.0, 2.5, 3.0], index=1, key="st_pthr")
        with st.popover("More options", icon=":material/tune:", width="stretch"):
            Bt = st.select_slider("Topology permutations", [200, 500, 1000, 2000], value=1000, key="st_Bt")
            moderated = st.checkbox("Moderated t (empirical Bayes)", True, key="st_mod")
            seed = st.number_input("Random seed", value=20261030, step=1, key="st_seed")
            _all_norm = bool(raw) and all(k in norm for k in raw)
            transform = st.selectbox("Log transform", ["auto", "log2", "none"], key="st_tf", disabled=_all_norm,
                                     help="Layers that are already normalized (Omnix exports) are never "
                                          "transformed again. auto: log2 if values look like raw intensities.")
            if _all_norm:
                st.caption("Skipped: every layer is already normalized on log2 scale.")
            max_missing = st.slider("Max missing per group", 0.1, 0.9, 0.5, 0.05, key="st_miss")
            imp_lbl = st.selectbox("Missing values", ["Minimum − 1 (below detection limit)",
                                                      "k-nearest neighbours (missing at random)"], key="st_imp",
                                   index=1,
                                   help="Minimum − 1 suits values missing because they were too low to detect "
                                        "(typical for metabolomics). kNN suits values missing at random (common "
                                        "in DIA proteomics).")

    # ---- step 4: engines
    _eng = [n for n, k in (("Network", "en_net"), ("MOFA+", "en_fac"), ("Signature", "en_sig"))
            if st.session_state.get(k, True)]
    with st.expander(f"**Engines** · {', '.join(_eng) or 'cascade only'}", key="step_engines", expanded=False,
                     icon=":material/check_circle:", on_change=_accordion, args=("step_engines",)):
        do_net = st.toggle("Network, propagation and targets", True, key="en_net")
        do_fac = st.toggle("Latent factors (MOFA+)", True, key="en_fac")
        n_fac = st.slider("Maximum factors", 3, 15, 8, key="en_nfac", disabled=not do_fac)
        do_sig = st.toggle("Integrated biomarker signature", True, key="en_sig")
        sig_rep = st.slider("Outer CV repeats", 3, 20, 10, key="en_rep", disabled=not do_sig)
        sig_perm = st.select_slider("Signature permutations", [20, 50, 100, 200], value=50, key="en_perm",
                                    disabled=not do_sig)
        drug = st.file_uploader("Optional: druggable genes (one symbol per line)", type=["txt", "csv"], key="up_drug")
        st.session_state["fig_fmt"] = st.selectbox("Figure download format (camera icon)", ["png", "svg"], key="fg_fmt")
        st.session_state["fig_scale"] = st.select_slider("PNG scale (× screen size)", [4, 6, 8, 10], value=6,
                                                         key="fg_scale")

    ok = data_ok and case is not None and ref is not None
    n_ready = int(data_ok) + int(ok) + 2
    _setup_head.markdown(UI.setup_head(n_ready, 4), unsafe_allow_html=True)
    st.write("")
    run = st.button("Run analysis", type="primary", icon=":material/play_arrow:", disabled=not ok,
                    width="stretch", key="run_btn")
    st.caption(f"This takes up to a minute.  \nKnowledge base {kb.version} · Human-GEM {kb.manifest['sources']['Human-GEM'].get('version', '')}  \n"
               "Light / dark: ⋮ menu → Settings")

# ===================================================================================== top bar
_top = st.container(key="topbar")


def _go(section: int, tab: int | None = None):
    st.session_state["nav_sec"] = section
    if tab is not None:
        st.session_state[f"tabs_{section}"] = TAB_NAMES[tab]
        st.session_state["_nav_nonce"] = st.session_state.get("_nav_nonce", 0) + 1


TAB_NAMES = [":material/space_dashboard: Overview", ":material/table_chart: Data & mapping",
             ":material/compare_arrows: Concordance", ":material/bolt: Reactions", ":material/linear_scale: Cascade",
             ":material/share: Pathways", ":material/hub: Pathway maps & network", ":material/my_location: Targets",
             ":material/tune: Factors", ":material/fingerprint: Signature", ":material/download: Export",
             ":material/menu_book: Methods"]
SECTIONS = [[0, 3, 4, 5], [6, 7], [2, 8, 9], [1, 10, 11]]


SECTION_PAGES = ["Overview · Reactions · Cascade · Pathways", "Pathway maps & network · Targets",
                 "Concordance · Factors · Signature", "Data & mapping · Export · Methods"]


def draw_topbar(done: bool):
    """Option 1: brand header, then (after a run) four section cards and a row with the study chips and Export.
    The current section's pages follow as pill tabs. HTML + st.columns only, so it renders the same on every
    supported Streamlit version; each card is a real button stretched invisibly over the card."""
    tr = st.session_state.get("trace")
    sec = st.session_state.get("nav_sec", 0)
    kb_txt = f"Human-GEM {kb.manifest['sources']['Human-GEM'].get('version', '')}"
    with _top:
        if OMNIX:      # inside Omnix the portal bar shows the brand: a 'Section / Page' note and the page title instead
            if done:
                _title = (st.session_state.get(f"tabs_{sec}") or TAB_NAMES[SECTIONS[sec][0]]).split(": ", 1)[-1]
                _crumb = SECTION_NAMES[sec]
            else:
                _crumb, _title = "Setup", "Analysis setup"
            st.markdown(f"<div class='omx-crumb'>{_crumb} / {_title}</div>", unsafe_allow_html=True)
            st.header(_title)
            return
        with st.container(key="brandbar"):
            st.markdown(UI.topbar(__version__, LOGO, UI.status_pill(tr is not None, tr.get("seconds") if tr else None),
                                  kb_txt, with_nav=False), unsafe_allow_html=True)
        if not done or OMNIX:          # inside Omnix the sections are listed in the side panel
            return
        with st.container(key="seccards"):
            cols = st.columns(4, gap="small")
            for i, nm in enumerate(SECTION_NAMES):
                with cols[i].container(key=f"seccard_{i}"):
                    st.markdown(UI.section_card(nm, SECTION_PAGES[i], i, i == sec), unsafe_allow_html=True)
                    st.button(f"Open {nm}", key=f"secbtn_{i}", on_click=_go, args=(i,), width="stretch")


if not ok:
    draw_topbar(False)
    st.markdown('<div class="trace-warn">Load a metabolomics layer, at least one of transcriptomics or proteomics, '
                'and a sample metadata table (step 1 in the setup panel), or choose a demo dataset.</div>',
                unsafe_allow_html=True)
    _stop()

# ===================================================================================== run
@st.cache_resource(show_spinner=False)
def _result_store() -> dict:
    """One store for the whole server: finished analyses keyed by a hash of the input data and every setting."""
    return {}


def _run_key(raw: dict, meta, *settings) -> str:
    import hashlib
    h = hashlib.sha256()
    for k in sorted(raw):
        df = raw[k]
        h.update(k.encode())
        h.update(repr((list(map(str, df.columns)), df.shape)).encode())
        h.update(pd.util.hash_pandas_object(df.astype(str), index=True).values.tobytes())
    if meta is not None:
        h.update(pd.util.hash_pandas_object(meta.astype(str), index=True).values.tobytes())
    h.update(repr(settings).encode())
    h.update(str(kb.version).encode())
    return h.hexdigest()


def _analyse(raw: dict, meta: pd.DataFrame, id_cols: dict, transform: str, max_missing: float, bc: dict,
             sample_col: str, group_col: str, batch_col, case: str, ref: str, params: dict, druggable: tuple,
             do_net: bool, do_fac: bool, n_fac: int, do_sig: bool, sig_rep: int, sig_perm: int, _prog=None):
    """The full Integrix analysis. Results are kept for the whole server (all sessions) under a hash of the input data
    and every setting, so a repeated run - the demo, a page reload, a second user - returns at once.
    (Not an st.cache_* function: those replay the progress-bar calls and fail across sessions.)"""
    import threading
    store = _result_store()
    key = _run_key(raw, meta, id_cols, transform, max_missing, bc, sample_col, group_col, batch_col, case, ref,
                   params, druggable, do_net, do_fac, n_fac, do_sig, sig_rep, sig_perm)
    if key in store:
        return store[key]
    pending = store.setdefault("_pending", {})
    with store.setdefault("_lock", threading.Lock()):
        ev = pending.get(key)
        mine = ev is None
        if mine:
            ev = pending[key] = threading.Event()
    if not mine:                             # the same analysis is already running (e.g. the background demo)
        t_w = time.time()
        while not ev.wait(0.5):
            if _prog is not None:
                _prog.progress(min(0.95, 0.1 + (time.time() - t_w) / 60), text="Finishing the analysis already in progress")
        if key in store:
            return store[key]
        with store["_lock"]:                 # it failed: compute here
            pending[key] = ev = threading.Event()
    try:
        out = _compute(raw, meta, id_cols, transform, max_missing, bc, sample_col, group_col, batch_col, case, ref,
                       params, druggable, do_net, do_fac, n_fac, do_sig, sig_rep, sig_perm, _prog)
        results = [k for k in store if not k.startswith("_")]
        if len(results) >= (3 if OMNIX else 8):   # keep memory bounded (less inside the shared Omnix server)
            store.pop(results[0])
        store[key] = out
        return out
    finally:
        ev.set()
        pending.pop(key, None)


def _compute(raw, meta, id_cols, transform, max_missing, bc, sample_col, group_col, batch_col, case, ref, params,
             druggable, do_net, do_fac, n_fac, do_sig, sig_rep, sig_perm, _prog=None):
    say = (lambda f, m: _prog.progress(f, text=m)) if _prog is not None else (lambda f, m: None)
    say(0.02, "Preparing layers")
    layers = {}
    for k, df in raw.items():
        layers[k] = make_layer(k, df, meta, kb, id_col=id_cols.get(k),
                               transform="none" if k in params.get("normalized_input", {}) else transform,
                               max_missing=max_missing, batch_correct=bc[k], sample_col=sample_col,
                               impute=params.get("impute", "min"),
                               group_col=group_col, batch_col=batch_col)
    design = build_design(meta, layers, case, ref, sample_col, group_col, batch_col)
    res = run_trace(layers, design, kb, params, progress=lambda f, m: say(min(0.6, 0.6 * f), f"Engine A · {m}"))
    extras = {"kb_version": kb.version}
    seed = params["seed"]
    if do_net:
        from trace_core import network as N
        say(0.65, "Engine C · network propagation and targets")
        drug_set = set(druggable) if druggable else None
        if drug_set:
            drug_set = {g for g in {x for gs in res.cascade.rx_genes.values() for x in gs} if g.upper() in drug_set}
        nodes, edges = N.build_graph(res, kb)
        prop = N.propagate(nodes, edges, B=200, seed=int(seed))
        extras["targets"] = N.rank_targets(res, kb, prop, drug_set)
        extras["subnet"] = N.subnetwork(res, prop, edges)
        extras["prop"] = prop
        extras["network_params"] = {"restart": 0.3, "null": "degree-binned seed permutation, B=200"}
    if do_fac:
        from trace_core import factors as FA
        say(0.75, "Engine B · MOFA+ latent factors")
        try:
            extras["factors"] = FA.run_factors(layers, design, n_factors=n_fac, seed=int(seed) % 10000)
            extras["factor_enrichment"] = FA.factor_enrichment(extras["factors"], layers, kb)
            extras["factor_params"] = {"max_factors": n_fac, "top_var": 1500, "restarts": 3}
        except Exception as e:  # keep the rest of the results
            extras["factors_error"] = str(e)
    if do_sig:
        from trace_core import signature as SG
        say(0.85, "Engine D · nested cross-validated signature")
        extras["signature"] = SG.run_signature(layers, design, n_repeats=sig_rep, n_perm=sig_perm,
                                               seed=int(seed) % 10000)
        extras["signature_params"] = {"keepX_grid": [5, 10, 20], "repeats": sig_rep, "permutations": sig_perm}
    return res, layers, extras


def _run_args():
    druggable_ = ()
    if drug is not None:
        druggable_ = tuple(sorted({x.strip().upper() for x in drug.getvalue().decode(errors="ignore")
                                   .replace(",", "\n").split("\n") if x.strip()}))
    params_ = dict(B=int(B), B_topo=int(Bt), alpha=float(alpha), pattern_thr=float(pthr), moderated=moderated,
                   seed=int(seed), impute="knn" if imp_lbl.startswith("k") else "min",
                   max_missing=float(max_missing), batch_corrected=[k for k in raw if bc.get(k) and batch_col],
                   normalized_input={k: v for k, v in norm.items() if k in raw}, layers_present=sorted(raw))
    return (raw, meta, dict(id_cols), transform, float(max_missing), dict(bc), sample_col, group_col, batch_col,
            case, ref, params_, druggable_, bool(do_net), bool(do_fac), int(n_fac), bool(do_sig), int(sig_rep),
            int(sig_perm))


# the progress bar's slot is drawn on every load (empty when nothing runs): if it existed only during a run, the
# next click would shift the page tabs and older Streamlit would reopen the first page (Overview)
_prog_slot = st.empty()
if run:
    t0 = time.time()
    prog = _prog_slot.progress(0.0, text="Starting the analysis")
    try:
        res, layers, extras = _analyse(*_run_args(), _prog=prog)
    except Exception as e:  # a problem in the uploaded data: explain it instead of showing a Python traceback
        prog.empty()
        draw_topbar(False)
        msg = str(e) if isinstance(e, (ValueError, KeyError)) and len(str(e)) < 400 else \
            "The analysis stopped because of a problem in the input data."
        st.error(f"**The analysis could not run.** {msg.strip(chr(34) + chr(39))}", icon=":material/error:")
        st.markdown("Common causes: sample names in the data tables do not match the metadata `Sample` column; "
                    "the `Group` column has fewer than 2 samples per group; the identifier column is wrong; "
                    "a table has no numeric sample columns.")
        with st.expander("Technical details (for support)"):
            import traceback
            st.code("".join(traceback.format_exception(type(e), e, e.__traceback__))[-3000:], language=None)
        st.session_state.pop("trace", None)
        _stop()
    st.session_state["trace"] = {"res": res, "layers": layers, "extras": extras, "inputs": inputs,
                                 "run_id": time.time(), "seconds": time.time() - t0}
    prog.empty()

if "trace" not in st.session_state:
    draw_topbar(False)
    _own = ("" if (OMNIX and DEMO_ONLY) else
            " To analyze <b>your own data</b>, open step 1 <b>Data</b> and choose <b>From this Omnix session</b> "
            "(data sent with <b>Send to Integrix</b> from Transcriptomics, Proteomics, Metabolomics or Lipidomics), "
            "<b>Upload Omnix exports</b> or <b>Upload my data</b>.")
    st.markdown('<div class="trace-note">Check the four setup steps in the <b>Analysis setup</b> panel, then press '
                '<b>Run analysis</b>. The demo takes about 20-40 seconds.' + _own + ' On a small screen, open the '
                'panel with the <b>&raquo;</b> button at the top left.</div>', unsafe_allow_html=True)
    _stop()

draw_topbar(True)
S = st.session_state["trace"]
res, layers, extras = S["res"], S["layers"], S["extras"]
rx, pw = res.reactions, res.pathways
lname = {"rna": "mRNA", "prot": "protein", "met": "metabolite"}
from trace_core import ipa as IPA  # noqa: E402
from trace_core import netview as NV  # noqa: E402
from trace_core import pubfig as PF  # noqa: E402
from trace_core.layout import LAYOUTS  # noqa: E402

if "canonical" not in S:
    S["canonical"] = IPA.canonical_pathways(res, kb, res.params["alpha"])
cp_tbl = S["canonical"]

SEC = st.session_state.get("nav_sec", 0)
SEC_IDX = SECTIONS[SEC]
_tkey = f"tabs_{SEC}"
if OMNIX:          # the result pages in the side panel, under their section (navigation only)
    _here = st.session_state.get(_tkey) or TAB_NAMES[SEC_IDX[0]]
    with st.sidebar.container(key="omx_navtree"):
        for _si, _nm in enumerate(SECTION_NAMES):
            st.markdown(f"<div class='omx-navgrp'>{_nm}</div>", unsafe_allow_html=True)
            for _ti in SECTIONS[_si]:
                st.button(TAB_NAMES[_ti].split(": ", 1)[-1], key=f"navtree__sec{_si}__tab{_ti}", on_click=_go,
                          args=(_si, _ti), type="primary" if (_si == SEC and TAB_NAMES[_ti] == _here) else "secondary",
                          use_container_width=True)
with st.container(key="pagebar"):
    _pb = st.columns([5, 1], vertical_alignment="center")
    _pb[0].markdown(UI.study_line(res.design.tier_text, res.design.case, res.design.ref, res.params["alpha"]),
                    unsafe_allow_html=True)
    _pb[1].button("Export", icon=":material/download:", key="top_export", on_click=_go, args=(3, 10),
                  width="stretch")
st.markdown(UI.layer_note(res), unsafe_allow_html=True)
_gm = res.design.meta["Group"].value_counts()
_nc, _nr = int(_gm.get(res.design.case, 0)), int(_gm.get(res.design.ref, 0))
import math as _math  # noqa: E402
_nperm = _math.comb(_nc + _nr, _nc) if _nc and _nr else 0
if 0 < _nperm < 200:
    st.markdown(f'<div class="trace-warn"><b>Small study ({_nc} vs {_nr} samples).</b> Only {_nperm:,} distinct '
                f'label permutations exist, so the permutation null is coarse and few reactions can reach q &lt; '
                f'{res.params["alpha"]:g}. Treat results as exploratory; more samples per group give real power.</div>',
                unsafe_allow_html=True)
try:  # Streamlit >= 1.50: remember the active tab and render only that tab (much faster)
    with st.container(key="ix_pages_lazy"):
        _sub = st.tabs([TAB_NAMES[i] for i in SEC_IDX], key=_tkey, on_change="rerun")
    LAZY = True
except TypeError:              # older Streamlit: every tab renders; open the one a button asked for (Export)
    # after a button asked for a page (Export), the labels change invisibly (zero-width space) so the browser
    # builds the tabs anew and opens that page; ordinary reruns keep the page the user is on
    _zw = "\u200b" * (st.session_state.get("_nav_nonce", 0) % 2)
    _labels = [TAB_NAMES[i] + _zw for i in SEC_IDX]
    _want = st.session_state.get(_tkey)
    _want = _want + _zw if _want else None
    with st.container(key="ix_pages"):
        try:
            _sub = st.tabs(_labels, default=_want if _want in _labels else None)
        except TypeError:
            _sub = st.tabs(_labels)
    LAZY = False
if OMNIX:      # the pages are listed in the side panel (and reached with Next): no second row of page tabs
    st.markdown("<style>[class*='st-key-ix_pages'] [data-testid='stTabs'] [role='tablist']{display:none !important}"
                "[class*='st-key-ix_pages'] [role='tabpanel'] [role='tablist']{display:flex !important}</style>",
                unsafe_allow_html=True)
_hidden = st.container()          # tabs of other sections render nothing
tabs = [(_sub[SEC_IDX.index(i)] if i in SEC_IDX else _hidden) for i in range(len(TAB_NAMES))]


def on(i: int) -> bool:
    if i not in SEC_IDX:
        return False
    if not LAZY:
        return True
    return st.session_state.get(_tkey, TAB_NAMES[SEC_IDX[0]]) in (TAB_NAMES[i], None)


# ------------------------------------------------------------------------------------------- network UI
def _preset_cb(prefix):
    name = st.session_state.get(f"{prefix}_preset")
    for k, v in NV.PRESETS.get(name, {}).items():
        st.session_state[f"{prefix}_{k}"] = v


SCOPE_LABEL = {"Pathway": "Pathway map (IPA-style)", "Coherent reactions": "Coherent reactions",
               "Enzyme neighbourhood": "Enzyme / metabolite neighbourhood", "Molecule list": "My molecule list",
               "Top targets": "Top prioritised targets"}


def _img(data):
    try:
        st.image(data, width="stretch")
    except TypeError:  # older Streamlit
        st.image(data, use_container_width=True)


def network_ui(prefix: str):
    """One studio for IPA-style pathway maps and every other network view. Scope controls sit on top (they
    apply immediately); style settings live in a collapsible, tabbed panel whose two-column grids stack on
    narrow screens; the figure always spans the full width and scales with the window."""
    d = NV.NetOptions()
    K = lambda f: f"{prefix}_{f}"
    pchoices = pw.loc[pw["n_reactions_observed"] > 0].sort_values("p_combined")["pathway"].tolist()
    if K("scope") not in st.session_state:
        st.session_state[K("scope")] = "Pathway"
    if st.session_state.get(K("pathway")) not in pchoices:
        st.session_state[K("pathway")] = pchoices[0] if pchoices else ""
    for f, v in vars(d).items():   # initialise widget state once
        if K(f) not in st.session_state and not isinstance(v, (list, tuple)):
            st.session_state[K(f)] = v
    c1, c2 = st.columns([1, 1.6])
    c1.selectbox("Draw", list(SCOPE_LABEL), format_func=SCOPE_LABEL.get, key=K("scope"))
    scope = st.session_state[K("scope")]
    with c2:
        if scope == "Pathway":
            st.selectbox("Pathway (most significant first)", pchoices, key=K("pathway"))
        elif scope in ("Enzyme neighbourhood", "Molecule list"):
            cc = st.columns([3, 1])
            cc[0].text_input("Genes / metabolites (comma-separated)", key=K("seeds_txt"))
            if scope == "Enzyme neighbourhood":
                cc[1].number_input("Hops", 1, 2, key=K("hops"))
        elif scope == "Coherent reactions":
            st.caption("Reactions whose enzyme and metabolite changes agree (label and topology q below the FDR).")
        else:
            st.caption("The highest-ranked targets from Engine C and the reactions around them.")
    st.selectbox("Style preset (sets layout, shapes, colors and edges at once; fine-tune below)",
                 ["(custom)"] + list(NV.PRESETS), key=K("preset"), on_change=_preset_cb, args=(prefix,))
    with st.expander("Layout, colors, labels and edges", expanded=False):
        with st.form(f"{prefix}_form", border=False):
            t = st.tabs(["Content", "Layout", "Nodes & colors", "Labels", "Edges", "Canvas"])
            with t[0]:
                c = st.columns(2)
                c[0].slider("Maximum reactions", 5, 120, key=K("max_reactions"))
                c[1].slider("Isozymes shown per reaction", 1, 8, key=K("isozymes_per_reaction"))
                c[0].selectbox("Representation", ["Molecule view (IPA-style)", "Reaction view (bipartite)"], key=K("representation"))
                c[1].selectbox("Gene value", ["Best available (protein > mRNA)", "Protein", "mRNA"], key=K("value"))
                c[0].checkbox("Show unmeasured molecules", key=K("include_unmeasured"))
                c[1].checkbox("Hide isolated molecules", key=K("hide_isolated"))
                c[0].checkbox("Add protein-protein edges (STRING add-on)", key=K("include_ppi"),
                              disabled="ppi_string" not in kb.extra)
            with t[1]:
                c = st.columns(2)
                c[0].selectbox("Layout", LAYOUTS, key=K("layout"))
                c[1].selectbox("Subcellular bands", ["IPA (5 bands)", "Organelles (10 bands)"], key=K("band_detail"))
                c[0].slider("Spacing", 0.8, 2.5, key=K("spacing"), step=0.1)
                c[1].text_input("Focus molecules (radial layout)", "", key=K("focus_txt"))
                c[0].number_input("Layout seed", 1, 9999, key=K("seed"))
                c[1].checkbox("Auto-enlarge canvas so nothing overlaps", key=K("auto_size"))
                c[1].checkbox("Shade network modules (communities)", key=K("hulls"))
            with t[2]:
                c = st.columns(2)
                c[0].selectbox("Color by", NV.COLOR_BY, key=K("color_by"))
                c[1].selectbox("Shape scheme", NV.SHAPES, key=K("shape_scheme"))
                c[0].selectbox("Palette", list(NV.PALETTES), key=K("palette"))
                c[1].slider("Color scale |z| max", 1.0, 8.0, key=K("zmax"), step=0.5)
                c[0].selectbox("Node size by", ["Constant", "|z|", "Degree", "Target score"], key=K("size_mode"))
                c[1].slider("Node size (pt)", 6.0, 36.0, key=K("node_size"), step=1.0)
                c[0].checkbox("Predict unmeasured metabolites (MAP-like)", key=K("prediction"))
                c[1].slider("Prediction threshold |E|", 0.5, 3.0, key=K("pred_thr"), step=0.1)
                c[0].text_input("Highlight (comma-separated)", "", key=K("highlight_txt"))
                c[1].color_picker("Highlight outline", key=K("highlight_color"))
            with t[3]:
                c = st.columns(2)
                c[0].selectbox("Label", ["All", "Measured", "Significant", "Top N by |z|", "None"], key=K("label_mode"))
                c[1].slider("Top N", 5, 100, key=K("label_top_n"))
                c[0].slider("Label size (pt)", 5.0, 16.0, key=K("font_size"), step=0.5)
                c[1].selectbox("Font family", _fonts(), key=K("font_family"))
                c[0].checkbox("Label background", key=K("label_box"))
                c[1].checkbox("Gene names inside boxes (KEGG scheme)", key=K("label_inside"))
                c[0].checkbox("Hide labels that cannot be placed cleanly", key=K("hide_unplaced"))
            with t[4]:
                c = st.columns(2)
                c[0].selectbox("Edge color", ["IPA relationship", "By relation type", "Single color"], key=K("edge_color"))
                c[1].color_picker("Single color", key=K("edge_single"))
                c[0].selectbox("Routing", NV.ROUTES, key=K("edge_route"))
                c[1].selectbox("Width by", ["Constant", "|enzyme evidence|"], key=K("edge_width_by"))
                c[0].slider("Edge width", 0.3, 4.0, key=K("edge_width"), step=0.1)
                c[1].slider("Edge opacity", 0.1, 1.0, key=K("edge_alpha"), step=0.05)
                c[0].checkbox("Arrow heads", key=K("arrows"))
                c[1].checkbox("Line style encodes effect (print)", key=K("effect_linestyle"))
            with t[5]:
                c = st.columns(2)
                c[0].text_input("Title", key=K("title"))
                c[1].selectbox("Theme", ["Light", "Dark"], key=K("theme"))
                c[0].slider("Width (in)", 5.0, 30.0, key=K("width_in"), step=0.5)
                c[1].slider("Height (in)", 4.0, 30.0, key=K("height_in"), step=0.5)
                c[0].checkbox("Legend", key=K("show_legend"))
                c[1].checkbox("Draw compartment bands", key=K("show_bands"))
                c[0].selectbox("Background", ["white", "transparent"], key=K("background"))
            st.form_submit_button("Apply", type="primary")
    vals = {f: st.session_state.get(K(f), v) for f, v in vars(d).items() if not isinstance(v, (list, tuple))}
    split = lambda t_: [x.strip() for x in str(st.session_state.get(K(t_), "")).split(",") if x.strip()]
    o = NV.NetOptions(**vals)
    o.seeds = split("seeds_txt") or ["ASS1"]
    o.highlight = split("highlight_txt")
    o.focus = split("focus_txt")
    o.size_range = d.size_range
    if o.scope == "Pathway" and not o.title:
        o.title = o.pathway
    ck = (S.get("run_id"), repr(o))
    cache = S.setdefault("net_cache", {})
    if ck not in cache:
        with st.spinner("Laying out the network"):
            N, E = NV.build_view(res, kb, o, extras.get("targets"))
            if N.empty:
                cache[ck] = (N, E, None, None)
            else:
                g = NV.geometry(N, E, o)
                import matplotlib.pyplot as plt
                buf = io.BytesIO()
                with plt.rc_context({"svg.fonttype": "path"}):
                    NV.render_mpl(N, E, o, g).savefig(buf, format="svg", bbox_inches="tight", pad_inches=0.05)
                plt.close("all")
                cache[ck] = (N, E, g, buf.getvalue())
        if len(cache) > 12:
            cache.pop(next(iter(cache)))
    N, E, g, png = cache[ck]
    if g is None:
        st.info("Nothing to draw for this selection. Widen the scope or show unmeasured molecules.")
        return
    cl = g.get("collisions", {})
    hard = cl.get("label_label", 0) + cl.get("label_node", 0) + cl.get("edge_node", 0)
    if not cl.get("checked", True):
        status = "labels placed radially outside the layout"
    elif hard == 0:
        status = "no overlapping labels or nodes, no edges through nodes"
    else:
        status = f"{hard} overlap(s) remain: raise Spacing or lower Maximum reactions"
    st.caption(f"{len(N)} molecules, {len(E)} relationships · canvas {g['W']:.1f} × {g['H']:.1f} in · {status}")
    view = st.radio("View", ["Exact preview (identical to the export)", "Interactive (hover and zoom)"],
                    horizontal=True, key=f"{prefix}_view", label_visibility="collapsed")
    if view.startswith("Exact"):
        import base64
        svg_ = UI.fluid_svg(UI.darken_svg(png) if MODE == "dark" else png)
        st.markdown(f'<img src="data:image/svg+xml;base64,{base64.b64encode(svg_).decode()}" '
                    f'style="width:100%;height:auto;display:block"/>', unsafe_allow_html=True)
        st.caption("Vector preview, sharp at any zoom (Ctrl/Cmd + to enlarge). Same layout as the export; "
                   "downloads use the white publication style.")
    else:
        # same rendering as the export (untrimmed, so it maps 1:1 onto the plot), with hover on every molecule
        ik = ("int", ck, MODE)
        if ik not in cache:
            import matplotlib.pyplot as plt
            buf = io.BytesIO()
            with plt.rc_context({"svg.fonttype": "path"}):
                NV.render_mpl(N, E, o, g).savefig(buf, format="svg")
            plt.close("all")
            cache[ik] = UI.darken_svg(buf.getvalue()) if MODE == "dark" else buf.getvalue()
        h_px = int(min(900, max(420, 1150 * g["H"] / g["W"])))
        fig = NV.render_interactive(N, E, o, g, cache[ik], h_px)
        cfg = F.config(st.session_state.get("fig_fmt", "png"), st.session_state.get("fig_scale", 6), f"{prefix}_network")
        cfg.update(scrollZoom=True, displaylogo=False)
        st.plotly_chart(fig, width="stretch", config=cfg, theme=None, key=f"{prefix}_iplot")
        st.caption("Hover a molecule for its details. Drag to zoom into a region, scroll to zoom, double-click to "
                   "reset; the map stays sharp at any zoom and matches the export exactly.")

    def build(w, h, fs_):
        import dataclasses
        o2 = dataclasses.replace(o, width_in=w, height_in=h or o.height_in * w / o.width_in, font_size=fs_)
        return NV.render_mpl(N, E, o2)
    download_panel(f"{prefix}_pub", build, o.width_in, o.height_in, o.font_size, f"{prefix}_network")
    nt, et = NV.export_tables(N, E)
    c = st.columns(4)
    with c[0]:
        dl_card("Cytoscape", f"{prefix}_network", "Open in Cytoscape (.cyjs)", key=f"{prefix}_cycard",
                buttons=[("CYJS", NV.to_cyjs(N, E, g), f"{prefix}_network.cyjs", "application/json", f"{prefix}_cy")])
    with c[1]:
        dl_card("GraphML", f"{prefix}_network", "For Gephi, yEd and other tools", key=f"{prefix}_gmlcard",
                buttons=[("GraphML", NV.to_graphml(N, E), f"{prefix}_network.graphml", "application/xml",
                          f"{prefix}_gml")])
    with c[2]:
        dl_card("Nodes", f"{prefix}_nodes", "One row per molecule", key=f"{prefix}_nodescard",
                buttons=[("CSV", nt, f"{prefix}_nodes.csv", "text/csv", f"{prefix}_nodes")])
    with c[3]:
        dl_card("Edges", f"{prefix}_edges", "One row per link", key=f"{prefix}_edgescard",
                buttons=[("CSV", et, f"{prefix}_edges.csv", "text/csv", f"{prefix}_edges")])


# ===================================================================================== Overview
def _findings():
    """(kind, html) rows for the Key findings card, from the actual results."""
    E = UI._esc
    out = []
    sig = rx[rx["significant"]].drop_duplicates(subset=["enzymes", "E", "M"])
    for _, r in sig.head(8).iterrows():
        kind = ("Bottleneck" if r["pattern"] == "Candidate bottleneck" else
                "Throughput" if r["pattern"].startswith("Throughput") else "Reaction")
        out.append((kind, f"<b>{E(r['enzymes'].split(',')[0])}</b> &middot; {E(r['name'])} &mdash; <i>{E(r['pattern'])}</i>. "
                          f"{E(PATTERN_HELP.get(r['pattern'], ''))} <span class='stat'>E = {r['E']:.1f}, M = {r['M']:.1f}, "
                          f"q = {r['q_label']:.1e}</span>"))
    md = rx[rx["metabolite_driven"] & ~rx["significant"]].sort_values("p_balance")
    for _, r in md.head(3).iterrows():
        out.append(("Metabolite-driven", f"<b>{E(r['name'])}</b> &mdash; metabolite balance changes while the enzyme "
                    f"({E(r['enzymes'].split(',')[0])}) is unchanged. <span class='stat'>M = {r['M']:.1f}, "
                    f"q = {r['q_balance']:.1e}</span>"))
    for _, r in pw[pw["class"].str.startswith("Shared") | pw["class"].str.endswith("-specific")].head(5).iterrows():
        out.append(("Pathway", f"<b>{E(r['pathway'])}</b> &mdash; {E(r['class'])}. <span class='stat'>combined z = "
                    f"{r['z_combined']:.1f}, q = {r['q_combined']:.1e}</span>"))
    if len(cp_tbl):
        for _, r in cp_tbl[cp_tbl["prediction"].isin(["Increased activity", "Decreased activity"])].head(3).iterrows():
            out.append(("Canonical", f"<b>{E(r['pathway'])}</b>: {E(r['prediction'].lower())}. <span class='stat'>"
                        f"directional z = {r['z_directional']:.2f}, p = {r['p']:.1e}</span>"))
    if extras.get("targets") is not None:
        for _, r in extras["targets"].head(3).iterrows():
            out.append(("Target", f"<b>{E(r['gene'])}</b> &mdash; {E(r['best_reaction'])}. "
                        f"{E(r['evidence'] or 'ranked by rank aggregation')}. <span class='stat'>RRA = {r['rra_score']:.1e}</span>"))
    if res.concordance is not None:
        vc = res.concordance["class"].value_counts()
        out.append(("Concordance", f"Gene&ndash;protein concordance: <b>{vc.get('Coupled', 0)}</b> coupled, "
                    f"<b>{vc.get('Protein-only (post-transcriptional)', 0)}</b> protein-only, "
                    f"<b>{vc.get('mRNA-only (buffered)', 0)}</b> mRNA-only, <b>{vc.get('Opposing', 0)}</b> opposing genes."))
    return out


with tabs[0]:
  if on(0):
      shared = pw[pw["class"].str.startswith("Shared")]
      spec = pw[pw["class"].str.endswith("-specific")]
      cpz = res.coupling.get("pearson_z", np.nan) if res.concordance is not None else np.nan
      n_md = int((rx["metabolite_driven"]).sum())
      st.markdown(UI.section("Primary tri-omics metrics", "chart"), unsafe_allow_html=True)
      st.markdown(UI.metric_cards([
          dict(label="Observable reactions", value=f"{len(rx):,}", icon="db", accent="var(--faint)",
               sub="measured enzyme + metabolite"),
          dict(label="Coherent reactions", value=int(rx["significant"].sum()), icon="check", accent="var(--green)",
               highlight=True, sub=f"<b>label &amp; topology</b> q &lt; {res.params['alpha']:g}"),
          dict(label="Metabolite-driven", value=n_md, icon="flask", accent="var(--amber)", highlight=True,
               sub="<b>balance changes</b>, enzyme unchanged"),
          dict(label="Shared pathways", value=len(shared), icon="branch", accent="var(--purple)",
               sub=UI._esc(shared["pathway"].iloc[0][:34]) if len(shared) else "none at this FDR"),
          dict(label="Omics-specific", value=len(spec), icon="stack", accent="var(--sky)", sub="single-layer pathways"),
          dict(label="mRNA–protein r", value=f"{cpz:.2f}" if np.isfinite(cpz) else "n/a", icon="trend",
               accent="var(--brand)", sub=f"{res.coupling.get('n_genes', 0):,} genes in both layers"
               if res.concordance is not None else "needs mRNA and protein"),
      ]), unsafe_allow_html=True)
      g1, g2 = st.columns([1.85, 1])
      with g1.container(border=True):
          st.markdown(UI.card_head("Layer concordance across top pathways",
                                   "Per-layer pathway activity z (permutation-standardized Stouffer) for the six "
                                   "most significant pathways; dashed lines at |z| = 2.", "trend"), unsafe_allow_html=True)
          show(F.layer_bars(pw, 6), "layer_concordance", width_px=720)
      with g2.container(border=True):
          st.markdown(UI.card_head("Reaction classification", "All observable reactions by evidence category.",
                                   "target"), unsafe_allow_html=True)
          show(F.reaction_donut(rx), "reaction_classification", width_px=400)
      with st.container(border=True):
          h1, h2 = st.columns([1.5, 1.6])
          h1.markdown(UI.card_head("Key findings & biological interpretations",
                                   "Top-scoring reactions, pathways and targets, categorised by Integrix.", "spark"),
                      unsafe_allow_html=True)
          q_ = h2.text_input("Search", placeholder="Search gene, reaction or pathway", label_visibility="collapsed",
                             key="kf_q", icon=":material/search:")
          _seg = st.segmented_control if ST_HAS["segmented_control"] else \
              (lambda label, options, **kw: st.radio(label, options, horizontal=True, key=kw.get("key"),
                                                     label_visibility="collapsed"))
          kind_ = _seg("Filter", ["All", "Bottleneck", "Throughput", "Metabolite-driven", "Pathway",
                                                  "Target"], default="All", key="kf_kind", label_visibility="collapsed")
          rows_ = _findings()
          if kind_ and kind_ != "All":
              rows_ = [r for r in rows_ if r[0] == kind_]
          if q_:
              rows_ = [r for r in rows_ if q_.lower() in r[1].lower()]
          n_all = len(rows_)
          show_all = st.session_state.get("kf_all", False) or bool(q_) or (kind_ not in (None, "All"))
          st.markdown(UI.findings(rows_ if show_all else rows_[:6]), unsafe_allow_html=True)
          if n_all > 6 and not (q_ or (kind_ not in (None, "All"))):
              st.toggle(f"Show all {n_all} findings", key="kf_all")
          st.caption("Results describe evidence consistent with a mechanism, not proof of causation.")
      p1, p2 = st.columns(2)
      with p1.container(border=True):
          st.markdown(UI.card_head("Reaction pattern map", "Enzyme evidence E vs metabolite balance M.", "bolt"),
                      unsafe_allow_html=True)
          show(F.pattern_grid(rx), "reaction_patterns", notitle=True, width_px=620,
               pub=(lambda w, h, fs: PF.pattern_grid(rx, 12, w, h or 4.2, fs), 5.6, 4.2, 7))
      if len(cp_tbl):
          with p2.container(border=True):
              st.markdown(UI.card_head("Canonical pathways", "Right-tailed Fisher's exact test; color = directional z.",
                                       "branch"), unsafe_allow_html=True)
              show(F.canonical_bars(cp_tbl, 15), "canonical_pathways", notitle=True, width_px=620,
                   pub=(lambda w, h, fs: PF.canonical(cp_tbl, 15, w, h, fs), 7.6, None, 7))

# ===================================================================================== Data & mapping
with tabs[1]:
  if on(1):
      st.markdown("### Samples and evidence tier")
      ov = res.design.overlap
      table(pd.DataFrame({"layer": [LAYER_NAMES[k] for k in ov["per_layer"]],
                          "samples in design": list(ov["per_layer"].values())}).assign(
          **{"samples in all layers": ov["all_layers"]}), 200, key="tier_tbl")
      st.caption(res.design.tier_text)
      if "Batch" in res.design.meta.columns:
          v = cramers_v(res.design.meta["Batch"], res.design.meta["Group"])
          msg = "warning: batch and group are confounded" if v > 0.5 else "batch is not confounded with group"
          st.caption(f"Batch vs group Cramér's V = {v:.2f} ({msg}).")
      show(F.coverage_bars(layers), "mapping_coverage",
           pub=(lambda w, h, fs: PF.coverage(layers, w, h or 2.6, fs), 3.6, 2.6, 7))
      for k, L in layers.items():
          with st.expander(f"{L.name}: processing log and identifier mapping"):
              for line in L.log:
                  st.write("• " + line)
              f = L.feat.reset_index().rename(columns={"index": "feature"})
              table(f.head(500), 260, key=f"map_{k}")
              un = f[f["map_rule"].isin(["unmapped", "not in metabolic KB"])]
              csv_button(un, f"unmapped_{k}.csv", f"Download {len(un)} unmapped features", key=f"un_{k}")

# ===================================================================================== Concordance
with tabs[2]:
  if on(2):
      if res.concordance is None:
          st.info("Concordance needs both transcriptomics and proteomics.")
      else:
          cp = res.coupling
          c = st.columns(4)
          c[0].metric("Genes in both layers", f"{cp['n_genes']:,}")
          c[1].metric("Deming slope (protein on mRNA)", f"{cp['beta']:.2f}")
          c[2].metric("Correlation of z-scores", f"{cp['pearson_z']:.2f}")
          c[3].metric("Median within-sample rho", f"{cp.get('median_within_rho', np.nan):.2f}" if "median_within_rho" in cp else "n/a (unpaired)")
          nl = st.slider("Genes labeled", 0, 40, 15, key="conc_nl")
          show(F.concordance(res.concordance, cp["beta"], nl), "concordance",
               pub=(lambda w, h, fs: PF.concordance(res.concordance, cp["beta"], nl, w, h or 4.0, fs), 5.4, 4.0, 7))
          cls = st.multiselect("Show classes", sorted(res.concordance["class"].unique()),
                               default=[c for c in res.concordance["class"].unique() if c != "Unchanged"])
          d = res.concordance[res.concordance["class"].isin(cls)]
          table(d, 360, key="conc_tbl")
          csv_button(d, "gene_protein_concordance.csv", key="conc_dl")

# ===================================================================================== Reactions
with tabs[3]:
  if on(3):
      st.markdown("Each reaction combines **enzyme evidence E** (protein z through the GPR rule) with the "
                  "**metabolite balance M** (products minus substrates). **Coherence C = E × M.** A reaction is "
                  f"called coherent when both label-permutation and topology q < {res.params['alpha']}.")
      nl = st.slider("Enzymes labeled", 0, 40, 14, key="rx_nl")
      show(F.pattern_grid(rx, nl), "reaction_patterns",
           pub=(lambda w, h, fs: PF.pattern_grid(rx, nl, w, h or 4.2, fs), 6.0, 4.4, 7))
      c1, c2, c3 = st.columns([1.2, 2, 2])
      only = c1.radio("Rows", ["Coherent", "Metabolite-driven", "All"], horizontal=False)
      pats = c2.multiselect("Patterns", list(PATTERN_HELP), default=[])
      q = c3.text_input("Search enzyme, reaction or metabolite")
      d = rx
      if only == "Coherent":
          d = d[d["significant"]]
      elif only == "Metabolite-driven":
          d = d[d["metabolite_driven"]]
      if pats:
          d = d[d["pattern"].isin(pats)]
      if q:
          m = d[["name", "enzymes", "substrates", "products", "subsystem"]].apply(
              lambda s_: s_.astype(str).str.contains(q, case=False)).any(axis=1)
          d = d[m]
      cols = ["name", "enzymes", "subsystem", "pattern", "E", "S", "P", "M", "C", "q_label", "q_topology",
              "q_balance", "enzyme_evidence", "substrates", "products", "ec", "kegg_rxn", "reaction"]
      table(d[cols], 420, key="rx_tbl")
      csv_button(d, "reactions.csv", key="rx_dl")
      with st.expander("What the patterns mean"):
          for k, v in PATTERN_HELP.items():
              st.markdown(f"- **{k}**: {v}")

# ===================================================================================== Cascade
with tabs[4]:
  if on(4):
      st.markdown("Follow the chain **enzyme → reaction → metabolite**. Ribbon width = |z|, red = up, blue = down; "
                  "pale ribbons and italic names are substrates. Full reaction names are shown, each in its own box.")
      mode = st.radio("Select reactions by", ["Pathway", "Top coherent reactions", "Enzyme"], horizontal=True)
      if mode == "Pathway":
          choices = pw.loc[pw["n_reactions_observed"] > 0, "pathway"].tolist()
          p = st.selectbox("Pathway", choices, key="casc_pw")
          sel = rx[rx["subsystem"].str.split("; ").apply(lambda x: p in x)]
          sel = sel.sort_values("p_label").head(st.slider("Reactions shown", 3, 30, 12, key="casc_n"))
      elif mode == "Top coherent reactions":
          sel = rx.sort_values("p_label").head(st.slider("Reactions shown", 3, 30, 12, key="casc_n2"))
      else:
          genes = sorted({g.strip() for e in rx["enzymes"] for g in e.split(",")})
          g = st.selectbox("Enzyme", genes, index=genes.index("ASS1") if "ASS1" in genes else 0)
          sel = rx[rx["enzymes"].str.split(", ").apply(lambda x: g in x)].head(15)
      if len(sel):
          PF.style("DejaVu Sans", 8)
          with st.container(border=True):
              st.markdown(UI.card_head("Enzyme → reaction → metabolite cascade", "Ribbon width = |z|; red = up, "
                                       "blue = down; full reaction names.", "flow"), unsafe_allow_html=True)
              cw_, _ = show_mpl(lambda: PF.cascade_alluvial(sel, w=7.5, fs=8), max_scale=1.3,
                                cache_key=("cascade", S.get("run_id"), tuple(sel["reaction"])))
          download_panel("casc", lambda w, h, fs: PF.cascade_alluvial(sel, w, h, fs), round(cw_ * 2) / 2, None, 8, "cascade")
          table(sel[["name", "enzymes", "pattern", "E", "S", "P", "C", "q_label"]], 260, key="casc_tbl")

# ===================================================================================== Pathways
with tabs[5]:
  if on(5):
      t1, t2 = st.columns([3, 1], vertical_alignment="bottom")
      t1.markdown(UI.page_title("Pathways", "Which pathways change, in which omics layer, and whether the layers "
                                "agree."), unsafe_allow_html=True)
      ntop = t2.selectbox("Pathways shown", [10, 15, 25, 40, 60], index=2, key="pw_n",
                          format_func=lambda n: f"Top {n}")
      c1, c2 = st.columns(2, gap="medium")
      with c1.container(border=True):
          st.markdown(UI.card_head("Canonical pathways (IPA-style)",
                                   "Bar = −log10 p (right-tailed Fisher's exact) · color = directional z "
                                   "(orange up, blue down, grey &lt; 4 directional molecules) · line = ratio",
                                   "branch"), unsafe_allow_html=True)
          if len(cp_tbl):
              show(F.canonical_bars(cp_tbl, ntop), "canonical_pathways", notitle=True, width_px=620,
                   pub=(lambda w, h, fs: PF.canonical(cp_tbl, ntop, w, h, fs), 7.6, None, 7))
          else:
              st.info("No pathway has enough measured molecules.")
      with c2.container(border=True):
          st.markdown(UI.card_head("Pathway activity across layers",
                                   "Color = permutation-standardized z · size = −log10 q · tile = class "
                                   "(shared-concordant, shared-discordant or omics-specific)", "chart"),
                      unsafe_allow_html=True)
          show(F.pathway_dotplot(pw, ntop), "pathway_activity", notitle=True, width_px=620,
               pub=(lambda w, h, fs: PF.pathway_dotplot(pw, ntop, w, h, fs), 7.2, None, 7))
      with st.container(border=True):
          h1, h2 = st.columns([2, 1.4], vertical_alignment="center")
          h1.markdown(UI.card_head("Pathway table", "Per-layer z and q, combined z, heterogeneity (Cochran's Q, I²) "
                                   "and class.", "grid"), unsafe_allow_html=True)
          cls = h2.multiselect("Classes", sorted(pw["class"].unique()), default=[], key="pw_cls",
                               placeholder="Filter by class", label_visibility="collapsed")
          tabA, tabB = st.tabs(["Cross-layer activity", "Canonical pathways"])
          with tabA:
              d = pw[pw["class"].isin(cls)] if cls else pw
              table(d, 360, key="pw_tbl")
              csv_button(d, "pathways.csv", key="pw_dl")
          with tabB:
              table(cp_tbl, 360, key="can_tbl")
              csv_button(cp_tbl, "canonical_pathways.csv", key="can_dl")

# ===================================================================================== Pathway maps & network studio
with tabs[6]:
  if on(6):
      st.markdown("**Pathway maps and network studio.** Start from an IPA-style pathway map or any other view. "
                  "**Shape** = molecule type; **red / green** = measured up / down; **orange / blue outline** = "
                  "predicted increase / decrease of unmeasured metabolites; **edge color** = leads to increase "
                  "(orange), decrease (blue), inconsistent with the data (yellow) or not predicted (grey). Edges are "
                  "routed around nodes and every label is placed clear of other labels and nodes.")
      network_ui("net")

# ===================================================================================== Targets
with tabs[7]:
  if on(7):
      if extras.get("targets") is None:
          st.info("Enable 'Network, propagation and targets' in the sidebar and run again.")
      else:
          tg = extras["targets"]
          st.markdown("Targets are ranked by robust rank aggregation of reaction coherence, mRNA-protein agreement, "
                      "choke-point status and network propagation (plus druggability when a list is uploaded). "
                      "**They are hypotheses for experimental follow-up.**")
          nt = st.slider("Targets shown", 5, 50, 20, key="tg_n")
          PF.style("DejaVu Sans", 8)
          with st.container(border=True):
              st.markdown(UI.card_head("Prioritised targets", "Ranked by robust rank aggregation.", "target"),
                          unsafe_allow_html=True)
              st.markdown(UI.targets_list(tg, nt), unsafe_allow_html=True)
          # the export is the same ranked list, drawn as a publication figure (11 in wide by default)
          download_panel("tg", lambda w, h, fs: PF.targets_list(tg, nt, w, h, fs), 11.0, None, 8, "targets")
          table(tg, 320, key="tg_tbl")
          csv_button(tg, "targets.csv", key="tg_dl")
      st.markdown("### Reporter metabolites")
      st.caption("Metabolites around which the surrounding reactions are jointly coherent (Patil & Nielsen 2005). "
                 f"p-values come from a normal fit to {res.params.get('B', 1000)} label permutations; values far below "
                 f"1/{res.params.get('B', 1000)} are extrapolated from that fit and are best read as a ranking.")
      table(res.reporter_mets.head(200), 300, key="rep_tbl")
      if res.regulators is not None:
          st.markdown("### Transcription-factor regulators")
          table(res.regulators, 300, key="tf_tbl")
      else:
          st.caption("TF regulator inference becomes available when CollecTRI regulons are added to the knowledge "
                     "base (python -m trace_core.kb_build --skip-core --collectri).")

# ===================================================================================== Factors
with tabs[8]:
  if on(8):
      f = extras.get("factors")
      if not f:
          st.info(extras.get("factors_error", "Enable 'Latent factors (MOFA+)' in the sidebar and run again."))
      else:
          st.caption(f"Method: {f['method']}. Factors explaining ≥ 2% variance in at least one layer are shown.")
          show(F.variance_heatmap(f["r2"]), "factor_variance",
               pub=(lambda w, h, fs: PF.variance_heatmap(f["r2"], w, h or 1.8, fs), 0.6 * f["r2"].shape[1] + 1.8, 1.8, 7))
          cols = list(f["scores"].columns)
          if len(cols) >= 2:
              a, b = st.columns(2)
              fx = a.selectbox("x", cols, 0)
              fy = b.selectbox("y", cols, 1)
              lbl = st.checkbox("Label samples", False, key="fac_lbl")
              show(F.factor_scatter(f["scores"], res.design.meta["Group"], fx, fy), "factor_scores",
                   pub=(lambda w, h, fs: PF.factor_scatter(f["scores"], res.design.meta["Group"], fx, fy, w,
                                                           h or 3.2, fs, label_samples=lbl), 4.2, 3.4, 7))
          table(f["summary"], 220, key="fac_sum")
          st.markdown("**Factor associations with group and covariates**")
          table(f["assoc"], 240, key="fac_assoc")
          if extras.get("factor_enrichment") is not None and len(extras["factor_enrichment"]):
              st.markdown("**Pathways enriched in factor loadings**")
              table(extras["factor_enrichment"].head(300), 300, key="fac_enr")

# ===================================================================================== Signature
with tabs[9]:
  if on(9):
      s = extras.get("signature")
      if not s:
          st.info("Enable 'Integrated biomarker signature' in the sidebar and run again.")
      else:
          if s.get("status") != "ok":
              st.markdown(f'<div class="trace-warn">{s["status"]}</div>', unsafe_allow_html=True)
          else:
              c = st.columns(4)
              c[0].metric("AUC (nested CV)", f"{s['auc_mean']:.2f}", f"95% CI {s['auc_ci'][0]:.2f}-{s['auc_ci'][1]:.2f}",
                          delta_color="off")
              c[1].metric("Balanced accuracy", f"{s['balanced_accuracy']:.2f}")
              c[2].metric("Permutation p", f"{s['p_perm']:.3f}")
              c[3].metric("Features per block", s["keepX"])
              st.markdown('<div class="trace-warn">Cross-validated performance on one cohort is optimiztic for '
                          'clinical use. Validate the panel in an independent cohort.</div>', unsafe_allow_html=True)
              a, b = st.columns(2)
              with a:
                  show(F.roc_plot(s), "signature_roc", width_px=620,
                       pub=(lambda w, h, fs: PF.roc(s, w, h or w, fs), 3.4, 3.4, 7))
              with b:
                  show(F.perm_hist(s), "signature_permutation", width_px=620,
                       pub=(lambda w, h, fs: PF.perm_hist(s, w, h or 2.6, fs), 3.6, 2.6, 7))
              show(F.panel_bars(s["panel"]), "signature_panel",
                   pub=(lambda w, h, fs: PF.panel_bars(s["panel"], w, h, fs), 5.0, None, 7))
              table(pd.DataFrame({"model": ["Integrated (all layers)"] + [f"{k} only" for k in s["single_layer_auc"]],
                                  "AUC": [s["auc_mean"]] + list(s["single_layer_auc"].values())}).round(3), 200,
                    key="sig_auc_tbl")
              csv_button(s["panel"], "signature_panel.csv", key="sig_dl")
          lf = s.get("late_fusion", {})
          if lf.get("status") == "ok":
              st.markdown("**Late fusion (uses every sample available in each layer)**")
              table(pd.DataFrame({"model": ["Fused"] + list(lf["auc_per_layer"]),
                                  "AUC (CV)": [lf["auc_fused"]] + list(lf["auc_per_layer"].values())}).round(3), 200,
                    key="lf_tbl")

# ===================================================================================== Export
with tabs[10]:
  if on(10):
      st.markdown("### Download results")
      xls = R.excel(res, extras)
      man = R.manifest(res, S["inputs"], kb, extras)
      meth = R.methods_text(res, kb, extras)
      c = st.columns(3)
      with c[0]:
          dl_card("All tables (Excel, one file)", "Integrix_results", "Every results table, one sheet each",
                  key="ix_xls", buttons=[("Excel", xls, "Integrix_results.xlsx", "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", "ix_xls_xlsx")])
      with c[1]:
          dl_card("Run manifest", "Integrix_manifest", "Inputs, settings and versions of this run", key="ix_man",
                  buttons=[("JSON", json.dumps(man, indent=2, default=str), "Integrix_manifest.json",
                            "application/json", "ix_man_json")])
      with c[2]:
          dl_card("Methods text", "Integrix_methods", "Ready to paste into a manuscript", key="ix_meth",
                  buttons=[("Markdown", meth, "Integrix_methods.md", "text/markdown", "ix_meth_md")])
      st.markdown("### All figures at once")
      from trace_core import plotly_export as PX
      _same = PX.available()
      st.caption("Every chart as on screen, typeset for a journal (default column width, 8 pt text), with the "
                 "settings currently chosen on each page." if _same else "Publication-style versions of every chart (figures identical to the screen cannot "
                               "be made on this server right now).")
      c = st.columns(3)
      fmts = c[0].multiselect("Formats", list(PF.FORMATS), default=["PNG", "PDF", "SVG"], key="all_fmts")
      dpi = c[1].selectbox("Resolution (DPI)", [300, 600, 1200], index=1, key="all_dpi")
      font = c[2].selectbox("Font family (charts drawn for print: targets, cascade, pathway maps)", _fonts(),
                            key="all_font")

      def _screen_set():
          """name -> (Plotly figure as on screen, on-screen width in px), with each page's current settings."""
          ss = st.session_state
          figs = {"layer_concordance": (F.layer_bars(pw, 6), 720), "reaction_classification": (F.reaction_donut(rx), 400),
                  "reaction_patterns": (F.pattern_grid(rx, ss.get("rx_nl", 14)), 1100),
                  "pathway_activity": (F.pathway_dotplot(pw, ss.get("pw_n", 25)), 620),
                  "mapping_coverage": (F.coverage_bars(layers), 1100)}
          if len(cp_tbl):
              figs["canonical_pathways"] = (F.canonical_bars(cp_tbl, ss.get("pw_n", 25)), 620)
          if res.concordance is not None:
              figs["concordance"] = (F.concordance(res.concordance, res.coupling["beta"], ss.get("conc_nl", 15)), 1100)
          fx_ = extras.get("factors")
          if fx_:
              figs["factor_variance"] = (F.variance_heatmap(fx_["r2"]), 1100)
              if fx_["scores"].shape[1] >= 2:
                  a_, b_ = list(fx_["scores"].columns[:2])
                  figs["factor_scores"] = (F.factor_scatter(fx_["scores"], res.design.meta["Group"], a_, b_), 1100)
          sg_ = extras.get("signature")
          if sg_ and sg_.get("status") == "ok":
              figs["signature_roc"] = (F.roc_plot(sg_), 620)
              figs["signature_permutation"] = (F.perm_hist(sg_), 620)
              figs["signature_panel"] = (F.panel_bars(sg_["panel"]), 1100)
          return figs

      def _screen_files(exts, dpi_):
          from trace_core import printfig as PRT
          out_ = {}
          for nm_, (fg_, wpx_) in _screen_set().items():
              sf_ = PX.screen_figure(fg_, UI.theme_figure, _blacken())
              tf_, w_in, h_in = PRT.typeset(sf_, nm_)        # journal width, 8 pt text, height from content
              for e_ in exts:
                  out_[f"{nm_}.{e_}"] = PX.figure_bytes(tf_, e_, dpi_, w_in, h_in)
          return out_

      if st.button("Build figure set (ZIP)"):
          with st.spinner("Rendering all figures"):
              PF.style(font, 7)
              builders = {"reaction_patterns": lambda: PF.pattern_grid(rx),
                          "pathway_activity": lambda: PF.pathway_dotplot(pw),
                          "mapping_coverage": lambda: PF.coverage(layers)}
              if len(cp_tbl):
                  builders["canonical_pathways"] = lambda: PF.canonical(cp_tbl)
              if res.concordance is not None:
                  builders["concordance"] = lambda: PF.concordance(res.concordance, res.coupling["beta"])
              if extras.get("targets") is not None:
                  builders["targets"] = lambda: PF.targets_list(extras["targets"])
              if extras.get("factors"):
                  builders["factor_variance"] = lambda: PF.variance_heatmap(extras["factors"]["r2"])
              sg = extras.get("signature")
              if sg and sg.get("status") == "ok":
                  builders["signature_roc"] = lambda: PF.roc(sg)
                  builders["signature_panel"] = lambda: PF.panel_bars(sg["panel"])
              top_pw = pw.loc[pw["n_reactions_observed"] >= 2, "pathway"].head(3).tolist()
              if top_pw:   # enzyme -> reaction -> metabolite cascade of the top pathway (same default as the page)
                  _cs = rx[rx["subsystem"].str.split("; ").apply(lambda x: top_pw[0] in x)].sort_values("p_label").head(12)
                  if len(_cs):
                      builders["cascade"] = lambda _cs=_cs: PF.cascade_alluvial(_cs, w=7.5, fs=7)
              if extras.get("factors") and extras["factors"]["scores"].shape[1] >= 2:
                  _f = extras["factors"]
                  _fx, _fy = list(_f["scores"].columns[:2])
                  builders["factor_scores"] = lambda: PF.factor_scatter(_f["scores"], res.design.meta["Group"], _fx, _fy)
              if sg and sg.get("status") == "ok":
                  builders["signature_permutation"] = lambda: PF.perm_hist(sg)
              for k_, p_ in enumerate(top_pw):
                  o_ = NV.NetOptions(scope="Pathway", pathway=p_, max_reactions=30, title=p_)
                  builders[f"pathway_map_{k_ + 1}"] = (lambda o_=o_: NV.render_mpl(*NV.build_view(res, kb, o_, extras.get("targets")), o_))
              out = {}
              if _same:            # charts that are Plotly on screen: the screen figure itself
                  out.update(_screen_files([PX.FORMATS[fl] for fl in fmts], dpi))
                  done = {k.rsplit(".", 1)[0] for k in out}
                  builders = {k: v for k, v in builders.items() if k not in done}
              for nm, b in builders.items():
                  for fl in fmts:
                      ext = PF.FORMATS[fl]
                      out[f"{nm}.{ext}"] = PF.export(b(), ext, dpi)
              st.session_state["figzip"] = R.bundle(out)
      if "figzip" in st.session_state:
          dl_card("Figure set", "Integrix_figures", "Every figure, ready for a manuscript", key="ix_figzip",
                  buttons=[("ZIP", st.session_state["figzip"], "Integrix_figures.zip", "application/zip", "ix_figzip_zip")])
      st.markdown("### Everything in one archive")
      if st.button("Build full bundle"):
          files = {"Integrix_results.xlsx": xls, "Integrix_manifest.json": man, "Integrix_methods.md": meth}
          figs_ = _screen_files(["png"], 600) if _same else R.pub_figures(res, extras, "png", 600)
          files.update({f"figures/{k}": v for k, v in figs_.items()})
          if len(cp_tbl):
              files["canonical_pathways.csv"] = cp_tbl.to_csv(index=False)
          if extras.get("subnet"):
              from trace_core.network import to_graphml
              files["Integrix_subnetwork.graphml"] = to_graphml(*extras["subnet"])
          st.session_state["fullzip"] = R.bundle(files)
      if "fullzip" in st.session_state:
          dl_card("Full bundle", "Integrix_bundle", "Tables, manifest, methods and figures in one archive", key="ix_fullzip",
                  buttons=[("ZIP", st.session_state["fullzip"], "Integrix_bundle.zip", "application/zip", "ix_fullzip_zip")])

# ===================================================================================== Methods
with tabs[11]:
  if on(11):
      st.markdown(R.methods_text(res, kb, extras))
      st.markdown("**IPA-style displays.** Canonical-pathway bars, molecule shapes, measured (red/green) and predicted "
                  "(orange/blue) colors, relationship-colored edges and subcellular layouts follow the visual "
                  "conventions of Ingenuity Pathway Analysis. They are computed from Integrix statistics and the "
                  "Human-GEM knowledge base; they do not use IPA's proprietary knowledge base or its activation "
                  "z-score algorithm.")
      with st.expander("Run log and timings"):
          for l in res.log:
              st.write("• " + l)
          st.write({k: f"{v:.1f} s" for k, v in res.timings.items()})
      with st.expander("Knowledge base"):
          st.json(kb.manifest)


# ===================================================================================== next page
# 'Next: <page> →' at the bottom of every result page, in the side-panel order (the same as the other platforms)
_ORDER = [(_s, _t) for _s, _ts in enumerate(SECTIONS) for _t in _ts]
for _t in SEC_IDX:
    with tabs[_t]:
        if on(_t) and (SEC, _t) in _ORDER and _ORDER.index((SEC, _t)) + 1 < len(_ORDER):
            _ns, _nt = _ORDER[_ORDER.index((SEC, _t)) + 1]
            st.divider()
            st.button(f"Next: {TAB_NAMES[_nt].split(': ', 1)[-1]} →", key=f"next_page_footer_{_t}", on_click=_go,
                      args=(_ns, _nt))

# ===================================================================================== footer
st.markdown(UI.footer(__version__, kb.version), unsafe_allow_html=True)
