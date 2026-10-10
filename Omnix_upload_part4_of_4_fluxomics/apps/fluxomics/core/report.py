"""
core.report — reproducible analysis report (self-contained HTML) and settings files.
"""
from __future__ import annotations

import base64
import hashlib
import html
import json
import platform
from datetime import datetime

import pandas as pd

from . import __version__
from .interpret import markdown_to_html


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def versions() -> dict:
    import numpy, scipy, matplotlib, plotly, openpyxl
    v = {"Omnix Fluxomics": __version__, "Python": platform.python_version(), "pandas": pd.__version__,
         "numpy": numpy.__version__, "scipy": scipy.__version__, "matplotlib": matplotlib.__version__,
         "plotly": plotly.__version__, "openpyxl": openpyxl.__version__}
    try:
        import streamlit
        v["streamlit"] = streamlit.__version__
    except Exception:
        pass
    return v


def settings_json(settings: dict) -> str:
    return json.dumps(settings, indent=2, default=str)


def _table(df: pd.DataFrame, n: int = 30) -> str:
    if df is None or df.empty:
        return "<p class='muted'>Not computed.</p>"
    d = df.head(n).copy()
    for c in d.columns:
        if pd.api.types.is_float_dtype(d[c]):
            d[c] = d[c].map(lambda x: "" if pd.isna(x) else (f"{x:.3g}" if abs(x) < 1e4 else f"{x:.3e}"))
    more = f"<p class='muted'>Showing {min(n, len(df))} of {len(df)} rows (full table in the Excel export).</p>" if len(df) > n else ""
    return d.to_html(index=False, escape=True, border=0, classes="tbl") + more


def build_report(*, title: str, inputs: list, design: pd.DataFrame, settings: dict, log: list,
                 tables: dict, figures: list, interpretation_md: str = "") -> str:
    """
    inputs: [(file name, sha256, rows)]; tables: {section title: DataFrame}; figures: [(caption, png bytes)].
    """
    css = """
    body{font-family:Inter,Segoe UI,Helvetica,Arial,sans-serif;color:#1d1d22;max-width:1100px;margin:0 auto;padding:32px 20px;line-height:1.5;background:#fff}
    h1{font-size:26px;margin:0 0 4px}h2{font-size:19px;margin-top:34px;border-bottom:1px solid #e4e3e8;padding-bottom:6px}
    h3{font-size:15px;margin-top:20px}.muted{color:#000000;font-size:13px}
    table.tbl{border-collapse:collapse;font-size:12px;margin:8px 0;width:100%}
    table.tbl th{background:#f3f2f7;text-align:left;padding:5px 7px;border-bottom:1px solid #dcdbe3}
    table.tbl td{padding:4px 7px;border-bottom:1px solid #eeedf2;font-variant-numeric:tabular-nums}
    figure{margin:18px 0}figure img{max-width:100%;border:1px solid #eeedf2}figcaption{font-size:12px;color:#000000}
    blockquote{border-left:3px solid #6f5fc7;margin:10px 0;padding:4px 12px;background:#f7f6fc}
    pre{background:#f7f6fc;padding:10px;font-size:11.5px;overflow-x:auto}
    """
    H = [f"<!doctype html><html><head><meta charset='utf-8'><title>{html.escape(title)}</title><style>{css}</style></head><body>"]
    H.append(f"<h1>{html.escape(title)}</h1><p class='muted'>Generated {datetime.now():%Y-%m-%d %H:%M} · "
             f"Omnix Fluxomics {__version__}</p>")
    H.append("<h2>1. Input files</h2>")
    H.append(_table(pd.DataFrame(inputs, columns=["File", "SHA-256", "Rows"]), 50))
    H.append("<h2>2. Processing log</h2><ul>" + "".join(f"<li>{html.escape(x)}</li>" for x in log) + "</ul>")
    H.append("<h2>3. Experimental design</h2>")
    H.append(_table(design, 200))
    H.append("<h2>4. Settings</h2><pre>" + html.escape(settings_json(settings)) + "</pre>")
    H.append("<h2>5. Results</h2>")
    for name, df in tables.items():
        H.append(f"<h3>{html.escape(name)}</h3>" + _table(df))
    if figures:
        H.append("<h2>6. Figures</h2>")
        for cap, png in figures:
            b = base64.b64encode(png).decode()
            H.append(f"<figure><img src='data:image/png;base64,{b}'><figcaption>{html.escape(cap)}</figcaption></figure>")
    H.append("<h2>7. Software</h2>" + _table(pd.DataFrame(list(versions().items()), columns=["Component", "Version"]), 50))
    H.append("<h2>8. Methods (template)</h2><p>" + html.escape(methods_text(settings)) + "</p>")
    H.append("</body></html>")
    return "\n".join(H)


def _correction_text(c: dict) -> str:
    from .correction import METHODS, METHOD_REFS, INSTRUMENTS
    mode = c.get("mode", "already")
    if mode != "apply":
        return {"already": "Natural-abundance correction: isotopologue areas were used as supplied (already corrected "
                           "by the vendor software). ",
                "none": "No natural-abundance correction was applied. "}.get(mode, f"Natural-abundance correction: {mode}. ")
    m = c.get("method", "isocor_hr")
    extra = ""
    if m == "accucor":
        extra = (f" ({INSTRUMENTS.get(c.get('instrument', 'orbitrap'), '')}, resolving power "
                 f"{float(c.get('resolving_power', 1e5)):,.0f} at m/z {float(c.get('ref_mz', 200)):g})")
    elif m == "lsq":
        extra = " (" + ("tracer element only" if c.get("resolution", "high") == "high" else "all elements") + ")"
    return (f"Natural-abundance correction was applied in this app: {METHODS.get(m, m)}{extra}, including tracer "
            f"isotopic purity [{METHOD_REFS.get(m, '')}]. ")


def methods_text(s: dict) -> str:
    nrm = s.get("normalization") or {}
    el = (s.get("tracer") or {}).get("element", "13C")
    src = ("Compound Discoverer stable-isotope-labeling exports were reorganized so that each metabolite is represented "
           "by one chromatographic peak from one acquisition method (ranking: not a likely in-source fragment, "
           "detection in labeled samples, number of labeled isotopologues, signal, peak-quality warnings). "
           if s.get("input") == "cd" else "A table of isotopologue peak areas with sample metadata was imported. ")
    return (
        src + "Every isotopologue (M+0…M+n, n = number of " + el + " atoms) was analysed as a separate channel; "
        "isotopologue areas were not summed. "
        + _correction_text(s.get("correction") or {}) +
        "Isotopologue fractions were calculated as the area of each isotopologue divided by the summed area of all "
        "isotopologues M+0…M+n of that metabolite in the same sample (equivalent to Compound Discoverer's “Exchange "
        "Rate in Percent”); incorporation = (1 − fraction M+0) × 100; mean label number = Σ i·fraction(M+i); "
        "fractional enrichment = Σ i·fraction(M+i)/n × 100. "
        + "For intensity analyses, isotopologue rows with a zero value in any labeled sample were removed (a zero "
        "M+0 removed the whole compound); values below the minimum value (default 1,000) in any labeled sample were treated the same way; "
        "metabolites left with only the M+0 isotopologue were excluded; missing values were not imputed. "
        + ("Areas were log2-transformed and normalized per sample by a robust interquartile-range method "
           "(x − median)/IQR × median IQR + median of medians (labeled samples only; 12C reference and blank "
           "samples excluded), computed over "
           f"{'M+0 channels' if nrm.get('reference') == 'm0' else 'all channels'}; "
           f"intensity analyses used the {'log2' if s.get('intensity_dataset') == 'log2' else 'log2 + IQR-normalized'} "
           "data. " if nrm else "")
        + "Two groups were compared with t-tests (Welch, Student or paired, as stated), more than two with ANOVA "
        "(one-way, Welch's, Kruskal–Wallis or two-way group × time, as stated) and post hoc tests; p values and "
        "Benjamini–Hochberg FDR were reported for every isotopologue channel. Label arrival was shown as group-mean "
        "isotopologue fractions (± SEM or 95 % CI) for every group and time point; labeling describes label "
        "appearance and is not absolute flux.")
