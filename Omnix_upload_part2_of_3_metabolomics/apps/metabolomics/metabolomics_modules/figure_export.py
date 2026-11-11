"""
Figure downloads: one panel and one set of options for every Omnix platform.

Under every figure: Format (PNG, TIFF, JPEG, SVG, PDF, EPS or all six in one ZIP) · Resolution (300, 600 or
1200 DPI; 600 by default) · Download. The file is built only when the button is pressed.

* Raster formats (PNG, TIFF, JPEG) are written at the chosen DPI. Very large figures are capped at
  60 million pixels so a 1200 DPI export cannot exhaust the server's memory; the panel says when that
  happens. PDF, SVG and EPS are vector files (sharp at any size); cells drawn as images (heatmaps) are
  embedded at the chosen DPI.
* TIFF is a flat RGB image with lossless LZW compression and the DPI tag set; JPEG is written on white at quality 95. Text in SVG, PDF and EPS
  stays editable.

The same file is kept in each platform's module folder, so the platforms behave identically.
"""
from __future__ import annotations

import functools
import io
import threading
import zipfile

FORMATS = {               # label -> (extension, MIME type)
    "PNG": ("png", "image/png"),
    "TIFF": ("tiff", "image/tiff"),
    "JPEG": ("jpg", "image/jpeg"),
    "SVG": ("svg", "image/svg+xml"),
    "PDF": ("pdf", "application/pdf"),
    "EPS": ("eps", "application/postscript"),
}
ALL_FORMATS = "All formats (ZIP)"
VECTOR = ("svg", "pdf", "eps")
DPI_OPTIONS = (300, 600, 1200)
DEFAULT_DPI = 600
MAX_MEGAPIXELS = 60.0
_LOCK = threading.Lock()        # downloads are built in a worker thread; one savefig per figure at a time


def effective_dpi(fig, dpi: int) -> int:
    """The requested DPI, lowered if the figure would exceed MAX_MEGAPIXELS at that resolution."""
    w, h = fig.get_size_inches()
    mp = w * dpi * h * dpi / 1e6
    if mp <= MAX_MEGAPIXELS:
        return int(dpi)
    return max(72, int(dpi * (MAX_MEGAPIXELS / mp) ** 0.5))


def figure_bytes(fig, fmt: str = "png", dpi: int = DEFAULT_DPI) -> bytes:
    """The figure as a file: fmt is png, tiff, jpg/jpeg, svg, pdf or eps."""
    ext = {"jpeg": "jpg", "tif": "tiff"}.get(fmt.lower(), fmt.lower())
    kw = {"format": "jpeg" if ext == "jpg" else ext, "dpi": effective_dpi(fig, dpi),
          "bbox_inches": "tight", "pad_inches": 0.06}
    if ext == "jpg":
        kw.update(facecolor="white", pil_kwargs={"quality": 95})
    elif ext == "tiff":
        kw.update(format="png", facecolor="white")     # converted below to a flat RGB TIFF
    elif ext == "eps":
        kw["facecolor"] = "white"          # EPS has no transparency
    buf = io.BytesIO()
    import matplotlib
    # vector files keep real, editable text (Illustrator, Inkscape): SVG <text>, TrueType fonts in PDF / EPS
    with _LOCK, matplotlib.rc_context({"svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42}):
        fig.savefig(buf, **kw)
    if ext == "tiff":
        # journals, Word, Prism and Illustrator expect a flat RGB TIFF (no alpha) with the DPI tag set
        from PIL import Image
        im = Image.open(io.BytesIO(buf.getvalue()))
        if im.mode in ("RGBA", "LA", "P"):
            im = im.convert("RGBA")
            bg = Image.new("RGB", im.size, "white")
            bg.paste(im, mask=im.split()[-1])
            im = bg
        out = io.BytesIO()
        im.convert("RGB").save(out, format="TIFF", compression="tiff_lzw", dpi=(kw["dpi"], kw["dpi"]))
        return out.getvalue()
    return buf.getvalue()


def all_formats_zip(fig, base_filename: str, dpi: int = DEFAULT_DPI) -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for ext, _ in FORMATS.values():
            z.writestr(f"{base_filename}.{ext}", figure_bytes(fig, ext, dpi))
    return buf.getvalue()


def render_download(st_module, fig, base_filename: str, key_prefix: str):
    """The download row under a figure: Format · Resolution · Download (built on click)."""
    c1, c2, c3 = st_module.columns([1.5, 1.1, 2.2], vertical_alignment="bottom")
    fmt = c1.selectbox("Format", list(FORMATS) + [ALL_FORMATS], key=f"{key_prefix}_fmt",
                       help="PNG, TIFF and JPEG are high-resolution images; SVG, PDF and EPS are vector files "
                            "(sharp at any size, editable text). The ZIP holds all six.")
    dpi = c2.selectbox("Resolution (DPI)", list(DPI_OPTIONS), index=DPI_OPTIONS.index(DEFAULT_DPI),
                       key=f"{key_prefix}_dpi",
                       help="Pixels per inch for PNG, TIFF and JPEG (and for heatmap cells inside vector files).")
    if fmt == ALL_FORMATS:
        name, mime = f"{base_filename}_all_formats.zip", "application/zip"
        data = functools.partial(all_formats_zip, fig, base_filename, dpi)
    else:
        ext, mime = FORMATS[fmt]
        name = f"{base_filename}.{ext}"
        data = functools.partial(figure_bytes, fig, ext, dpi)
    c3.download_button(f"Download {'ZIP' if fmt == ALL_FORMATS else fmt}", data, file_name=name, mime=mime,
                       icon=":material/download:", key=f"{key_prefix}_dl", use_container_width=True)
    w, h = fig.get_size_inches()
    eff = effective_dpi(fig, dpi)
    note = (f"{name} · {w:.1f} × {h:.1f} in" + ("" if fmt in ("SVG", "PDF", "EPS") else
                                       f" · {int(w * eff):,} × {int(h * eff):,} px at {eff} DPI"))
    if eff < dpi and fmt not in ("SVG", "PDF", "EPS"):
        note += f" (lowered from {dpi} DPI for this large figure; use PDF or SVG for full quality)"
    st_module.caption(note)
