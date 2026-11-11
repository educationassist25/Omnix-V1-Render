"""Screen-identical figure files (Kaleido). Skipped where no Chrome is available to Kaleido."""
import io
import os
import sys
import zipfile

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from trace_core import plotly_export as PX  # noqa: E402

pytestmark = pytest.mark.skipif(not PX.available(), reason="Kaleido cannot start Chrome here")


def _fig():
    import plotly.graph_objects as go
    f = go.Figure(go.Bar(x=["a", "b"], y=[1, 2]))
    f.update_layout(height=300)
    return PX.screen_figure(f)


def test_every_format_at_the_requested_size():
    from PIL import Image
    f = _fig()
    for ext in ("png", "tiff", "jpg"):
        im = Image.open(io.BytesIO(PX.figure_bytes(f, ext, 300, 4.0, 3.0)))
        assert im.size == (1200, 900) and round(im.info["dpi"][0]) == 300
    assert PX.figure_bytes(f, "svg", 300, 4.0, 3.0).lstrip().startswith(b"<svg")
    assert PX.figure_bytes(f, "pdf", 300, 4.0, 3.0).startswith(b"%PDF")
    eps = PX.figure_bytes(f, "eps", 300, 4.0, 3.0)
    assert eps.startswith(b"%!PS") and b"%%BoundingBox: 0 0 288 216" in eps and len(eps) < 2e6


def test_large_raster_is_capped_and_zip_has_all_six():
    assert PX.effective_dpi(20, 20, 1200) < 1200
    z = zipfile.ZipFile(io.BytesIO(PX.all_formats_zip(_fig(), "fig", 300, 3, 2)))
    assert sorted(n.rsplit(".", 1)[1] for n in z.namelist()) == sorted(PX.FORMATS.values())


def test_webgl_dots_are_drawn_in_files():
    """Scattergl draws nothing in the server's Chrome at print resolution: files must use SVG dots."""
    import plotly.graph_objects as go
    from PIL import Image
    f = go.Figure(go.Scattergl(x=[0, 1, 2], y=[0, 1, 2], mode="markers", marker=dict(color="red", size=40)))
    f.update_layout(height=300, xaxis_visible=False, yaxis_visible=False, showlegend=False)
    sf = PX.screen_figure(f)
    assert [t.type for t in sf.data] == ["scatter"]
    im = Image.open(io.BytesIO(PX.figure_bytes(sf, "tiff", 600, 6.0, 3.0))).convert("RGB")
    red = sum(1 for r, g, b in im.resize((300, 150)).getdata() if r > 200 and g < 80 and b < 80)
    assert red > 50


def test_typeset_for_print():
    """Journal layout: column width, real point sizes, no title, height from content, tight to the chart."""
    import plotly.graph_objects as go
    from PIL import Image
    from trace_core import printfig as PRT
    n = 12
    f = go.Figure(go.Bar(x=list(range(1, n + 1)), y=[f"Pathway number {i} with a long name" for i in range(n)],
                         orientation="h"))
    f.update_layout(title="screen title", height=900)
    t, w, h = PRT.typeset(PX.screen_figure(f), "canonical_pathways", 3.5, 8)
    assert w == 3.5 and not t.layout.title.text and 1.5 < h < 6
    assert abs(t.layout.font.size - 8 * 4 / 3) < 0.01 and t.layout.yaxis.tickfont.size >= 6 * 4 / 3
    im = Image.open(io.BytesIO(PX.figure_bytes(t, "png", 300, w, h))).convert("L")
    assert im.size == (1050, int(round(h * 300))) or abs(im.size[1] - h * 300) < 3
    bbox = Image.eval(im, lambda v: 255 - v).getbbox()      # drawn content reaches near every edge
    assert bbox[0] < 30 and bbox[1] < 60 and im.size[0] - bbox[2] < 60 and im.size[1] - bbox[3] < 60
    assert PRT.widths_for("reaction_patterns") and 3.5 not in PRT.widths_for("reaction_patterns").values()


def test_eps_of_a_heatmap_colour_bar():
    import plotly.graph_objects as go
    f = go.Figure(go.Heatmap(z=[[1, 2], [3, 4]], colorbar=dict(title="v")))
    eps = PX.figure_bytes(PX.screen_figure(f), "eps", 300, 4, 3)
    assert eps.startswith(b"%!PS") and len(eps) < 2e6
