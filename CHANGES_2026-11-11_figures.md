# Omnix: Integrix chart error, overlapping size legend, black figure text (2026-11-11)

## StreamlitDuplicateElementId in Integrix (Results)
Choosing 15 pathways on the Pathways tab made its canonical-pathway chart identical to the one in Overview. On
Streamlit 1.52 both tabs are drawn in the same run, and two identical charts got the same automatic ID.
Every Integrix chart and table now has its own key. New test: `test_same_chart_twice_on_one_page`.
Files: `apps/integrix/app.py` (part 5), `tests/test_integrix_portal.py` (part 1).

## Size legend circles overlapping (Integrix pathway activity dot plot)
Each legend row is now as tall as the largest circle, so the q 0.05 / 0.001 / 1e-6 circles never touch. The class
labels next to the dots are black (were teal).
File: `apps/integrix/trace_core/pubfig.py` (part 5).

## All figure text in pure black, every platform
New `omnix_portal/figure_style.py`, switched on in `portal.py`:
- Matplotlib figures (all platforms, on screen and in every download): titles, axis labels, tick labels, legends,
  colour bars and annotations are drawn in #000000.
- Interactive (Plotly) charts in Integrix and Fluxomics: the same.
- Text that has to stay light (white labels on dark heatmap cells or coloured boxes) and Integrix's dark theme are
  left unchanged.
Files: `omnix_portal/figure_style.py` (new), `omnix_portal/portal.py` (part 1).

## Checks
- Portal tests: 77 passed, 1 skipped on Streamlit 1.52.2 and on 1.67.
- Integrix tests: 49 passed.
- Browser: the Proteomics QC page shows all figure text in black.

## Update: stacked colour bars removed
- Charts without a colour scale (bar charts, the reaction pattern map, …) showed several empty colour bars on
  top of each other at the right, with overlapping numbers. The black-text step created them; it now restyles
  only colour bars a chart really has (`omnix_portal/figure_style.py`).
- Integrix overview "pathway activity z" chart: the first pathway name is no longer cut at the left edge on
  narrow screens (`trace_core/figures.py`).
- Integrix reaction pattern map: the four quadrant captions ("coherent: enzyme up …", "enzyme down, balance up" …)
  now sit in empty bands above and below the data, so dots and enzyme labels never cover them
  (`trace_core/figures.py`).

## Comparisons read "case over control"
Every two-group comparison is now written as **case over control** (for example "Resistant over Parental"),
in the order of the fold change, in all six platforms: result notes, figure titles, volcano plot titles,
table and file names (`Resistant_over_Parental_Statistics.csv`), contrast lists, the Integrix contrast banner,
report and demo descriptions, Help and tutorials.
- Metabolomics/Lipidomics and Fluxomics: Log2FC = Group B − Group A, written "B over A".
- Proteomics and Transcriptomics: Log2FC = Group A − Group B, written "A over B".
- Integrix: "Case over Reference".
Sample counts such as "(3 vs 4 samples)" are unchanged.

## Same group pickers in every platform
Every two-group comparison now has the same two pickers, in the same order: **Case** (the group of interest)
first, then **Control (reference)**. Results always read *case over control*, and a positive fold change always
means higher in the case. Before, "Group A / Group B" meant opposite things in different platforms.
- Statistics in Metabolomics/Lipidomics, Proteomics and Transcriptomics; Fluxomics Step 7; Metabolomics
  Correlation Rewiring; Integrix Study design; the Transcriptomics multi-group "Control (reference)".
- Defaults pick a non-control group as the case (Transcriptomics demo: Type 2 Diabetes over Healthy Control).
- Fluxomics volcano legend: "Higher in case" / "Lower in case".
- Tutorials and Help use the same words; fold-change definitions read "case − control".
- The statistics themselves are unchanged.

## Integrix: every figure downloads exactly as shown
- Under every chart: **Download this figure as shown**. The file is the on-screen chart itself (same Plotly
  figure, light style, black text), as PNG, TIFF, JPEG, SVG, PDF, EPS or all six (ZIP), at 300 / 600 / 1200 DPI;
  the size starts at the chart's size on screen and can be changed. Before, downloads were a separately drawn
  "publication" version that looked different.
- Layer concordance and Reaction classification can now be downloaded too.
- Export tab: the figure set (ZIP) and the full bundle use the same screen figures, with each page's current
  settings. The targets list, cascade and pathway maps are drawn for print and keep their publication layout.
- How: Plotly's image engine Kaleido (new `trace_core/plotly_export.py`; `kaleido` and `svglib` added to
  `requirements.txt`). Kaleido needs a headless Chrome: `render_start.sh` downloads it in the background at
  start-up. If it cannot run on the server, each panel says so and offers the publication-style figure instead.
- Raster files are capped at 60 million pixels; SVG, PDF and EPS are vector files.
- Fix: charts with many dots (reaction pattern map, gene-protein concordance) are drawn on screen with WebGL,
  which the server's Chrome cannot draw at print resolution, so their files showed axes and legend but no dots.
  Files now draw those dots as ordinary SVG dots (same look; vector in SVG, PDF and EPS).

## Integrix: downloads typeset as journal figures
Review of the first version: the files copied the browser layout (a large canvas with tall rows and wide
margins), so the text became too small once the figure was reduced to a journal column, with empty space
around the chart. New `trace_core/printfig.py` keeps the chart (data, marks, colours) and sets it for print:
- **Width:** single column 3.5 in (89 mm), 1.5 column 5.5 in (140 mm) or double column 7.2 in (183 mm). Dense
  charts offer only the widths their labels fit: reaction pattern map and concordance 1.5 or double column,
  layer concordance double column.
- **Type:** 8 pt by default (7–12 pt to choose); tick labels and legend 1 pt smaller, never below 6 pt;
  Arial/Liberation Sans, black; gene and enzyme labels at most 7.5 pt; long axis titles wrap instead of being cut.
- **Height from the content** (rows, wrapped labels, legend lines); tight margins; no title inside the image.
- **Chart details:**
  - The canonical pathways chart has no empty band above the first and below the last bar.
  - The pathway activity dot plot has a key across the figure whose size circles are drawn at exactly the dot
    sizes of the grid.
  - Pathway names in the layer chart are written in full, wrapped under each group.
  - Heatmap values are black on light cells and white on dark ones.
  - The y-axis titles "Features" and "Permutations" are added.
- **Reference lines:** the template drew shapes at 30 % opacity, which hid the p = 0.05 threshold and the
  observed-AUC line. They are now drawn solid, on screen too.
- **EPS:** gradient colour bars are converted to 64 solid bands, so heatmap EPS files work and stay vector.
- **Checked:** every chart in all six formats, at the default and at the narrowest width with 10 pt text. No
  blank files, no clipped or overlapping text.

## Integrix: no overlapping text with long names
Every chart was stress-tested with very long names (110-character pathways, 27-character enzyme and gene
names, 50-character group names, long panel features) and with the demo names: 288 renderings (every width,
7 / 8 / 10 / 12 pt, and the screen charts at half and full width). An automatic check measured every text
element (rotation included). Result: no overlaps and no text outside the figure.
- **Exact text widths:** measured from the font outlines (Liberation Sans / Arial) instead of a character-count
  estimate.
- **Long names:** wrapped within their column, breaking long words at - _ / . when needed. When a name needs
  more than the allowed lines (3 for pathways), the last line ends with "…". Before, the overflow was merged
  into one overlong line.
- **Scatter maps:** gene and enzyme labels are placed again for each print size, with the quadrant captions
  treated as occupied. A label that cannot be placed cleanly is left out rather than drawn over another; the
  screen charts follow the same rule.
- **Bar charts with wrapped names:** every row gets the height of the tallest name.
- **Legends and axis titles:** legend entries and long axis titles wrap to the figure width.
- **Dot plot headers:** the column headers stand vertical when 45 degrees would let them touch.
- **Narrow widths:** the factor variance heatmap is offered at 1.5 and double column only.
- **Screen charts:** two-line chart titles are no longer clipped at the top. The panel chart wraps long feature
  names and its axis title.
