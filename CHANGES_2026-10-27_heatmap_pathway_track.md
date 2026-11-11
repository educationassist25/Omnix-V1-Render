# Omnix update — 27 Oct 2026: MSEA pathway track on the Metabolomics/Lipidomics heatmap

Changed files (Omnix_upload_part2_of_3_metabolomics/apps/metabolomics/). Nothing else changed.
- metabolomics_modules/pathway_analysis.py — the MSEA page keeps its latest result for the heatmap
- app.py — Heatmap: new "MSEA pathway" row track, its settings, and a Pathway column in the Z-score table

How it works
1. Run Statistics, then open Pathway Analysis → MSEA (choose KEGG, SMPDB or LIPID MAPS).
2. On the Heatmap, "Pathway track from MSEA" appears under Row (Metabolite) Annotation:
   rank pathways by FDR / P-value / Holm, number of pathways shown (default 8), min metabolites per pathway (default 3).
3. Tick "MSEA pathway" in Row annotation tracks; optionally "Group metabolites into module bands by: MSEA pathway".
Each metabolite gets its best-ranked enriched pathway; others are labeled "Other". Your own row-annotation
columns stay available as before. If Statistics is rerun, the track waits until MSEA is rerun for that result.

Tests: metabolomics 88 passed, Omnix portal 61 passed.

## Update: band names
- metabolomics_modules/heatmap_module.py, app.py — when rows are grouped into bands, the group names are no longer
  written on the left; each band is identified by its color track and the legend. A checkbox
  "Write group names beside the bands" (off by default) brings the names back. When grouping, the grouping column
  is shown as a color track automatically so the legend always names the bands.
