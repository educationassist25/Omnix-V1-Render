# BulkRNAAI Pro

A standalone Streamlit app for bulk RNA-seq differential expression analysis,
built around raw RNA-seq counts (or pre-normalized logCPM) and **pyDESeq2**.

## Pipeline (11 tabs)

1. **Data Upload** — upload a genes x samples expression matrix + clinical
   metadata (or use the built-in Rich Clinical Demo). Choose input type:
   *Raw Counts* (integers — enables pyDESeq2) or *Pre-normalized logCPM*.
2. **Data Cleaning** — the standard pyDESeq2/DESeq2 workflow, with no
   imputation step for either input mode: for raw counts, low-expression
   gene filtering (CPM-based, the standard `filterByExpr`-style rule, with
   user-adjustable thresholds), with a downloadable filter report and
   filtered count matrix recording every gene removed; any stray non-numeric
   cell is treated as zero (undetected), matching the convention that an
   unmapped gene is a true zero rather than missing data. For logCPM input,
   the supplied matrix is used as-is (no imputation), with any missing values
   excluded gene-by-gene, automatically, by downstream statistics.
3. **QC Validation** — library size & genes-detected per sample; per-sample
   count/expression distribution boxplot; sample-sample Pearson correlation
   across all biological samples (plus QC-replicate CV and correlation, if QC
   replicates are present); automatic outlier-sample flagging (low mean
   correlation to the cohort and/or PCA-centroid distance); all figures and
   tables downloadable.
4. **Normalization** — pyDESeq2 median-of-ratios size factors + a
   variance-stabilizing transform (VST) for raw counts; pass-through for
   logCPM input. Raw counts, linear-scale normalized counts, and the VST/log2
   transformed matrix are each kept separate and separately downloadable, and
   a PCA-before-vs-after-normalization comparison is shown once normalization
   completes.
5. **PCA** — dynamic group selection (any subset — 2, 3, 4 or more groups),
   preset palettes or per-group custom colors, per-group marker styles,
   95% confidence ellipses, score/loading/variance plots, and downloads for
   scores, explained variance, loadings, and the input matrix.
6. **Statistics** — choose the comparison type:
   - **Two-group**: pyDESeq2 Wald test (raw counts) or a linear-model Wald
     test (logCPM) — same statistic and output schema either way, with
     explicit `Group_A`/`Group_B`/`Higher_In` columns. For raw counts, the
     model is fit with the comparison group as the reference level so
     pyDESeq2's own recommended follow-up step, **apeGLM LFC shrinkage**
     (`DeseqStats.lfc_shrink`), applies automatically: `Log2FC`/`Linear_FC`
     report the shrunk estimate (pyDESeq2's recommendation for ranking/
     volcano plots), with the original unshrunk MLE kept alongside in
     `Log2FC_MLE_Unshrunk`/`Linear_FC_MLE_Unshrunk` — shrinkage never changes
     the p-value/FDR. An optional supplementary **Welch's two-sample t-test
     on the VST/log2 values** can also be run for the same comparison —
     reported in its own table, kept clearly separate from the DESeq2
     model-based result.
   - **3+ groups**: either classical **one-way ANOVA on transformed data**
     (VST/logCPM, with tukey/dunnett/pairwise post-hoc), or, for raw counts,
     **PyDESeq2 count-based contrasts** — a single DESeq2 GLM fit across the
     selected groups with pairwise Wald contrasts drawn from it, each with
     its own Benjamini-Hochberg FDR. "Every group vs a reference" (Dunnett-
     style) fits with that reference as the model baseline, so LFC shrinkage
     applies to every contrast from the one shared fit; "all pairwise" has no
     single reference coefficient to shrink, so those contrasts report the
     unshrunk MLE estimate. Ordinary ANOVA is never applied directly to raw
     counts.
   Covariates supported throughout. Comparison-named CSV downloads for every
   result and post-hoc/contrast table.
7. **Volcano Plot** — full publication-grade control panel: thresholds,
   colorblind-friendly palettes and per-category color overrides, point
   style, label modes (top-N / all significant / manual) with auto-repel,
   legend, axes, title, gridlines, gene highlighting, background, figure
   size, stats box, themes; exports to PNG/PDF/SVG/JPEG/TIFF plus
   up/down/complete data tables and a settings JSON.
8. **Biomarker Discovery** — candidate gene shortlist by p-value/FDR/combined
   criteria, with pathway breakdown when row annotations are loaded.
9. **Heatmap** — dynamic group subsetting, stacked column annotation tracks
   (categorical or continuous), row (pathway) tracks, per-track color
   customization, live FDR/p-value feature filtering, clustering and distance/
   linkage options, custom gradients or preset palettes, discrete breakpoints,
   figure sizing, fonts, multi-format export, and the Z-score table.
10. **Boxplot** — per-gene expression across any subset of groups, per-group
    colors, mean/median/point toggles, panel layout and figure sizing,
    statistics table, and multi-format export.
11. **GSEA** — three distinct methods: preranked **GSEA** (Subramanian et al.
    2005), hypergeometric **ORA** with an explicit background/universe, and
    **STRING** enrichment. Gene sets come from an organism gene-set database
    or your own `.gmt` file. Enrichment (running-score) plots are annotated
    with ES, NES, p-value and FDR, and a reproducibility metadata panel
    records the library and version used.

## Rich Clinical Demo dataset

`generate_rich_demo_data.py` builds `sample_counts_matrix_richdemo.csv`
(2,178 genes x 46 samples), `sample_metadata_richdemo.csv`, and
`gene_row_annotations_richdemo.csv`:

- 4 diagnosis groups — Healthy Control, Prediabetic, Type 2 Diabetes,
  Metabolic Syndrome (10 biological samples each) — plus 6 pooled QC
  replicates.
- Clinical metadata: Diagnosis, Age, Gender, Treatment, Ethnicity,
  Body Weight, Batch.
- Counts are simulated from a negative-binomial model (as DESeq2/edgeR
  themselves assume), with realistic library sizes (~9-28M reads), a
  mean-dependent dispersion trend, batch effects, and genuine
  diagnosis-driven differential expression across 15 curated real-gene-symbol
  pathways (insulin signaling, inflammation, OXPHOS, adipogenesis, hypoxia,
  etc. — the same pathways in the bundled `gene_sets_demo.gmt`), so
  Statistics/Volcano/Heatmap/GSEA all return genuinely meaningful,
  reproducible results out of the box.

Regenerate with:
```
python3 generate_rich_demo_data.py
```

## Running locally

```
pip install -r requirements.txt
streamlit run app.py
```

## Notes on pyDESeq2

- pyDESeq2 requires **integer raw counts** (samples x genes internally —
  this app transposes for you) and fits its own size factors and
  dispersions; feeding it already-normalized data violates its model, so
  the pre-normalized logCPM input path uses classical statistics instead.
- Design: list covariates to control for (e.g. Batch) first, and the primary
  variable of interest last (`~ Batch + Diagnosis`), matching DESeq2/R
  convention.
- Contrasts are Wald tests: `contrast=[factor, level_A, level_B]` gives
  log2FC of level_A relative to level_B (level_A is the numerator).

## Consistency of the two input modes

Both input modes run a **Wald test** — coefficient ÷ standard error, referenced
to a normal distribution — and return the identical result schema
(`baseMean`, `Log2FC`, `lfcSE`, `stat`, `p-value`, `FDR`, `Linear_FC`,
`Significant`), so every downstream tab behaves the same either way.

They are not the same underlying model, because the input does not allow it:

| | Raw Counts | Pre-normalized logCPM |
|---|---|---|
| Model | Negative-binomial GLM (pyDESeq2) | Gaussian linear model (OLS) |
| Dispersion shrinkage | Yes (DESeq2) | No |
| Test | Wald | Wald |
| Covariates | Yes | Yes |

DESeq2's model is defined on counts and cannot be fitted to values that are no
longer counts, so the logCPM path applies the same Wald construction to a
linear model — what limma/voom-style workflows do with already-normalized
expression values. Relative to a Student's t-test on the same logCPM data, the
Log2FC is identical and the p-value differs only by referencing the normal
rather than a t-distribution (slightly anti-conservative in very small samples).

## Offline gene sets

`gene_sets_demo.gmt` (15 pathways, matched to the Rich Clinical Demo's real
gene symbols) still ships with the app. It is no longer a built-in GSEA option,
but you can select it through **Upload custom .gmt file** to run GSEA fully
offline.
