## Overview

Transcriptomics takes one bulk RNA-seq dataset (a gene expression matrix plus a clinical metadata table) through cleaning, quality control, normalization, differential expression, visualization, biomarker selection and pathway analysis, all in one place.

It accepts two kinds of input:

- **Raw Counts**: integer read counts per gene per sample. These are analyzed with pyDESeq2 (the Python implementation of DESeq2): median-of-ratios size factors, a variance-stabilizing transform (VST), and negative-binomial Wald tests.
- **Pre-normalized logCPM**: values that are already log2(CPM+1)-style, for example exported from edgeR or limma-voom. These are used as supplied. Two-group tests use a linear-model Wald test, and multi-group tests use one-way ANOVA.

The sidebar lists seven workflow steps. Clicking a step shows its pages as buttons at the top of the main area. Work through the steps in order, because each page uses earlier results. Loading a new dataset clears all downstream results. Most figures have a download row underneath (**Format**: PNG, JPEG, TIFF, SVG or PDF; **DPI**: 150, 300 or 600 for raster formats).

## Before you start: input files

### Expression matrix (required)

- CSV, XLSX or XLS, with genes as rows and samples as columns (column headers = sample names).
- The first column holds gene identifiers, used exactly as written (no ID conversion). Duplicates are rejected. HGNC symbols work better than Ensembl IDs for pathway analysis.
- For **Raw Counts**, negative values cause an error, fractional counts (e.g. Salmon/tximport) are rounded, and non-numeric cells are treated as zero during cleaning.

### Clinical metadata (required)

- CSV, XLSX or XLS format, with one row per sample.
- **Sample** (required): must match the expression-matrix column names exactly. Every sample in the matrix must appear here.
- **Group** (required): the main experimental or clinical group.
- **IsQC** (optional): marks technical QC replicates. The values True, 1, yes and qc (in any capitalization) count as QC. If the column is missing, no sample is treated as QC. QC samples are used only on the QC Validation page and are left out of normalization and statistics.
- **Batch** (optional): if missing, every sample gets batch "1".
- Any other clinical columns (Diagnosis, Gender, Treatment, and so on) can be added. A column can be used for grouping, coloring or as a covariate if it is text, or if it is numeric with 15 or fewer distinct values. Continuous columns such as Age are not offered as groups, but they can still be shown as heatmap annotation tracks.

### Gene row annotations (optional)

A CSV/XLSX/XLS file whose first column holds gene identifiers that match the matrix, followed by columns such as **Pathway**. This file adds:
- row annotation tracks on the Heatmap, and
- extra columns plus a pathway breakdown chart on Biomarker Discovery.

### Demo dataset

Tick **Use built-in Rich Clinical RNA-seq Demo dataset** to load a simulated raw-count study (it always loads in Raw Counts mode): 2,178 genes x 46 samples, with four diagnosis groups (Healthy Control, Prediabetic, Type 2 Diabetes, Metabolic Syndrome; 10 each) plus 6 pooled QC replicates. It includes Diagnosis, Age, Gender, Treatment, Ethnicity, Body Weight and Batch metadata, and a Pathway annotation file.

**Example files in Excel (Help → Data preparation)**
Ready-made workbooks of the demo data, one per file type. Each has a **Data** sheet (the table to upload; Omnix reads the first sheet) and a **Read me** sheet describing every column.
- Matrix: raw counts, 2,178 genes × 46 samples.
- Metadata: Sample, Group, IsQC, Batch and the clinical variables.
- Row annotation: Gene and Pathway.

Replace the values with your own and upload the workbook as it is, or save it as CSV.

## Setup

### Data Upload

1. Choose **Input data type**: **Raw Counts** or **Pre-normalized logCPM**.
2. Either tick the demo checkbox or upload:
   - **Expression matrix (CSV/XLSX)**
   - **Clinical metadata (CSV/XLSX)**
   - optionally, the gene row-annotation file
3. Click **Load Dataset**. A message reports the number of genes and samples, split into biological and QC samples.
4. Check the **Loaded Data Preview**, which shows the first 10 genes and the full metadata table. The caption lists the grouping variables the app detected.
5. Optional: click **Download Raw_Expression_Matrix.csv** to save the matrix as loaded.

## Preprocessing

### Data Cleaning

The app does not impute missing values in either mode.

**Raw Counts**

1. Set **Minimum CPM** (0.1 to 10, default 1.0).
2. Set **Minimum fraction of samples meeting CPM threshold** (0.05 to 1.0, default 0.2).
3. Click **Run Gene Filtering**. CPM is computed per biological sample, and a gene is kept if its CPM is at or above the minimum in at least that fraction of samples (rounded up; at least one sample).
4. Review the retained-gene count, the **Gene Filtering Summary** figure (kept vs. filtered, and a log10 mean-CPM histogram with the threshold line) and the per-gene report (**Mean_CPM**, **N_Samples_Detected**, **Kept**).
5. Download **Filter_Report.csv** and **Filtered_Counts_Matrix.csv**.

**Pre-normalized logCPM**

The page reports how many values are missing, previews the matrix, and offers **Download Cleaned_LogCPM.csv**. Just visiting the page is enough to move on. Missing values are dropped gene by gene later, wherever a test cannot use them.

### QC Validation

This page runs automatically once cleaning is done and shows up to eight figures:

- **Library Size per Sample** and **Genes Detected per Sample** (raw counts only; computed on the full, unfiltered gene set).
- **QC CV Distribution**, **Feature Count by CV Quality** (CV ≤ 20% = "Acceptable", otherwise "Variable") and **QC Sample Correlation Matrix**. These need samples marked IsQC.
- **Count / Expression Distribution Across Samples** (per-sample boxplots), **Sample Correlation Analysis (All Biological Samples)** (Pearson) and **Gene Expression Distribution**.

Tables can be downloaded as **Library_Size_Summary.csv**, **QC_CV_Table.csv**, **QC_Sample_Correlation_Matrix.csv** and **Sample_Correlation_Matrix.csv**.

Under **Potential Outlier Samples**, a sample is flagged if its mean correlation to the other samples is more than 2.5 robust (MAD-based) deviations below the cohort median, or if its distance from the PC1/PC2 centroid is more than 3 SD above average. **QC_Outlier_Report.csv** gives the reason for each flag. **Nothing is removed automatically**, so review flagged samples before normalizing.

## Normalization

**Raw Counts**

1. Choose the **Primary grouping variable (design factor)**. This is usually your main biological variable, such as Diagnosis.
2. Optionally, choose **Additional covariates to control for (optional, e.g. Batch)**. The design becomes `~ covariates + primary variable`.
3. Click **Run DESeq2 Normalization**. pyDESeq2 estimates median-of-ratios size factors, fits dispersions and the GLM, and computes the VST. If the installed pyDESeq2 version lacks VST, log2(normalized counts + 1) is used instead.
4. Review the figures: **DESeq2 Size Factors** (reference line at 1.0), **Before: Raw Counts Distribution** / **After: VST-Normalized Distribution**, and **PCA — Before Normalization** (log2(CPM+1)) / **PCA — After Normalization** (VST). The before/after PCA needs at least 4 samples. If samples cluster by batch or depth before normalization but by biology after, normalization has worked.
5. Download **Normalized_Counts.csv** (linear, size-factor corrected), **Normalized_Expression_Matrix.csv** (VST/log2, the matrix used by PCA, Heatmap and Boxplot) and **PCA_Before_After_Scores.csv**.

**Pre-normalized logCPM:** click **Use logCPM Matrix As-Is**. No transformation is applied, and **Normalized_Counts.csv** is disabled.

## Statistical Analysis

### Statistics

1. Pick the **Grouping variable:**. Diagnosis is the default if present, otherwise Group.
2. Choose the **Comparison type**: **Two-group comparison** or **ANOVA (≥3 groups)**.

#### Two-group comparison

1. Select **Group A** and **Group B**. They must be different; if they are the same, the run button does nothing. **Log2FC is Group A relative to Group B**, so a positive value means the gene is higher in Group A.
2. Under **Significance Threshold**, set **Metric** (**FDR** or **P-value**) and **Threshold** (default 0.05). By default the table and CSV show all genes, with a Significant column; tick **Apply significance filter** to show only genes with metric ≤ threshold.
3. Set **Covariates to control for (optional, e.g. Batch)** if needed.
4. Click **Run Two-Group Test**:
   - **Raw Counts** (**pyDESeq2 Wald test (recommended for raw counts)**): DESeq2 is refit on the filtered counts with this page's design and Group B as the reference. Fold changes are shrunk with apeGLM, so **Log2FC** and **Linear_FC** are shrunken values (unshrunk ones are kept as **Log2FC_MLE_Unshrunk** / **Linear_FC_MLE_Unshrunk**). Shrinkage does not change p-values or FDR.
   - **Pre-normalized logCPM** (**Wald test (linear model on logCPM)**): a per-gene linear model with the same Wald construction and output columns. It needs at least 3 samples; genes with missing values are dropped.
5. The table shows **Gene, p-value, FDR, Log2FC, Linear_FC, Significant**. Download it as **\<A\>_vs_\<B\>_Statistics.csv**.
6. Optional (raw counts): tick the Welch's t-test checkbox and click **Run Welch t-test (transformed data)** for a separate, descriptive test on the VST values (**..._Welch_Transformed_Statistics.csv**). It does not replace the DESeq2 result.

#### Three or more groups

With raw counts, choose a **Multi-group comparison type**. With logCPM input, only ANOVA is available.

**PyDESeq2 count-based contrasts (model-based, recommended for raw counts)**

1. Select the groups to include (at least 3).
2. Choose **Contrasts**:
   - **Every group vs a reference (Dunnett-style)**: pick a **Reference group**. LFC shrinkage is applied to every contrast.
   - **Every pairwise combination**: fold changes are unshrunk MLE estimates.
3. Set covariates and the **Significance Threshold** controls.
4. Click **Run PyDESeq2 Contrasts**. One DESeq2 fit is made, and each contrast gets its own Benjamini-Hochberg FDR.
5. Use **View contrast:** to choose which contrast to display. The contrast you are viewing becomes the active result for Volcano, Biomarker, Heatmap and GSEA.
6. Download **\<contrast\>_DESeq2_Contrast.csv** or **All_DESeq2_Contrasts.csv** (all contrasts, with a Contrast column).

**ANOVA on transformed data (VST/logCPM)**

1. Select at least 3 groups.
2. Choose a **Post-hoc test**: **tukey** (Tukey HSD), **dunnett** (Welch t-tests of each group against the first group in sample order, BH-corrected) or **pairwise** (Welch t-tests for every pair, BH-corrected).
3. Set the threshold controls and click **Run ANOVA**. The table shows the F-statistic, ANOVA p-value, FDR and Significant.
4. Post-hoc results are calculated only for significant genes. Browse them with **View post-hoc results for gene:**.
5. Downloads: **ANOVA_\<groups\>.csv**, the post-hoc table for a single gene, and **..._Posthoc_All_Genes.csv**.

## Exploratory Analysis

### PCA

1. Choose **Color/group samples by:** and the groups to include. At least 3 samples are needed.
2. Set **Number of components** (default 5, maximum 10).
3. Under **Customize Appearance**, choose a **Color mode** (**Preset palette** or **Custom colors (pick each group)**), **Show 95% confidence ellipses** (on by default) and a marker shape for each group.
4. You get a score plot, a loading plot and a variance-explained plot. **Top Contributing Genes (Loadings)** lists the 20 genes with the largest loadings.
5. Download the scores (with group labels), explained variance, loadings, and the exact input matrix. File names include the selected groups.

### Volcano Plot

The volcano plot uses the active Statistics result. The settings are grouped into collapsible panels:

- **1 — Statistical Thresholds**: **Y-axis metric** (**p-value** or **FDR**), **Significance cutoff** (0.05), **Fold change threshold (log2 units)** (default 0), optional **Also require FDR <**, and **Show linear FC on x-axis instead of log2FC**. A gene is "Up" if its metric is below the cutoff and its log2FC is above the threshold, and "Down" if its metric is below the cutoff and its log2FC is below the negative threshold.
- **4 — Gene Labels**: **Top N significant** (default 10), **All significant**, **Manually selected** or **None**, with automatic label repelling.
- The remaining panels control colors (including colorblind-safe palettes), point style, legend, axes, title, gridlines, gene highlighting, background, figure size and a publication theme (a stylistic approximation only).

Export the figure as PNG, PDF, SVG, JPEG or TIFF (300, 600 or 1200 DPI for raster formats). Below the plot are the **Top 20 Biomarkers** table and downloads for the upregulated, downregulated and complete volcano data, plus a **Volcano_Figure_Settings.json** for reproducibility.

### Heatmap

1. Choose **Filter samples by:** and the groups to include.
2. Pick the **Column annotation tracks:** from any metadata columns. Categorical columns become color blocks; numeric columns become gradients.
3. If you loaded an annotation file, also pick **Row annotation tracks:**. You can set per-category colors in **Customize Annotation Track Colors** (for tracks with 12 or fewer distinct values).
4. **Feature Selection**: choose **Filter by** (FDR or p-value) and set the slider (≤ 0.05 by default). The number of genes that pass is shown live. Above 150 genes, row labels are hidden.
5. **Clustering Options**: **Clustering** (Rows only by default; rows and/or columns, or none), **Distance metric** (euclidean or correlation) and **Linkage method** (ward, average or complete).
6. **Color Scale Customization**: a predefined palette or a custom low/mid/high gradient, **Z-score minimum/maximum** (-2.5 / 2.5) and optional discrete breakpoints. Panel size and fonts can also be adjusted.
7. The heatmap shows row-scaled z-scores of the normalized matrix. Export it as PNG, PDF, SVG, JPEG or TIFF; cells are drawn as an image, so use 300 DPI or more even for PDF/SVG. **Zscore_Table.csv** matches the figure's order.

### Boxplot

1. Search for and select one or more genes, then choose the **Grouping variable:** and at least 2 groups.
2. Adjust the appearance if needed: **Panels per row** (default 3), data points, mean and median lines, and group colors.
3. Statistics are calculated on the normalized matrix: a Welch's t-test for 2 groups, or a one-way ANOVA for 3 or more. The FDR is corrected only across the genes you selected, so use the Statistics page for FDR values across the whole gene panel.
4. Download the statistics CSV and the figure in PNG, PDF, SVG, JPEG or TIFF.

## Biomarker Discovery

1. Set the **p-value cutoff** (default 0.05) and the **FDR cutoff** (default 0.05).
2. Choose the **Selection criterion**: both FDR and P value below their cutoffs (default), P value only, or FDR only.
3. Click **Run Biomarker Discovery**. The table lists p-value, FDR, Linear_FC, Log2FC and Direction (Up/Down) for each candidate gene, plus any columns from your gene annotation file. If the annotation file has a Pathway column, a bar chart shows how many candidates fall in each pathway.
4. Download **\<A\>_vs_\<B\>_Biomarkers.csv**.

## Pathway Analysis

### GSEA

This page analyzes the active Statistics result.

1. Set the significant-gene cutoff with **Filter by** (FDR or p-value) and **\< threshold** (default 0.05). The caption shows how many genes pass. STRING and ORA use this gene list. GSEA ignores it and ranks every gene.
2. Choose the **Analysis type**:

**STRING Enrichment Analysis** (needs internet access)
- Choose the **Organism**: human, mouse, rat, yeast, zebrafish, or one of five plant species.
- Click **Run STRING Enrichment Analysis**. Your significant genes are sent to string-db.org, and enrichment is computed there against STRING's own GO, KEGG, Reactome, Pfam and InterPro annotations.

**Over-Representation Analysis (ORA)** (needs internet access to download the library from Enrichr)
- Choose the **Organism** and **Gene-set database**. Human has GO (BP/MF/CC), KEGG, Reactome, WikiPathways and MSigDB (Hallmark). Mouse has KEGG and WikiPathways; yeast and zebrafish have GO and WikiPathways; rat and plants have none (use STRING). Unavailable databases are marked "unavailable", with the reason.
- Choose the **Background / universe**: **All detected genes in this study (recommended)** or **Entire gene-set library**.
- Set **Min gene set size** (1) and **Max gene set size (0 = no limit)**.
- Click **Run Over-Representation Analysis**. A one-sided hypergeometric test is run locally. The results include P-value, FDR, overlap counts, Fold Enrichment and the overlapping genes.

**Gene Set Enrichment Analysis (GSEA)**
- Under **Gene sets from:**, pick one of:
  - **Organism gene-set database**, then click **Fetch gene set library** (needs internet access); or
  - **Upload custom .gmt file**, which runs fully offline. The app ships with a 15-pathway **gene_sets_demo.gmt** that matches the demo dataset.
- Choose the **Ranking metric**: **signed_neglogp** (default; sign(Log2FC) × −log10(p-value)), **log2fc**, or **statistic** (see Tips).
- Set **Min gene set size** (15), **Max gene set size** (500) and **Permutations (per gene set)** (100 to 2000, default 500).
- Click **Run Gene Set Enrichment Analysis (GSEA)**. The results give ES, NES, P-value and BH FDR for each gene set, plus the leading-edge genes. Significance comes from gene-set permutation.

**Results (all methods)**
- A results table (**Download GSEA_Results.csv**, the same file name for all three methods) and an **Analysis metadata (for reproducibility)** panel (**GSEA_Analysis_Metadata.csv**).
- A **Pathway Dot Plot**: x-axis is NES (GSEA), Fold Enrichment (ORA) or Gene Ratio (STRING); dot size is gene count; color is significance. Options include **Show top** and **Selection** (**Top up + down** for GSEA only). The figure and its data can be downloaded.
- GSEA only: an **Enrichment Plot** for a chosen gene set, drawn with **Generate Enrichment Plot** and annotated with NES, p-value and FDR.

## Tips & troubleshooting

- **Upload errors**: both files are required unless you use the demo. "Metadata missing required column(s)" means you need **Sample** and **Group** columns. "Samples present in expression matrix but missing from metadata" means the names must match exactly. "Duplicate gene identifiers found" means duplicates must be collapsed or renamed first.
- **"Negative values found in a matrix marked 'Raw Counts'"**: the data is probably already log-transformed. Reload it as **Pre-normalized logCPM**.
- **A page says to run an earlier step first**: follow the order Upload → Cleaning → Normalization → Statistics. Volcano, Heatmap, Biomarker Discovery and GSEA all need a Statistics result.
- **No QC CV or QC correlation figures**: the metadata has no samples with IsQC = True.
- **Nothing happens when you click Run Two-Group Test**: Group A and Group B are the same group.
- **"Need at least 3 values of '...'"** for multi-group tests: choose a grouping variable that has 3 or more groups.
- **"Not enough samples ... for the requested design"** (logCPM): remove a covariate or add samples.
- **After an ANOVA**, the active result has placeholder fold changes (Log2FC = 0). Volcano Up/Down calls, biomarker Direction, and GSEA rankings are therefore not meaningful. For directional analysis, use a two-group test or a DESeq2 contrast.
- **Heatmap: "Only N significant gene(s) at this cutoff"**: clustering needs at least 2 genes, so relax the threshold.
- **GSEA "statistic" ranking**: the result table stores the Wald statistic in a column named "stat", which this option does not recognize. As a result, "statistic" currently gives the same ranking as **signed_neglogp**.
- **"Could not reach ..." errors**: STRING, ORA and database-backed GSEA need internet access. To work offline, use GSEA with an uploaded .gmt file.
- **ORA returns no gene sets / GSEA finds few matches**: use gene symbols that match the organism, such as HGNC symbols for human. Gene matching ignores capitalization, but Ensembl IDs will not match symbol-based libraries.

## Glossary

- **CPM (counts per million)**: a gene's count divided by the sample's total counts, times one million.
- **Size factor**: DESeq2's median-of-ratios scaling factor for each sample, correcting for sequencing depth and RNA composition. Near 1.0 means little rescaling.
- **VST (variance-stabilizing transform)**: DESeq2's log2-like transform that makes variance roughly independent of expression level. Recommended for PCA, clustering and heatmaps.
- **Design / covariate**: the model formula, e.g. `~ Batch + Diagnosis`. Covariates are variables whose effects are controlled for.
- **Wald test**: an estimated coefficient divided by its standard error, compared with a normal distribution.
- **log2 FC (Log2FC)**: log2 of Group A expression relative to Group B (+1 = twice as high in A). **Linear_FC** = 2^Log2FC.
- **apeGLM shrinkage**: pulls noisy fold changes for low-count genes toward zero without changing p-values.
- **padj / FDR**: the Benjamini-Hochberg adjusted p-value (DESeq2 calls it padj), which controls the false discovery rate across tested genes.
- **CV (coefficient of variation)**: SD divided by the mean, as a percentage, across QC replicates.
- **Z-score (heatmap)**: each gene's values scaled to mean 0 and SD 1 across the displayed samples.
- **ORA**: a hypergeometric test of whether a gene set contains more significant genes than expected by chance.
- **ES / NES**: GSEA's enrichment score and its size-normalized version. A positive NES means the set sits toward the top of the ranking (higher in Group A with signed_neglogp or log2fc).
- **Leading edge**: the genes that drive a gene set's enrichment score.
