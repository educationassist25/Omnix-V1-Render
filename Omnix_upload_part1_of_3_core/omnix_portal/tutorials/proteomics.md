## Overview

**Proteomics** analyzes one LC-MS/MS protein quantification dataset at a time: one protein intensity matrix and one sample metadata table. You take it through cleaning, QC, normalization and statistics, then on to visualization, biomarker lists, pathway enrichment and protein-protein interaction networks. It supports two quantification methods:

- **Label-Free Proteomics**: each sample is acquired separately and normalized across samples.
- **TMT Proteomics**: samples are multiplexed, usually with a pooled/bridge reference channel, and every other channel can be normalized as a ratio to that reference.

The sidebar lists the workflow steps in order: **Setup**, **Preprocessing**, **Normalization**, **Statistical Analysis**, **Exploratory Analysis**, **Biomarker Discovery**, **Pathway Analysis** and **Network Analysis**. Click a step and its pages appear as a row of buttons at the top of the main area. Once you start working, a collapsible **Processing Log** appears under the navigation. It records every step you applied, with its settings, which is useful for writing methods sections.

## Before you start: input files

All tables can be CSV or XLSX. CSV files saved in Windows/Excel encodings are handled automatically.

### Protein intensity matrix (required)

- **First column:** protein identifier. It is renamed to "Protein" internally, whatever its header, and must be unique. Duplicate protein names cause a validation error.
- **Optional second column:** gene symbol. It is recognized if the header is `Gene`, `Gene Name`, `Gene_Name`, `GeneSymbol`, `Gene Symbol`, `Gene Names` or `Genes` (case-insensitive). The app keeps it as a protein-to-gene mapping rather than as data. It is strongly recommended for GSEA and P-P Interaction, because STRING and Enrichr gene sets use gene symbols.
- **Remaining columns:** one column per sample with raw intensities, including QC replicates and, for TMT, the reference channel. Non-numeric entries become missing values. Blank cells and zeros are both treated as "not detected" by default.

### Sample metadata (required when you upload your own data)

| Column | Required | Meaning |
|---|---|---|
| `Sample` | Yes | Must match the sample column names in the matrix. Every matrix column must appear here. |
| `Group` | Yes | Experimental group or condition. |
| `IsQC` | No | `True`/`1`/`yes`/`QC` marks QC replicates. If the column is absent, all samples are treated as biological. |
| `Batch` | No | Batch label. If absent, every sample gets batch "1". Used by batch-based IQR normalization. |

Any additional columns are kept, for example Diagnosis, Gender or Treatment. Text columns, and numeric columns with 15 or fewer distinct values, become selectable grouping and coloring variables in Statistics, PCA, Heatmap and Boxplot. Continuous columns such as Age can still be shown as heatmap annotation tracks.

### Protein row annotations (optional)

This table needs a `Protein` column plus any annotation columns, such as Method or Pathway. It is used only for the Heatmap's row annotation tracks. If the file has no `Protein` column, it is ignored with a warning.

### Demo datasets

Tick **Use built-in demo dataset instead** (on by default when no file is uploaded), then pick a **Demo dataset**:

- **Standard**: the app picks the demo that matches your quantification method.
  - Label-Free: 243 proteins with no reference channel.
  - TMT: 60 proteins plus a pooled reference channel (`Reference_Pool`).
  - Both have 24 biological samples in four groups (Control/Mild/Moderate/Severe) and 6 QC replicates.
- **Rich Clinical Demo**: 200 proteins and 43 samples.
  - 36 biological samples in four diagnosis groups: Healthy Control, Prediabetic, Type 2 Diabetes and Metabolic Syndrome.
  - 6 QC replicates and 1 pooled reference channel.
  - Metadata: Diagnosis, Age, Gender, Treatment, Ethnicity and Body Weight.
  - Also loads Method and Pathway row annotations.

All demo matrices include a Gene column and contain zeros, so run Data Cleaning & Imputation before normalizing.

**Example files in Excel (Help → Data preparation)**
Ready-made workbooks of the demo data, one per file type. Each has a **Data** sheet (the table to upload; Omnix reads the first sheet) and a **Read me** sheet describing every column.
- Matrix: Label-Free (Protein, Gene, 30 samples) and TMT (with the `Reference_Pool` channel).
- Metadata: Label-Free and TMT (the TMT file lists `Reference_Pool` with Group = Reference).
- Row annotation: Protein, Gene, Method and Pathway.

Replace the values with your own and upload the workbook as it is, or save it as CSV.

## Setup

### Data Upload

1. Choose the **Quantification method**: **Label-Free Proteomics** or **TMT Proteomics**. This choice sets the default normalization later.
2. Under **Protein Intensity Matrix**, use **Upload protein intensity matrix (CSV/XLSX)**, or leave the demo option ticked.
3. Under **Sample Metadata**, use **Upload metadata table (CSV/XLSX)**.
4. Optionally, under **Protein Row Annotations (optional)**, use **Upload row annotation table (CSV/XLSX)**.
5. Click **Load & Validate Data**. A success message reports protein and sample counts (QC vs. study), any row annotation columns, and whether a Gene column was detected.
6. Check the **Preview** of the first 10 proteins, the full **Metadata** table and, if loaded, the **Row Annotations**.

Loading new data clears all earlier statistics, ANOVA and biomarker results.

## Preprocessing

### Data Cleaning & Imputation

This page is optional, but you need it if your data has blanks or zeros. The next steps require strictly positive values. It works on biological samples only; QC replicates are never included.

1. Decide on **Treat exact-zero values as missing (common LC-MS/MS convention for 'not detected')**. It is on by default.
2. **Step 1: Assess & Filter by Missingness**. Counters show how many proteins fall in each band (**Keep (<20%)**, **Careful imputation (20-50%)**, **Usually remove (50-70%)**, **Remove unless essential (>70%)**), next to a **Missingness Distribution** histogram. **View missingness table** gives per-protein counts.
3. Set **Remove features with missingness above this threshold (%)**. The default is 50, in steps of 5. The caption shows how many proteins would be retained.
4. Click **Apply Missingness Filter**.
5. **Step 2: Choose an Imputation Method**. If no missing values remain, the app tells you imputation isn't needed. Otherwise pick an **Imputation method**; the default is K-Nearest Neighbors.

   | Method | What it does | Extra settings |
   |---|---|---|
   | **Half-Minimum (LOD/2)** | Replaces a missing value with half of that protein's minimum observed value. | None |
   | **Mean** / **Median** | Replaces a missing value with that protein's observed mean or median. | None |
   | **K-Nearest Neighbors (KNN)** | Estimates the value from the most similar samples. | **Number of neighbors (k)**: 2–15, default 5 |
   | **Random Forest (MissForest-style)** | Iteratively predicts the value from all other proteins. | **Number of trees**: 5–50, default 10. **Iterations**: 1–10, default 3 |
   | **BPCA (approximated)** | Bayesian-PCA-style imputation, approximated with Bayesian ridge regression. | **Iterations**: 1–15, default 5 |
   | **QRILC (approximated)** | Left-censored imputation that draws values from the lower tail of each protein's distribution. | None |

6. Click **Apply Imputation**. Random Forest and BPCA can take a moment.
7. Review **Result: Complete Post-Imputation Data** and download **Imputed_Data.csv**.

**Downloads:** Missingness Distribution figure and Imputed_Data.csv.

### QC Validation

This page needs at least 2 QC samples. With fewer, the module shows a warning and is skipped.

1. Open the page. The coefficient of variation (CV) is computed for each protein across the QC replicates, using raw intensities.
2. Read the counters **Total Features**, **Acceptable (CV≤20%)** and **Variable (CV>20%)**.
3. Review the four panels. Each has its own figure download.
   - **QC CV Distribution**, with a line at 20%.
   - **Feature Counts by CV Quality**.
   - **Protein Intensity Distribution**, on a log10 scale.
   - **QC Sample Correlation**: Pearson correlation of the log2-transformed QC replicates only.
4. Inspect the **CV Filtering Table** (Mean, SD, CV(%), Quality) and download it as **single_CV_Table.csv**.
5. Optionally tick **Filter out 'Variable' features (CV>20%) from downstream analysis**. It is off by default, so all proteins are kept unless you opt in.
6. Click **Confirm QC & Proceed**. The cleaned and imputed data (or, if you skipped imputation, the raw biological samples) moves forward, and QC samples are excluded from everything that follows.

## Normalization

### Normalization

If you skipped **Confirm QC & Proceed**, this page uses the raw biological samples. Imputed data is only carried forward through the QC confirmation.

1. **Step 1: Log2 Transformation**. Click **Apply Log2 Transformation**.
   - This is a strict log2 with no pseudo-count. It fails if any zero, negative or missing value remains.
   - Download **Log2_Transformation.csv**.
2. **Step 2: Normalization**. Choose a **Normalization method**.
   - **Reference-Channel Normalization** (TMT only; the TMT default).
     - Select the pooled/bridge channel in **Select reference/pooled channel (sample column)**, then click **Apply Reference-Channel Normalization**.
     - Each channel becomes log2(channel) − log2(reference).
     - The reference column is dropped from the output.
   - **IQR Normalization** (Label-Free only).
     - Robust scaling: (X − median) ÷ IQR, applied to the log2 values.
     - Choose a **Normalization axis**:
       - **sample** (default): each sample is scaled across proteins.
       - **feature**: each protein is scaled across samples.
       - **batch**: uses the metadata `Batch` column.
     - Click **Apply IQR Normalization**.
   - **Median Centering Normalization** (both modes; the Label-Free default).
     - Subtracts each sample's median and adds back the grand median.
     - Click **Apply Median Centering Normalization**.
   - **Global MAD-based Variance Scaling** (both modes).
     - Divides all values by one dataset-wide MAD × 1.4826.
     - Click **Apply Global MAD-based Variance Scaling**. The success message shows the scale factor.
   - In TMT mode, the Median Centering and MAD options also ask you to **Select reference/pooled channel (sample column) to exclude from this normalization**. That channel is left out of the calculation and dropped from the output.
3. Review the **Distribution Diagnostics** (raw vs. log2, then log2 vs. normalized, as distributions and box plots).
4. Review the **Current Working Matrix — Log2 & Normalized**. This is the matrix every later page uses.

Changing the method clears any earlier normalized result, so results are never mixed between methods.

**Downloads:** Log2_Transformation.csv, Normalization_Table.csv, Normalized_Log2.csv and both diagnostic figures.

## Statistical Analysis

### Statistics

1. Choose the **Grouping variable:** (defaults to `Group`).
2. Choose the **Comparison type**: **Two-group comparison** or **ANOVA (≥3 groups)**.

**Two-group comparison**

1. Pick **Group A**, **Group B** and a **Method**: **Welch's t-test** (unequal variances) or **Wilcoxon rank-sum**.
   - In the current version, the results table is always computed with an unequal-variance (Welch) t-test, whichever Method you pick.
   - The Method you selected is still written to the Processing Log.
2. Click **Run Two-Group Test**. The table contains:
   - **Log2FC**: mean of Group A minus mean of Group B on the normalized log2 data. Positive means higher in Group A.
   - **p-value**.
   - **FDR**: Benjamini-Hochberg.
   - **Significant**: true when p < 0.05 and FDR < 0.25.
3. Download **<GroupA>_vs_<GroupB>_Statistics.csv**.

This two-group result drives the Volcano Plot, Heatmap, Biomarker Discovery, GSEA and P-P Interaction pages.

**ANOVA**

1. This mode needs at least 3 values of the grouping variable. Select 3 or more in the multiselect.
2. Choose a **Post-hoc test**:
   - **tukey**: Tukey HSD.
   - **dunnett**: t-tests of each group against the first group, with BH correction.
   - **pairwise**: all pairs by t-test, with BH correction.
3. Click **Run ANOVA**. You get F-statistic, ANOVA p-value and FDR for each protein.
4. Post-hoc results are computed only for proteins with FDR < 0.25. Use **View post-hoc results for feature:** to inspect one protein.

**Downloads:** ANOVA table, post-hoc results for one protein, and post-hoc results for all proteins.

## Exploratory Analysis

### PCA

This page uses the normalized matrix, biological samples only.

1. Choose **Color/group samples by:** and select which groups to include. You need at least 3 samples in total.
2. Set the **Number of components**: 2 up to 10, default 5.
3. Under **Customize Appearance**, choose a **Color mode**: **Preset palette** or **Custom colors (pick each group)**.
4. Toggle **Show 95% confidence ellipses** (on by default) and pick a marker for each group.
5. You get a score plot, a loading plot (top 20 proteins) and a variance-explained plot, plus the **Top Contributing Proteins (Loadings)** table.

**Downloads:** each figure; PCA_Scores.csv, PCA_Explained_Variance.csv, PCA_Loadings.csv and PCA_Input_Data.csv, each prefixed with the selected groups.

### Volcano Plot

Run a two-group test in Statistics first.

1. **1 — Statistical Thresholds**:
   - **Y-axis metric**: **p-value** (cutoff default 0.05) or **FDR** (cutoff default 0.25).
   - **Fold change threshold (log2 units)**: default 0.
   - Optional **Also require FDR <**.
   - **Show linear FC on x-axis instead of log2FC**.
   - **Threshold line style**.
   - A protein counts as Up or Down when it passes the cutoff and its Log2FC is beyond ± the fold-change threshold.
2. Optional styling sections cover colors and point style (including colorblind-safe palettes), protein labels (**Top N significant**, **All significant**, **Manually selected** or **None**), legend, axes and title, gridlines, highlighting specific proteins, background and figure size, and a publication theme (for example Nature, Cell or PNAS).
3. Export the figure as PNG, PDF, SVG, EPS, JPEG or TIFF. DPI options are 300, 600 and 1200, for raster formats only.
4. Review the **Top 20 Biomarkers** table.

**Downloads:** Upregulated_Proteins.csv, Downregulated_Proteins.csv, Complete_Volcano_Data.csv, and Volcano_Figure_Settings.json (for reproducibility).

### Heatmap

This page requires the two-group Statistics result.

1. Choose **Filter samples by:** and select the groups to show.
2. **Column (Sample) Annotation**: pick any metadata columns as annotation tracks.
3. **Row (Protein) Annotation**: shown only if you loaded row annotations. The caption reports how many proteins matched.
4. Optionally open **Customize Annotation Track Colors**.
5. **Feature Selection**: **Filter by** FDR (default ≤ 0.25) or p-value (default ≤ 0.05). **Proteins at this cutoff** updates live; above 150 proteins, row labels are hidden.
6. **Clustering Options**: **Clustering** (**Rows only** by default, **Both rows and columns**, **Columns only** or **No clustering**), **Distance metric** (euclidean or correlation) and **Linkage method** (ward, average or complete).
7. **Color Scale Customization**: a preset palette or custom low/mid/high gradient, and **Z-score minimum**/**Z-score maximum** (default −2.5 and 2.5). **Figure Size** and **Title & Font Size** adjust the layout.

The heatmap shows row-scaled z-scores.

**Downloads:**
- The heatmap as PNG, PDF, SVG, JPEG or TIFF, at 150, 300 or 600 DPI. DPI matters even for PDF and SVG here.
- Zscore_Table.csv, in the clustered order.

### Boxplot

1. Use **Search and select one or more proteins**.
2. Choose the **Grouping variable:** and at least 2 groups.
3. Customize font, **Panels per row**, figure size, group colors, and **Show individual data points** / **Show mean (dashed line)** / **Show median (solid line)**.
4. The **Statistical Results** table uses Welch's t-test for 2 groups and one-way ANOVA for 3 or more. Its FDR is corrected only across the proteins shown. For panel-wide FDR, use Statistics.

**Downloads:** Boxplot_Statistics.csv, and the figure as PNG, PDF, SVG, JPEG or TIFF.

## Biomarker Discovery

### Biomarker Discovery

1. Choose the **Filtering criterion**: **p<0.05 AND FDR<0.25**, **p<0.05 only** or **FDR<0.25 only**.
2. Adjust the **p-value cutoff** (default 0.05) and **FDR cutoff** (default 0.25). The criterion you choose decides which cutoffs apply.
3. Click **Discover Biomarkers**. The table lists p-value, FDR, Log2FC and **Direction** (Up/Down), sorted by p-value.
4. Download **<GroupA>_vs_<GroupB>_Biomarkers.csv**.

## Pathway Analysis

### GSEA

This page uses the two-group Statistics result. If a Gene column was provided, protein IDs are translated to gene symbols.

1. Set **Filter by** (FDR or p-value) and a threshold. The defaults are FDR < 0.25 and p < 0.05. This builds the significant list used by STRING and ORA; GSEA ranks all proteins instead.
2. Choose an **Analysis type**.

**STRING Enrichment Analysis** (needs internet)

1. Pick the **Organism**. The page shows the STRING database version.
2. Click **Run STRING Enrichment Analysis**. Enrichment is computed on STRING's servers, against STRING's own GO, KEGG, Reactome, Pfam and InterPro annotations.

**Over-Representation Analysis (ORA)** (needs internet to fetch the gene-set library)

1. Pick the **Organism** and a **Gene-set database**:
   - GO Biological Process, GO Molecular Function or GO Cellular Component;
   - KEGG, Reactome, WikiPathways or MSigDB (Hallmark).
   - Combinations that don't exist for an organism are marked "unavailable" and explained. Human has all of them; rat and the plant species have none, so use STRING for those.
2. Choose the **Background / universe**: **All detected proteins in this study (recommended)** or **Entire gene-set library**.
3. Set **Min gene set size** (default 1) and **Max gene set size (0 = no limit)**.
4. Click **Run Over-Representation Analysis**. This runs a one-sided hypergeometric test locally, with BH-FDR correction.

**Gene Set Enrichment Analysis (GSEA)** (preranked)

1. Choose where gene sets come from under **Gene sets from:**:
   - **Organism gene-set database**, then click **Fetch gene set library** (needs internet); or
   - **Upload custom .gmt file**, which runs fully offline.
2. Choose the **Ranking metric**:
   - **signed_neglogp**: sign(Log2FC) × −log10(p). This is the default.
   - **log2fc**.
   - **statistic**: the Statistics table has no test-statistic column, so this currently falls back to signed_neglogp.
3. Set **Min gene set size** (default 15), **Max gene set size** (default 500) and **Permutations (per gene set)** (100–2000, default 500).
4. Click **Run Gene Set Enrichment Analysis (GSEA)**.
5. Under **Enrichment Plot**, pick a gene set and click **Generate Enrichment Plot** to see the running-score plot. This plot is only available for GSEA runs.

**Results for all three methods**

- A results table, downloadable as **GSEA_Results.csv**.
- **Analysis metadata (for reproducibility)**: database, version, gene counts and so on, downloadable as **GSEA_Analysis_Metadata.csv**.
- **Pathway Dot Plot**: choose the **Number of top pathways to show** (3–30, default 15), then click **Generate Pathway Dot Plot**. The x-axis is NES (GSEA), Fold Enrichment (ORA) or Gene Ratio (STRING); color shows −log10(FDR).

## Network Analysis

### P-P Interaction

This page needs internet access, because it calls the public STRING API.

1. Set the significance filter, the same way as on the GSEA page. At least 2 proteins must pass.
2. Choose the settings:
   - **Species**: Human, Mouse, Rat, Yeast or Zebrafish.
   - **Confidence threshold**: Low (0.15), Medium (0.4, the default), High (0.7) or Highest (0.9).
   - **Layout**: spring, circular or kamada_kawai.
3. Click **Fetch STRING Network**.
4. Read the **Network Enrichment Statistics**: Nodes, Edges, Avg. Degree and **PPI Enrichment p-value**.
5. Adjust **Network Visualization**: **Node size by** (pvalue or degree), **Show protein labels**, and the Log2FC node color scale (preset palette or custom gradient; **Log2FC color scale (±)** defaults to 2).
6. Review the **Interaction Table** and **Hub Proteins (Degree, Betweenness, Clustering)**.

**Downloads:** the network figure, PPI_Interactions.csv and PPI_Hub_Table.csv.

## Tips & troubleshooting

- **"Load data in Setup → Data Upload first."**: go to **Setup → Data Upload** and click **Load & Validate Data**.
- **Validation errors on load:**
  - "Metadata missing required column(s)": add `Sample` and `Group`.
  - "Samples present in peak matrix but missing from metadata": every matrix column needs a metadata row.
  - "Duplicate protein names found": make the first column unique.
- **"Please upload a metadata table…"**: uploading your own matrix also requires a metadata file.
- **Log2 fails with "requires every value to be strictly positive"**: zeros or missing values remain. Apply the missingness filter and imputation, then click **Confirm QC & Proceed** so the imputed data reaches Normalization.
- **Fewer than 2 QC samples:** the QC module is skipped, and Normalization then uses the raw, unimputed biological samples. Your data must already be free of zeros and blanks.
- **"Could not compute the QC correlation matrix"**: the QC replicates contain zeros or blanks, which log2 cannot handle. The CV results are still shown.
- **Statistics, PCA and Boxplot say "Complete the Normalization page (Log2 transform) first."**: you must finish Step 2 (Normalization), not just Step 1.
- **Volcano, Heatmap, Biomarker, GSEA and P-P Interaction ask for a two-group comparison:** running ANOVA alone is not enough.
- **ANOVA warnings:** you need at least 3 values of the grouping variable, and at least 3 selected.
- **Heatmap "clustering needs at least 2" features:** relax the FDR/p-value cutoff.
- **"No proteins pass the current cutoff"** (GSEA) or **"Need at least 2 proteins"** (P-P Interaction): relax the threshold.
- **"Could not reach STRING (string-db.org)…"** or **"Could not reach … Enrichr"**: check your internet connection. For offline pathway analysis, use GSEA with an uploaded .gmt file.
- **ORA returns no gene sets:** none of your genes were found in the universe. This is often because no Gene column was provided and the Protein IDs are not gene symbols. Check the metadata panel for the mapping counts.
- **GSEA "Only N ranked genes available…"**: lower **Min gene set size** or use a larger dataset.

## Glossary

- **Log2 FC:** difference of group means on the log2 scale (Group A − Group B). A value of 1 means two-fold higher in Group A.
- **FDR:** false discovery rate; here the Benjamini-Hochberg adjusted p-value.
- **CV:** coefficient of variation (SD ÷ mean × 100). CV ≤ 20% across QC replicates is flagged Acceptable.
- **QC replicate:** repeated injection of a pooled sample, used to judge reproducibility. Excluded from normalization and statistics.
- **TMT reference channel:** pooled/bridge sample in each TMT plex. Other channels are expressed as log2 ratios to it.
- **IQR / MAD:** interquartile range and median absolute deviation, robust measures of spread used for scaling.
- **Z-score (heatmap):** a protein's value minus its mean, divided by its SD across the samples shown.
- **ORA:** over-representation analysis; a hypergeometric test of whether a gene set holds more significant genes than expected.
- **GSEA:** ranks all proteins and tests whether a gene set's members cluster at the top or bottom of the ranking.
- **NES:** normalized enrichment score. Positive means enriched among proteins ranked higher in Group A; negative means ranked lower.
- **Leading edge:** the gene-set members that drive a GSEA enrichment score.
- **Fold Enrichment (ORA):** observed ÷ expected overlap; 1.0 is what chance alone would give.
- **PPI enrichment p-value:** STRING's test of whether your list has more interactions than a random list of the same size.
- **Degree / Betweenness:** number of connections, and how often a node lies on shortest paths between others; high values flag hub proteins.
