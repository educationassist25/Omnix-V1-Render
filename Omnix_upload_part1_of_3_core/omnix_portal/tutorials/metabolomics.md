## Overview

The Metabolomics app takes LC-MS peak areas and runs them through a fixed pipeline: cleaning and imputation, QC review, normalization, statistics, figures, pathway analysis and correlation-network analysis. It supports four analysis types: **Untargeted Metabolomics**, **Targeted Metabolomics**, **Untargeted Lipidomics** and **Targeted Lipidomics**. You can analyze one dataset, or several combined (for example, two methods in positive and negative mode).

The sidebar lists the workflow steps in order. Clicking a step shows its pages as buttons across the top of the main area. Downstream pages need upstream results, and each page warns you when something is missing.

Most tables and figures have a download button. The standard figure download lets you pick a **Format** (PNG, JPEG, TIFF, SVG or PDF) and, for the raster formats (PNG, JPEG, TIFF), a **DPI** of 150, 300 (the default) or 600.

## Before you start: input files

All files can be CSV or XLSX. CSV encoding is detected automatically, so exports that contain characters like ± or µ load fine.

**Peak area matrix (required)**
- First column: the metabolite or lipid name. Names must be unique; duplicates stop the upload.
- Remaining columns: one per injection (biological samples and QC injections), holding peak areas. Non-numeric cells are read as missing.
- For targeted assays, include the dataset's internal standard (ISTD) as an ordinary row, for example `ISTD_D4-Alanine` (metabolomics) or `ISTD_PC(15:0/18:1-d7)` (lipidomics). A name starting with `ISTD` is selected automatically; any other name is selected when it appears in the `ISTD` column of the row annotations.

**Sample metadata (required)**
- **Sample**: must match the peak-matrix column headers exactly. Every sample column in the matrix must be listed here.
- **Group**: the experimental group.
- **IsQC**: `True`, `1`, `yes` or `qc` marks a QC injection. If the column is missing, all samples are treated as biological.
- **Batch** (optional): used by batch-wise Median-IQR normalization. If absent, every sample is set to batch 1.
- Extra columns (Diagnosis, Gender, Treatment and so on) are kept. Text columns, and numeric columns with at most 15 distinct values, become selectable grouping variables in Statistics, PCA, Heatmap and Boxplot.

**Metabolite row annotations (optional)**
- A **Metabolite** column followed by any annotation columns, such as `HMDB_ID`, `Class`, `Method`, `Pathway` or `ISTD`.
- The Heatmap uses them as row tracks. Pathway Analysis uses the HMDB IDs, or matches by name when there are none.
- **ISTD** (targeted assays): the internal standard used to normalize each metabolite or lipid, the same ISTD for the whole dataset (for example `ISTD_D4-Alanine` on every row). The ISTD row names itself. The column may also be called `Internal Standard`.
- A file without a `Metabolite` column is ignored, with a warning.

**Demo datasets**
- **Standard**: 24 biological samples in 4 groups (Control/Mild/Moderate/Severe) plus 6 QCs, one dataset per analysis type:
  - Targeted Metabolomics: 59 metabolites and one ISTD, `ISTD_D4-Alanine`.
  - Targeted Lipidomics: 39 lipids and one lipid ISTD, `ISTD_PC(15:0/18:1-d7)`.
  - Untargeted Metabolomics: 300 metabolites. Untargeted Lipidomics: 90 lipids. No ISTD.
- **Rich Clinical Demo**: 200 metabolites and one ISTD (`ISTD_D4-Alanine`), 36 samples plus 6 QCs, six clinical metadata variables, and Method, Pathway, HMDB ID and ISTD annotations. Targeted types normalize with the ISTD; untargeted types leave the ISTD row out.
- In every targeted demo, the `ISTD` column of the row annotations links each metabolite to the dataset's ISTD. The multi-dataset demo uses `ISTD_D4-Alanine` for Targeted Metabolomics and `ISTD_PC(15:0/18:1-d7)` for Targeted Lipidomics.
- **Demo Data** (multi-dataset mode): one dataset per analysis type, all sharing the same samples.

**Example files in Excel (Help → Data preparation)**
Ready-made workbooks of the demo data, one per file type. Each has a **Data** sheet (the table to upload; Omnix reads the first sheet) and a **Read me** sheet describing every column.
- Matrix: Untargeted metabolomics, Targeted metabolomics (with the ISTD row), Untargeted lipidomics, Targeted lipidomics (with the ISTD row).
- Metadata: one sample metadata file that fits all four matrices.
- Row annotation: one file for each of the four analysis types. The targeted files include the `ISTD` column that links every metabolite or lipid to the internal standard.

Replace the values with your own and upload the workbook as it is, or save it as CSV.

## Setup

### Data Upload

1. Under **How many datasets are you uploading?**, choose **Single dataset** (multi-dataset mode is described later).
2. Choose the **Analysis type**. It sets the normalization pipeline and how adjusted p-values are labelled.
3. Upload the **Peak Area Matrix**, the **Sample Metadata**, and optionally the **Metabolite Row Annotations (optional)**. To use a demo instead, tick **Use built-in demo dataset instead** (it is ticked automatically when no file is uploaded) and pick **Standard** or **Rich Clinical Demo**.
4. Click **Load & Validate Data**. The app reports how many features, QC samples and study samples it found, and shows previews.

Loading new data clears all earlier statistics and biomarker results.

## Preprocessing

### Data Cleaning & Imputation

This step is optional and works on biological samples only.

1. **Treat exact-zero values as missing** is on by default.
2. **Step 1: Assess & Filter by Missingness**
   - Counters sort features into Keep (<20%), Careful imputation (20–50%), Usually remove (50–70%) and Remove unless essential (>70%).
   - Set **Remove features with missingness above this threshold (%)** (default 30), then click **Apply Missingness Filter**.
3. **Step 2: Choose an Imputation Method** appears only if missing values remain. The **Imputation method** options are:
   - **Half-Minimum (LOD/2)**
   - **Mean**
   - **Median**
   - **K-Nearest Neighbors (KNN)** (the default; set k from 2 to 15, default 5)
   - **Random Forest (MissForest-style)**
   - **BPCA (approximated)**
   - **QRILC (approximated)**, which draws from the lower tail of each metabolite's values and assumes missing values fell below the detection limit

   The BPCA and QRILC options are Python approximations, not the original algorithms.
4. Click **Apply Imputation**, then save the result with **Download Imputed_Data.csv**.

### QC Validation

You need at least 2 QC injections.

1. The app computes each feature's CV across the QC injections (raw peak areas). Features at 20% or below are **Acceptable**; the rest are **Variable**.
2. Four figures appear, each downloadable: **QC CV Distribution**, **Feature Counts by CV Quality**, **Peak Area Distribution** and **QC Sample Correlation** (QC replicates only).
3. The **CV Filtering Table** (Mean, SD, CV(%), Quality) downloads as `single_CV_Table.csv`.
4. Optionally tick **Filter out 'Variable' features (CV>20%) from downstream analysis** (off by default).
5. Click **Confirm QC & Proceed**. QC is applied to the imputed data if you imputed, otherwise to the raw samples. From here on, QC injections are excluded from everything else.

## Normalization

### Normalization

The pipeline depends on the analysis type. If you did not confirm QC, the page continues with the biological samples only.

**Targeted Metabolomics / Targeted Lipidomics: ISTD → Log2**

1. **Step 1: ISTD Normalization** is required.
   - **Internal standard (ISTD)** is pre-selected with the dataset's ISTD: the one named in the row annotations' `ISTD` column, or the row whose name starts with `ISTD`. You can pick another row.
   - Before applying, check the number of features normalized, the ISTD's QC CV (after QC is confirmed) and any samples where the ISTD is zero or missing. Warnings flag a CV above 20% or missing ISTD values.
   - **Metabolite → ISTD** lists every feature with its ISTD. Download it as **ISTD_Assignments.csv**.
   - Click **Apply ISTD Normalization**. Every peak area is divided by the ISTD peak area in the same sample, and the ISTD row is removed. If the ISTD is zero or missing in a sample, that sample's values become missing.
   - The ISTD row is never removed by the missingness or QC CV filters, so it is always available here.
2. **Step 2**: click **Apply Strict Log2 Transformation**. This is log2(x) with no offset; zero or negative values become missing, with a warning.
3. Download **Normalized_NoLog2.csv** (the ISTD-normalized data) and **Normalized_Log2.csv** (the final matrix).

**Untargeted Metabolomics / Untargeted Lipidomics: Log2 → Median-IQR**

1. **Step 1**: click **Apply Strict Log2 Transformation**.
2. **Step 2: Median-IQR Normalization (Robust Scaling)**. Choose the **Normalization axis**, then click **Apply Median-IQR Normalization**:
   - **sample** (default): (value − sample median) / sample IQR. The sample statistics use all QC-passed features, even ones removed by the CV filter.
   - **feature**: each metabolite is scaled across samples.
   - **batch**: each metabolite is scaled within each batch, using the **Batch** column.
3. **Distribution Diagnostics** shows before/after plots. Download **Log2_Only.csv** and **Normalized_Log2.csv**.

Internal-standard rows are left out of the untargeted pipeline, since they are not biological features. **Advanced: add an optional ISTD normalization step before Log2** divides every feature by an ISTD first; Log2 → Median-IQR then runs on the ISTD-normalized data.

## Statistical Analysis

### Statistics

1. Choose the **Grouping variable:** (default Group) and a **Comparison type**.

**Two-group comparison**

2. Pick **Group A** and **Group B**, and a **Method**: **Welch's t-test** or **Wilcoxon rank-sum**.
3. Click **Run Two-Group Test**. The results table has four columns:
   - **p-value**
   - the Benjamini–Hochberg adjusted p-value, labelled **BH P Value** for untargeted data and **FDR** for targeted data
   - **Log2FC**: mean log2 of B minus mean log2 of A, so a positive value means higher in B
   - **Linear_FC**

   The success message counts metabolites with p < 0.05 and an adjusted p < 0.25. Download `<A>_vs_<B>_Statistics.csv`.

The Volcano Plot, Heatmap and Biomarker Discovery pages use the latest two-group result.

**ANOVA (≥3 groups)**

2. Select 3 or more groups and a **Post-hoc test**:
   - **tukey**
   - **dunnett**: Welch tests against the first group in your metadata
   - **pairwise**: Welch tests between all pairs
3. Click **Run ANOVA**. You get the F-statistic, ANOVA p-value and adjusted p-value per metabolite. Post-hoc tests run only for metabolites with an adjusted p < 0.25.
4. Download the ANOVA table, a single feature's post-hoc result, or **Posthoc_All_Features.csv**.

## Exploratory Analysis

### PCA

1. Choose **Color/group samples by:** and the groups to include (at least 3 samples in total).
2. Set **Number of components**, **Color palette**, **Show 95% confidence ellipses** (on by default) and a marker for each group.
3. You get a **PCA Score Plot**, a **PCA Loading Plot** (top 20 metabolites) and a **PCA Variance Plot**. Missing values are filled with the metabolite mean before PCA.
4. **Download PCA Data** saves scores, explained variance, loadings and the input matrix as CSV.

### Volcano Plot

This page needs a two-group result.

1. **Quick Style** sets the **Journal theme** (a stylistic approximation), **Color palette**, **Legend position** and **Figure size**.
2. **Statistical Thresholds**:
   - **Y-axis metric**: p-value or adjusted p-value.
   - **Significance cutoff**: default 0.05, or 0.25 for adjusted p.
   - **Fold change threshold (log2 units)**: default 0.
   - An optional second adjusted-p cutoff.
   - **Show linear FC on x-axis instead of log2FC**.

   A metabolite is Up or Down when it is below the cutoff and beyond ± the fold-change threshold.
3. Further expanders cover points, **Metabolite Labels** (Top N significant, All significant, Manually selected, None), legend, axes, gridlines, **Highlight Specific Metabolites**, background and a statistics box.
4. Export the figure as PNG, PDF, SVG, EPS, JPEG or TIFF (300, 600 or 1200 DPI). Under **Export Data**, download the upregulated, downregulated and complete tables as CSV.

### Heatmap

This page needs a two-group result, which decides which metabolites are shown.

1. Choose **Filter samples by:** and the groups. Add **Column annotation tracks:**. If you uploaded row annotations, you can also add row tracks and **Group metabolites into module bands by:** a category. To change track colors, open **Annotation Colors** and click **Apply Colors** (or **Reset Colors**).
2. **Significant feature cutoff**: adjusted p ≤ 0.25 (default), adjusted p ≤ 1, P-value ≤ 0.05, or P-value ≤ 1. Row labels are hidden above 150 metabolites.
3. **Clustering**: Rows only (default), Both rows and columns, Columns only, or No clustering. Also choose a **Distance metric** (euclidean or correlation) and a **Linkage method** (ward, average or complete).
4. Set the palette or a custom gradient, and **Z-score minimum** / **Z-score maximum** (default −2.5 to 2.5). Values are z-scored per metabolite.
5. Export as PNG, PDF, SVG, JPEG or TIFF. Use 300 DPI or higher, even for PDF and SVG, because the cells are drawn as an image. The **Z-score Data Table** also downloads as CSV.

### Boxplot

1. Select metabolites, a **Grouping variable:** and at least 2 groups. Adjust the font, **Panels per row**, points, mean and median lines, and group colors.
2. Tests run on the log2 data: Welch's t-test for 2 groups, one-way ANOVA for 3 or more. The adjusted p-value here covers only the selected metabolites; for panel-wide values, use Statistics.
3. Download the **Statistical Results** CSV and the figure.

## Biomarker Discovery

### Biomarker Discovery

This page needs a two-group result.

1. Choose a **Filtering criterion**: **p<0.05 AND FDR<0.25** (default), **p<0.05 only**, or **FDR<0.25 only**. For untargeted data, "FDR" reads "BH P Value". You can adjust both cutoffs.
2. Click **Discover Biomarkers**. The table shows p-value, adjusted p, Linear_FC, Log2FC and Direction (Up/Down).
3. Download `<A>_vs_<B>_Biomarkers.csv`.

## Pathway Analysis

Both pages reuse results you already ran on the Statistics page (two-group or ANOVA).

**Input Metabolite Selection**
1. Choose **Statistical results to use**.
2. Under **Select metabolites**, choose **Significant by P-value** (≤ 0.05), **Significant by FDR** (≤ 0.25) or **All metabolites (from t-test/ANOVA)**.
3. For two-group results, choose a **Direction** (both, up in B, or down in B).

### Metabolite-Set Enrichment Analysis (MSEA)

1. **ID Mapping**: metabolites are matched to HMDB IDs, then to SMPDB and KEGG pathway IDs.
2. **Reference Metabolome**: choose a **Metabolite-set library**: KEGG, SMPDB, LIPID MAPS, or **User-defined / custom metabolite sets** (CSV of Set_Name and Metabolite, or a .gmt file). By default the background is the whole library. Tick **Restrict reference to detected/measured metabolites only** for a more conservative test.
3. **Enrichment Analysis**:
   - **Statistical test**: Hypergeometric or Fisher's exact, both one-sided.
   - **Min set size**: default 2.
   - **Max set size**: 0 means no limit.

   Holm and FDR corrections are applied.
4. **Results** include a pathway table, an **Enrichment Dot Plot**, and a **Metabolite-Level UMAP** in which each dot is a metabolite colored by its enriched pathway. The UMAP needs about 15 or more metabolites.
5. Download **ID mapping (CSV)**, **Pathway results (CSV)**, **Complete results (Excel)** (which includes the parameters), and the UMAP coordinates.

### Metabolite–Gene Pathway Analysis (MGPA)

1. **HMDB Mapping** links selected metabolites to HMDB-associated genes. Metabolites without genes are listed and left out.
2. **Metabolite–Gene Mapping & Unique Gene List** shows the gene list that will be tested.
3. **Gene-Set ORA Settings**:
   - **Organism**: Human, Mouse, Rat, Yeast or Zebrafish.
   - **Gene-set database**: GO, KEGG, Reactome, WikiPathways or MSigDB (Hallmark). Hallmark works offline; the others need internet access (Enrichr).
   - **Background / universe**.
   - Gene-set size limits.
4. Click **Run Over-Representation Analysis**, then optionally **Generate Pathway Dot Plot**.
5. Download the results, analysis metadata, **HMDB mapping (CSV)**, **Metabolite–Gene mapping (CSV)**, **Unique gene list (CSV)** and **Complete results (Excel)**.

## Network Analysis

### Correlation Rewiring Map

This page finds metabolite pairs whose correlation changes between two groups: a coupling can be lost, gained or flipped. Before you have normalized data, you can try it with **Explore with the simulated rewiring demo (Control vs Tumor, 30 + 30 samples)**.

1. Choose a **Grouping variable**, a **Reference group (A)** and a **Comparison group (B)**. Each group needs at least 3 samples. With 4 or fewer, the test is exact and a note gives the smallest p-value that can be reached.
2. Main settings:
   - **Pathway library**: KEGG, SMPDB, Human-GEM subsystems or LIPID MAPS.
   - **Correlation**: Pearson or Spearman.
   - **Pairs to test**: all pairs, reaction-linked (≤ 2 steps), or direct reactions only.
   - **Rewired-pair FDR q-value cutoff (≤)**: default 0.05.
   - **Minimum |Δr|**: default 0.45.
3. **Advanced settings**:
   - features analysed (default 300)
   - missingness filter (default 0.3)
   - permutations (default 500) and bootstrap resamples (default 200)
   - median-centring and random seed
   - **FDR method**: Permutation FDR (recommended) or Benjamini–Hochberg
4. Click **Run Rewiring Analysis** to get a summary, an interactive network and hypothesis cards.
5. Downloads:
   - **Rewired pairs (CSV)** and **All pairs (CSV)**
   - **Hypothesis evidence (JSON)**
   - **Interactive report (HTML)**
   - **Pathway-level tests (CSV)**, **Reaction pairs (CSV)** and **Metabolite ID / pathway mapping (CSV)**

   The **How this works** expander explains the method.

## Multi-dataset mode

1. On **Data Upload**, choose **Multiple datasets (combine methods/modes)**. For the demo, choose **Demo Data** and click **Load Demo Datasets**.
2. For **Real Data**:
   - Set **Number of datasets** (1–20).
   - For each dataset, give it a unique label, choose its analysis type and upload its matrix.
   - Upload one shared metadata table covering all samples.
   - Click **Load All Datasets**. A warning flags any dataset whose samples are missing from the metadata.
3. Process each dataset separately, using **Select a dataset/method to configure** on the cleaning, QC and normalization pages. Or click **Auto-Process All Datasets (recommended defaults)**, which runs the following on every dataset:
   - a 30% missingness filter
   - Half-Minimum imputation
   - QC CV calculation
   - normalization: ISTD then Log2 for targeted datasets, Log2 then sample-wise Median-IQR for untargeted ones

   Details of the auto-processing: each targeted dataset is normalized with its ISTD (the one named in the row annotations' `ISTD` column, or the row named `ISTD…`); a targeted dataset without an ISTD stops with a message. Untargeted datasets leave the ISTD row out. The CV > 20% filter applies to targeted datasets only and never removes the ISTD.
4. **Combined QC (Multi-Dataset)** compares QC across datasets. Its **Apply to All Eligible Datasets** button sets the CV filter on all QC-confirmed datasets at once. Re-normalize afterwards.
5. On **Combined Normalization Data (Multi-Dataset)**, click **Generate Combined Normalized Data**. The button stays disabled until every dataset is normalized.
   - Features are prefixed by dataset (for example `Targeted Metabolomics::Glucose`), so the same compound measured by different methods is never merged.
   - Downloads: **Combined_Normalized_NoLog2.csv** and **Combined_Normalized_Strict_Log2.csv**.
   - All later pages then use this combined table. On screen the prefix is hidden, and repeated names get " (2)", " (3)" and so on.

## Tips & troubleshooting

- **"Metadata missing required column(s)"**: add the **Sample** and **Group** columns.
- **"Samples present in peak matrix but missing from metadata"**: every matrix column after the first must be listed under **Sample**.
- **"Duplicate metabolite names found"**: make the names in the first column unique.
- **QC module skipped**: at least 2 samples must be flagged in **IsQC**.
- **Log2 warns about non-positive values**: zeros became missing. Filter and impute on the cleaning page first.
- **Targeted Log2 is blocked**: apply ISTD normalization first.
- **"No internal standard was detected"**: name the ISTD row `ISTD_…` or list it in an `ISTD` column of the row annotations, or pick it by hand.
- **Volcano, Heatmap or Biomarker pages are empty**: run a **Two-group comparison**; an ANOVA alone does not feed these pages.
- **Heatmap says too few significant features**: clustering needs 2 or more, so relax the cutoff.
- **Pathway pages say no metabolites pass**: relax the threshold. If MGPA finds no genes, add an `HMDB_ID` column to your row annotations.
- **Multi-dataset downstream pages are empty**: normalize every dataset, then click **Generate Combined Normalized Data**.

## Glossary

- **Adjusted p-value (FDR q-value / BH P Value):** the p-value corrected for testing many metabolites at once (Benjamini–Hochberg). It is labelled "BH P Value" for untargeted data and "FDR" for targeted data; the calculation is the same.
- **Log2 FC:** mean log2 in Group B minus Group A. +1 means twice as high in B. **Linear_FC** = 2^Log2FC.
- **CV (coefficient of variation):** SD ÷ mean × 100 across QC injections. The app uses 20% as its cutoff.
- **QC sample:** a repeated pooled-QC injection, flagged in **IsQC** and excluded from analysis.
- **ISTD (internal standard):** a spiked-in reference compound. Dividing by it corrects for differences between injections.
- **Median-IQR normalization:** subtract the median, then divide by the interquartile range.
- **Strict log2:** log2(x) with no added constant.
- **Z-score:** standard deviations from a metabolite's mean. It sets the heatmap colors.
- **ORA / MSEA:** tests whether a pathway holds more of your selected metabolites (or their genes) than expected by chance.
- **Δr:** the change in a metabolite pair's correlation between two groups.
