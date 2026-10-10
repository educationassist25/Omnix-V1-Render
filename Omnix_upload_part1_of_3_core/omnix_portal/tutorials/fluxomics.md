## Overview

Fluxomics analyses stable-isotope tracing experiments measured by LC-MS (for example cells grown on [U-¹³C₆]-glucose or [U-¹³C₅]-glutamine and sampled at several time points). It takes the isotopologue peak areas of every metabolite through experimental design, natural-abundance correction, normalization, isotope incorporation, statistics, bar plots, heatmaps, pathway maps and a label-arrival figure, and ends with a complete export and report.

Two rules hold throughout the platform:

- **Isotopologues are never summed.** Every M+0 … M+n of a metabolite is its own channel in normalization, statistics, bar plots, heatmaps and maps. The only sum is the denominator of an isotopologue *fraction*.
- **Labeling is not flux.** Fractions, fractional enrichment and label arrival describe where and how fast label appears, not absolute pathway flux.

The sidebar lists seven workflow groups: **Setup**, **Normalization**, **Isotope Labeling**, **Statistical Analysis**, **Exploratory Analysis**, **Pathway Analysis** and **Export & Report**. Clicking a group shows its steps as buttons at the top of the main area; the steps keep their numbers (Step 1 to Step 10), and a ✓ marks a step that has been completed. Work through them in order: each step uses earlier results, and a step that needs an earlier one says so and offers a **Go to …** button. Under the groups, the sidebar shows the loaded **Dataset** and the **Intensity dataset for downstream steps** switch (log2 + IQR-normalized, the default, or log2), which every intensity analysis uses.

Every figure and table has the same **⬇ Downloads** box: figure format (PNG, TIFF, JPEG, PDF or SVG), resolution 300 / 600 (default) / 1200 DPI and width × height in inches → **Render figure** → download; tables as CSV or Excel.

## Before you start: input files

Fluxomics accepts two kinds of input. Choose one in Step 1 under **Input format**.

### Option A: Compound Discoverer SIL export(s)

- One export per LC-MS method (Excel or CSV), e.g. `…_RP_Pos.xlsx` and `…_HILIC_Neg.xlsx`. The workbook name must contain the method; the data must be on the **first sheet**.
- Required columns (exact names, any order): **Name**, **Formula**, **Molecular Weight**, **RT in min**, **Area**, **Status**, **Study File ID**, **Filename**, and **Relative Area 0** … **Relative Area 25** (Relative Area i = isotopologue M+i).
- Raw-file names carry the design: `<date>_<method>_<CellLine>_<Group>_<SampleNo>_<Time>.raw`, e.g. `08302026_RP_Pos_CellLine_GeneKD_1001_12Hr.raw`. Time `12C` marks the unlabeled reference, and a name containing `Blank` marks a blank. Other group names also work; check what was read in Step 2.

### Option B: Generic isotopologue table + sample metadata

**Isotopologue table (required)**
- CSV, TSV or Excel (first sheet). One row per metabolite isotopologue, one column per sample.
- **Metabolite**, **Formula** (recommended: sets n, the number of carbons and so the maximum M+n), **Method**, optional identifiers **KEGG**, **HMDB**, **PubChemCID**, **ChEBI**, **BiGG** (KEGG IDs improve pathway mapping), **Isotopologue** (M+0, M+1 … or 0, 1 …), then one column of peak areas per sample.
- A long table also works: **Metabolite**, **Isotopologue**, **Sample**, **Intensity**.
- Each isotopologue must appear once per metabolite and sample. A missing row means not detected; an empty cell is a missing value (it is not imputed).

**Sample metadata (optional but recommended)**
- **Sample** (exactly the sample column names of the table), **Type** (Labeled, Unlabeled reference, Blank or QC), **Group**, **Treatment** (optional), **Time** (a number; empty for the reference and blank), **TimeUnit** (min, h or d), **Replicate** (biological replicate number), **Include** (TRUE / FALSE; FALSE leaves a sample out of every analysis).

### Rules applied to the data

- **Zero and low values in labeled samples are not valid data.** An isotopologue with a value below the **minimum value** (default 1,000; 0 = off) in any labeled sample is removed, and a compound whose M+0 is below it is removed with all its isotopologues. The unlabeled reference and blank are not checked; empty cells do not count.
- With Compound Discoverer exports, peaks are matched across samples by retention time and each metabolite is kept **once** (best peak, name and method). Metabolites left with only M+0 are excluded.

### Demo datasets

- **Load synthetic example data** (Generic route): a fully simulated study, shControl vs GeneKD on [U-¹³C₆]-glucose at 2 / 8 / 24 h with 3 replicates, one unlabeled reference and one blank; 63 metabolites. The values come from no real experiment.
- **Load Compound Discoverer example** (Compound Discoverer route): the two template exports (RP_Pos and HILIC_Neg), processed exactly like uploaded files: shControl vs GeneKD at 12 / 18 / 40 h, 3 replicates, one 12C reference and one blank per method.

**Templates and example files (Help → Data preparation, and Step 1)**
- Generic: isotopologue table with and without the 12C reference column, sample metadata, and an instructions workbook.
- Compound Discoverer: the RP_Pos and HILIC_Neg template exports, each with a *Sample_List* and an *Instructions* sheet.

## Setup

### 1 · Upload & inspect

1. Choose the **Input format**: **Compound Discoverer SIL export(s)** or **Generic isotopologue table + metadata**.
2. Load a demo dataset, or upload your files (with a subscription) and click **Process files** (Compound Discoverer) or **Read table** (Generic).
3. Set the **Minimum value in labeled samples** if needed (default 1,000) before processing.
4. Inspect the result: metrics (source, metabolites, isotopologue channels, samples, labeled samples), the **Processing log**, and the tabs **Isotopologue channels**, **Samples** and **Missing values** (per sample and per isotopologue). Compound Discoverer data add **Combined table (clean)**, **Selection decisions** and **Removed peaks**.
5. Downloads: the combined workbook (Compound Discoverer), the channel table in the generic layout and the sample metadata.

No values are changed in this step; the analysis uses exactly the rows of the Inspect table.

### 2 · Experimental design

1. **Samples:** edit Group, Treatment, Time and Replicate if needed; untick **Include** to leave a sample out. When Treatment is given, each group × treatment combination becomes its own group. A check reports whether the design looks complete (≥ 3 replicates per group and time, ≥ 3 time points, an unlabeled reference).
2. **Tracer:** [U-¹³C₆]-glucose or [U-¹³C₅]-glutamine, and the isotopic purity. The tracer enrichment in the medium is recorded in the report (labeling values are not rescaled by it).
3. **Natural-abundance correction:** *Data already corrected by the vendor software*, *Apply natural-abundance correction here* or *No correction*. When an unlabeled reference is present, the app checks automatically whether the data look corrected. Methods: IsoCor high-resolution (tracer element only), IsoCor / IsoCorrectoR low-resolution (all elements), AccuCor resolution-dependent (Orbitrap, FT-ICR or Q-TOF, with resolving power), or the classical skewed-matrix least-squares inversion. Correction applies to isotopologue fractions and labeling metrics; intensity analyses use the measured areas.
4. Click **Confirm design**.

## Normalization

### 3 · Normalization

Runs automatically the first time you open this step after confirming the design (**Re-run log2 + IQR normalization** repeats it).

1. **log2:** L = log2(area). Values are not imputed; an empty cell stays missing.
2. **Robust IQR:** each labeled sample is centred on its median and scaled by its interquartile range over all isotopologue channels, so every sample gets the same median and IQR. Values stay in log2 units, so group differences remain log2 fold changes.
3. The 12C reference and blank never enter the median / IQR and are not part of the normalized data.

The page shows the spread of sample medians and IQRs and per-sample box plots before and after normalization, with downloads for the log2 plot, the log2 + IQR plot, the combined figure, the scale factors and both matrices. The sidebar switch **Intensity dataset for downstream steps** chooses which dataset later intensity analyses use.

## Isotope Labeling

### 4 · Isotope incorporation

Where the label went, for every metabolite:
- **Fraction M+i (%)** = area(M+i) ÷ Σ area(M+0 … M+n) × 100, over every isotopologue of the metabolite (as Compound Discoverer's "Exchange Rate in Percent").
- **Incorporation (%)** = 100 − fraction M+0.
- **Fractional enrichment (%)** = Σ i × fraction(M+i) ÷ n × 100.
- **Mean label number** = Σ i × fraction(M+i).

Choose a **Metabolite** to see its isotopologue distribution per sample (toggle **Hide M+0** to zoom into the labeled part) and its mean fractions and metrics per group × time. The **Isotope incorporation table** lists every metabolite, per sample or as group × time means (mean, SD, n), with a filter.

### 9 · Label arrival

Every labeled isotopologue (M+1 … M+n, never summed) of the metabolites in the chosen pathways, for every group and time point side by side, as a static publication figure. Settings: pathways, groups, time points, **Show isotopologues reaching ≥ (%)** (from 0 %) and **Max isotopologues per pathway** (up to 100). Needs at least two time points.

## Statistical Analysis

### 6 · Statistics

Tests run on every isotopologue channel separately, on the active intensity dataset (sidebar). Effects are log2 fold changes (B − A).

- **Two-Group (T-Test):** Welch's t-test, Student's t-test, paired t-test (pairs matched by replicate number) or Mann–Whitney U, at each time point separately.
- **ANOVA:** one-way ANOVA (post hoc Tukey–Kramer or pairwise Welch), Welch's ANOVA (Games–Howell or pairwise Welch), Kruskal–Wallis (pairwise Mann–Whitney), or two-way ANOVA (group × time, Type II sums of squares), comparing groups at each time point or time points within each group.

Click **Run statistics**. Every result table lists the p value and the Benjamini–Hochberg FDR for each channel, sorted by p value; there is no fixed significance threshold. Several runs are kept and can be chosen under **Result set** and **Comparison**.

## Exploratory Analysis

### 5 · Dynamic bar plots

One plot per metabolite (up to 8 at a time): **Values** = isotope incorporation (%), every isotopologue fraction (%) or every isotopologue's intensity (active dataset). One bar per group × time point (group colour, lighter = earlier time), bar = mean of the biological replicates, error bar **SD**, **SEM** or none, dots = replicates (**Show replicates**). Choose the groups and time points, and **Hide M+0** if needed.

### 8 · Heatmaps

Clustered heatmaps of isotopologue channels (row-wise Z-score, one column per sample) from the active intensity dataset.

1. **Conditions, time points and pathways:** choose any conditions and time points (the Z-score is recomputed across the samples shown), and either the **Entire dataset** or **Selected pathways** (one combined heatmap or one per pathway). With **Separate rows by pathway**, each pathway is its own block of rows with its name on the left; the gap between blocks is adjustable.
2. **Rows:** by significance from Step 6 (FDR, default ≤ 0.25; P value, default ≤ 0.05; or both), by variance, or all. Non-significant rows can be hidden or shown marked "(n.s.)". **Highlight rows** marks chosen channels.
3. **Column annotations, clustering and colours:** sample annotation tracks, clustering, colour scale, and row labels as *Metabolite | M+i* or *Grouped by metabolite*.
4. **Proportions:** *Balanced* (default), fixed cell size or a custom figure size, plus font sizes for row labels, metabolite names and pathway names.

## Pathway Analysis

### 7 · Pathway mapping

One curated pathway map at a time (25 maps: glycolysis, pentose phosphate pathway, TCA cycle, one-carbon metabolism, purine and pyrimidine synthesis, glutathione, glutamine, fatty acids, amino acids, serine/glycine, urea cycle and polyamines, and more). Every detected metabolite is a card with a Z-score heatmap of its isotopologue channels (rows M+0 … M+n; columns = groups × time points). Metabolites are matched by KEGG ID (name as fallback); undetected metabolites are shown in italics, and other detected members of the pathway are listed below the map.

- The map is a vector viewer: drag to move, mouse wheel or + / − to zoom, **Fit** for the whole width.
- Identifiers of the shown pathway are verified online (KEGG / PubChem) automatically. If those services cannot be reached, verification is skipped and can be retried.
- **Metabolite annotation:** review and edit identifiers and pathways (edited rows become *User-confirmed*), accept low-confidence rows, run online verification, and download or restore the annotation table. The annotation also supplies the pathways used in Heatmaps and Label arrival.

## Export & Report

### 10 · Export & report

- **Prepare Excel workbook:** processed data and every statistics table in one workbook. **Prepare CSV bundle (ZIP)** gives the same as CSV files.
- **Collected figures:** every figure rendered in a ⬇ Downloads box of Steps 4–9 is collected here.
- **Build report:** an HTML report with input checksums, processing log, design, all settings, key result tables, collected figures and a methods paragraph to copy into a manuscript.
- **Settings file:** save the settings as JSON and load them later (**Apply settings**), then re-run Step 3 and the analysis steps.

## Tips & troubleshooting

- **"Complete … first"**: the step needs an earlier one; click the **Go to …** button. The usual order is Upload → Design → Normalization → the analysis steps.
- **Uploads are disabled**: without a subscription only the two demo datasets can be used. Subscribe or enter your access key on the Subscription page.
- **"Could not process the files" (Compound Discoverer)**: check the exact column names (Relative Area 0 … 25), that the data are on the first sheet, and that the workbook name contains the method.
- **Generic table errors**: each isotopologue must appear once per metabolite and sample, and the metadata's Sample column must match the column names exactly.
- **Few metabolites after Step 1**: the minimum-value rule removes isotopologues (or whole peaks) below the threshold in any labeled sample, and metabolites with only M+0 are excluded. Lower the minimum value (0 = off) if appropriate, and check the Processing log and Removed peaks.
- **Labeling in the unlabeled reference**: a metabolite that looks labeled in the 12C reference usually has an isobaric interference; interpret it with care.
- **Correcting twice**: if the vendor software has already corrected the data, choose *Data already corrected*; correcting again distorts the fractions.
- **Label arrival is empty**: it needs at least two time points, and the threshold or pathway selection may exclude every isotopologue.
- **Online ID verification skipped**: KEGG / PubChem could not be reached; the maps still work with the built-in library.
- **Reset all steps** (sidebar → Session) clears the Fluxomics analysis only; your Omnix access is kept.

## Glossary

- **Isotopologue (M+i)**: a molecule carrying i heavy atoms (e.g. ¹³C) of the tracer element; M+0 is unlabeled.
- **Isotopologue channel**: one isotopologue of one metabolite, analysed as its own variable.
- **Isotopologue fraction**: area of M+i divided by the sum over every isotopologue of the metabolite in the same sample.
- **Incorporation**: 100 % minus the M+0 fraction: the share of the pool that carries any label.
- **Fractional enrichment (FE)**: mean label number divided by the number of tracer-element atoms (n), as a percentage.
- **Mean label number**: Σ i × fraction(M+i), the average number of heavy atoms per molecule.
- **Natural-abundance correction**: removal of the heavy isotopes that occur naturally (about 1.1 % ¹³C per carbon) so that fractions reflect the tracer only.
- **Unlabeled (12C) reference**: a sample grown without tracer; used to check whether data are already corrected and to spot interferences.
- **Robust IQR normalization**: per-sample scaling to a common median and interquartile range on log2 data.
- **Z-score (heatmaps, maps)**: each channel's values scaled to mean 0 and SD 1 across the samples or conditions shown.
- **FDR**: the Benjamini–Hochberg adjusted p value, which controls the false discovery rate across tested channels.
- **Label arrival**: how quickly label appears in each isotopologue across time points.
