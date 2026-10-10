# Omnix™ Fluxomics 2.5 — stable-isotope metabolomics and flux analysis

A standalone, step-by-step application (Python · Streamlit) for stable-isotope tracing experiments measured by
LC-MS. Everything runs on your computer; data are not uploaded anywhere.

## 1 · Install (Anaconda, once)

Open **Anaconda Prompt** (Windows) or a terminal (macOS/Linux) in this folder:

```bash
conda env create -f environment.yml      # creates the "fluxomics" environment (Python 3.11)
conda activate fluxomics
```

Or, into an existing Anaconda environment: `pip install -r requirements.txt`.

## 2 · Start

* Windows: double-click **Start_App.bat**  ·  macOS / Linux: `bash start_app.sh`
* or: `conda activate fluxomics` then `streamlit run app.py` → opens at http://localhost:8501

## 3 · Input

* **Compound Discoverer SIL exports** (“Relative Area 0 … 25” = M+0 … M+25), one file per method (e.g. RP_Pos +
  HILIC_Neg). Peaks are matched across samples by retention time and each metabolite is kept once (best peak,
  name and method). Sample names are parsed into group / time / replicate.
* **Generic isotopologue table** (any software): `Metabolite`, `Formula`, `Method`, optional `KEGG` `HMDB`
  `PubChemCID` `ChEBI` `BiGG`, `Isotopologue` (M+0, M+1 …), one column per sample — or long format
  (`Metabolite`, `Isotopologue`, `Sample`, `Intensity`). Plus a **metadata** file: `Sample`, `Group`, `Time`,
  `TimeUnit`, `Treatment`, `Replicate`, `Type` (Labeled / Unlabeled reference / Blank / QC), `Include`.
  Templates can be downloaded in Step 1. Each isotopologue must appear once per metabolite and sample.

**Zero values in labeled samples are not valid data** (rule of `sil_flux_pipeline.py`, unchanged): an
isotopologue row that is 0 in any labeled sample is removed, and a compound whose M+0 is 0 in any labeled sample is
removed with all its isotopologues. The analysis uses exactly the rows of the Step 1 Inspect table. Empty cells
(compound not reported in a sample) are not zeros; they stay missing (no imputation).

**Isotopologues are never summed.** Every M+0 … M+n is its own channel in normalization, statistics, bar plots,
heatmaps, maps and the label-arrival figure. The only sum is the denominator of an isotopologue *fraction*
(M+i ÷ Σ M+0…M+n of the same metabolite in the same sample — every isotopologue, as Compound Discoverer's
"Exchange Rate in Percent"; the zero / minimum-value rule only filters the intensity analyses).

## 4 · Workflow (10 steps)

| Step | What it does |
|---|---|
| 1 Upload & inspect | CD exports → combined table built exactly as `sil_flux_pipeline.py` stage 1, plus the **minimum-value rule** (< 1,000 in a labeled sample), then metabolites with **only M+0** are excluded. Download the combined workbook. Generic tables also supported. |
| 2 Experimental design | Groups, treatments, time points, replicates, include/exclude; tracer and natural-abundance correction (user choice; IsoCor high / low resolution, AccuCor resolution-dependent, classical least-squares). |
| 3 Normalization | log2, then robust IQR per labeled sample (median / IQR over all isotopologue channels; 12C reference and blank excluded). The log2 + IQR data are used by every later intensity analysis. |
| 4 Isotope incorporation | Isotopologue fractions, incorporation, FE, mean label number per sample and group × time. |
| 5 Dynamic bar plots | One plot per metabolite: bars = mean per group × time, error bar SD / SEM, dots = replicates. |
| 6 Statistics | t-tests and ANOVA per isotopologue channel on log2 + IQR data; p and FDR for every channel. |
| 7 Pathway mapping | 25 curated pathway maps, one at a time. Every detected metabolite is a card with a Z-score heatmap (row-wise Z of each isotopologue channel from the same log2 / log2 + IQR data as the Heatmap step): rows M+0…M+n; columns = groups (coloured bar) × time points (dot colour). Identifiers of the shown pathway are verified online (KEGG / PubChem) automatically. Undetected metabolites in italics; other detected pathway members listed below the map. |
| 8 Heatmaps | Samples × isotopologue channels, row-wise Z-score across the selected conditions × time points; entire dataset or selected pathways (combined or one heatmap per pathway); rows by FDR (default ≤ 0.25), P value (default ≤ 0.05), both, variance or all; non-significant rows hidden or shown marked “(n.s.)”; row labels “Metabolite | M+i” or grouped by metabolite (M+i labels, metabolite name joined to its rows by a vertical line); rows separated into pathway blocks with an adjustable gap and vertical pathway names on the left; balanced proportions by default, fixed cells or custom size; the on-screen heatmap is centred and shrinks to fit narrow windows. |
| 9 Label arrival | All groups and time points side by side; up to 100 isotopologues per pathway; threshold from 0 %. |
| 10 Export & report | Excel workbook, CSV bundle, collected figures, HTML report, methods paragraph, settings file. |

**Downloads.** Every figure and table has the same “⬇ Downloads” box: figure format (PNG / TIFF / JPEG raster,
PDF / SVG vector), resolution 300 / 600 (default) / 1200 DPI, width × height → Render figure → download; each table
as CSV or Excel, plus one multi-sheet workbook where a view has several tables. The camera icon on interactive
charts saves a 5× PNG.

### Changes in 2.5.11
- **Templates included** (`example_data/templates/`), downloadable in Step 1:
  - *Compound Discoverer SIL:* `CD_SIL_Template_RP_Pos.xlsx` and `CD_SIL_Template_HILIC_Neg.xlsx` — the SIL export
    on the first sheet (columns Name … Relative Area 0–25, raw files named
    `<date>_<method>_CellLine_<Group>_<SampleNo>_<Time>.raw`), plus *Sample_List* and *Instructions* sheets.
    New button **Load Compound Discoverer example** processes both templates like uploaded exports
    (shControl vs GeneKD, 12 / 18 / 40 h, 3 replicates, one 12C reference, one blank per method).
  - *Generic format:* isotopologue table with and without the 12C reference column, sample metadata (with
    Replicate and Include) and an instructions workbook. Same layout as the synthetic example.

### Changes in 2.5.10
- **Isotopologue distribution chart (Step 4): sample labels never overlap.** On screen the x axis has two levels —
  the time point under each bar and the group name centred under its bars (no slanted, clipped names); the legend
  sits in its own column on the right. When a group name is too long for its bars (e.g. many single-sample groups),
  every bar gets its own vertical “group · time” label instead. The exported figure uses the same rule.

### Changes in 2.5.9
- **Pathway blocks (Step 8):** with *Separate rows by pathway* (on by default for *Selected pathways → One combined
  heatmap*; optional for *Entire dataset*, using each metabolite's first assigned pathway plus “Other / not
  mapped”), every pathway is its own block of rows — pathways never mix, and clustering orders rows only within a
  pathway. A metabolite in two selected pathways appears in both blocks.
- **Gap between pathways:** adjustable blank space between blocks (in row heights, default 0.5; 0 = none).
- **Pathway labels on the left:** a vertical line per block with the pathway name written vertically beside it,
  mirroring the metabolite labels on the right; long names wrap (or shrink slightly) to fit their block. Font size
  under *Figure size and fonts → Pathway names*. Exported tables gain a *Pathway* column.
- **Always centred:** the heatmap stays in the middle of the page whatever conditions / time points are chosen
  (also when it shrinks on a narrow window); in the figure, extra width from a long title is split evenly on both
  sides so the heatmap block stays centred under the title.

### Changes in 2.5.8
- **Heatmap shrinks with the window (Step 8):** on a narrow or minimised page the heatmap fits the width instead of
  being cut off (cells get narrower; row labels, metabolite names, legend and colour scale keep their size and
  position). The exported figure keeps its own size. On a very narrow window the long title can still be clipped at
  its ends.
- **Non-significant rows:** with *Select rows by → Significance*, the new choice *Non-significant rows (n.s.)* under
  Criterion: *Hide (significant only)* (as before) or *Show, marked “(n.s.)”* — every channel is shown and those not
  passing the threshold in any selected comparison are labelled e.g. “M+2 (n.s.)”, on screen and in the export.
- **Row label style** (under *Column annotations, clustering and colours*): *Metabolite | M+i* (default) or
  *Grouped by metabolite* — rows of a metabolite stay together (M+0, M+1 …) with the metabolite name in its own
  column, joined to its rows by a plain vertical line (no overlap with the M+ labels). Metabolites follow the
  clustering order; no row dendrogram in this mode.
  The line sits a short, fixed distance right of the widest M+ label (measured from the font, on screen and in
  the export); highlighted rows are shown as bold red labels in this mode.

### Changes in 2.5.7
- **Step order:** 7 · Pathway mapping, 8 · Heatmaps.
- **Pathway mapping (Step 7):** Z-score heatmaps only (isotope-incorporation panels removed; incorporation stays in
  Steps 4, 5 and 9). Metabolite IDs of the selected pathway are verified online (KEGG / PubChem) automatically, once
  per metabolite and session; offline it is skipped with a retry button. Failed lookups are no longer cached.
- **Heatmaps (Step 8):** choose any conditions and time points (e.g. Control vs KD at 12 & 18 h); the Z-score is
  recomputed across the samples shown. *Entire dataset* or *Selected pathways* (from the annotation), as one combined
  heatmap or one heatmap per pathway, each with its own downloads.
- **Grouped row labels (default):** isotopologues of a metabolite stay together (M+0 … M+n), labelled “M+i”, with
  “(n.s.)” for channels not significant in the selected Step 6 result, a bracket and the metabolite name on the
  right, and thin lines between metabolites. Row clustering then orders whole metabolites (mean channel profile).
- **Balanced proportions (default):** cell width and height follow the number of rows and samples, so the heatmap
  body keeps publication proportions; *Fixed cell size* and *Custom figure size* remain available.

### Changes in 2.5.6
- **Bar plot export (Step 5) rebuilt.** Fixed layout in inches: bold metabolite title and subtitle centred over the
  plot, legend in its own column on the right, axis labels measured — the plot is never squeezed by the title or the
  legend. Default size gives every bar ≈ 0.15 in.
- **TIFF downloads** (all figures) are now flat RGB with the DPI tag set (no transparent alpha channel), so they open
  correctly in Word, PowerPoint, Prism, Illustrator and journal upload systems.
- **Heatmap (Step 7): landscape proportions.** Separate *Cell width* (default 0.30 in) and *Cell height* (default
  0.10 in): the heatmap is much wider and much shorter; row labels are sized to the row (≈ 6 pt) and stay clear.

### Changes in 2.5.5
- **Heatmap (Step 7): smaller, square cells.** New *Cell size* setting (default 0.12 in ≈ 3 mm, range 0.08–0.30 in);
  the figure size follows from the number of rows and samples, on screen and in the export. *Custom figure size*
  still lets you type an exact width/height. Every row and sample label is always shown (no automatic thinning), and
  the colour scale no longer shrinks the heatmap.
- **Pathway mapping (Step 8): compact heatmap cards.** Cells are about 30 % narrower and 25 % shorter with unchanged
  text sizes; panel titles never overlap. *Fit* no longer blows small maps up beyond 160 %.

### Changes in 2.5.4
* Fixed: dynamic bar plots crashed with “No error bar” for metabolites with more than one isotopologue.
* Very tall or wide figures: raster exports (PNG / TIFF / JPEG) are kept within 60 million pixels by lowering the DPI
  automatically (a note says so); PDF / SVG stay full vector quality. Prevents out-of-memory crashes on small servers
  (e.g. Streamlit Cloud) for “All channels” heatmaps.

### Changes in 2.5.3
* Isotopologue fractions use Compound Discoverer's “Exchange Rate in Percent” formula: Area(M+i) / Σ Area(M+0…M+n)
  over EVERY isotopologue of the metabolite (also those removed by the zero / minimum-value rule). Incorporation
  (100 − M+0), mean label number (Σ i·fraction) and fractional enrichment (mean label number / n) follow from it, and
  every labeling view (Step 4, bar plots, pathway maps, label arrival, exports) shows all isotopologues M+0…M+n, so
  distributions sum to 100 %. Intensity analyses (normalization, statistics, heatmaps, Z-scores) are unchanged.

### Changes in 2.5.2
* Example data replaced by a fully synthetic data set (generic format; shControl vs GeneKD). The reference raw-file
  naming pattern now recognises `…_shControl_<n>_<time>.raw` / `…_GeneKD_<n>_<time>.raw`; other names use the generic
  file-name parser.
* Heatmap colour scale: shorter and thicker bar, title on its own line above it (no overlap).

### Changes in 2.5.1
* Sidebar intensity dataset: Log2 intensity listed first, then Log2 + IQR-normalized.
* Normalization: separate figure downloads for the log2 plot and the log2 + IQR plot (plus the combined figure).
* Dynamic bar plots: three value types — isotope incorporation (%), isotopologue fraction (%), intensity (sidebar).
* Heatmap: centred title; colour scale bottom-right directly under the heatmap; ≈ half the previous width
  (0.19 in per sample) and thicker rows (0.24 in); black text.
* Pathway maps: shown in a vector viewer (drag / wheel / + − / Fit) — names, values and labels scale with the zoom;
  hover a cell for its value. Metabolite name on top of its incorporation + Z-score heatmaps; spacing tightened;
  larger, clearer cofactor, link and label text; smaller default export (6.5 pt per heatmap row).

### Changes in 2.5
* Pathway mapping redrawn from curated pathway layouts (the 10 central-carbon templates — glycolysis, PPP, TCA,
  one-carbon, purine and pyrimidine synthesis, glutathione, glutamine, fatty acids, amino acids → TCA — plus 15 more
  in the same style: serine/glycine, urea cycle & polyamines, creatine, hexosamine & nucleotide sugars, Kennedy
  pathway, glycerolipids, carnitine shuttle, ketone bodies, BCAA catabolism, GABA shunt & NAA, NAD⁺ & kynurenine,
  purine and pyrimidine salvage/degradation, phenylalanine/tyrosine, histidine). Layouts expand automatically so
  every metabolite carries its incorporation + Z-score card without overlaps; arrows, cofactors, links to other
  pathways, sub-panels and the TCA ring follow the expansion. Metabolites matched by KEGG ID.
* Old map code removed (built-in central-carbon network, Escher / KEGG import, card-map modules).

### Changes in 2.4.2
* Natural-isotope correction (Step 2): applied only when you choose it (data not yet corrected). Four methods —
  IsoCor high-resolution (tracer element only), IsoCor / IsoCorrectoR low-resolution (all elements), AccuCor
  resolution-dependent (Orbitrap / FT-ICR / TOF, resolving power), classical skewed-matrix least squares.
* Isotope incorporation figure: title, legend, time labels and group brackets in separate bands (no overlaps).
* Heatmap: FDR (default 0.25) and P value (default 0.05) thresholds, alone or combined; width / height, font family
  and font sizes; tighter legend spacing; colour scale underneath the heatmap.
* Pathway mapping: one pathway at a time, isotope incorporation + Z-score heatmap on every metabolite; Escher and
  KEGG (KGML) map options removed; vertical colour scales at the bottom with selectable colours.
* "Label arrival (race)" renamed to "Label arrival"; all text solid black.

### Changes in 2.4.1
* Sidebar box (highlighted) **Intensity dataset for downstream steps**: Log2 + IQR-normalized intensity (default)
  or Log2 intensity — used by bar plots, statistics, heatmaps and pathway abundance.

### Changes in 2.4
* Label arrival: up to 100 isotopologues per pathway, threshold default 0 %; publication layout measured in inches (label
  column, title, group headers, time labels, axis label and legend in separate bands; ≥ 0.15 in per row).
* Upload & inspect: Advanced options removed. Normalization: reference-channel choice and sample-amount factors
  removed (median / IQR over all channels).
* Heatmap: thresholds up to 1; only row-wise Z-score (no clipping); sample columns only.
* Pathway map: all time points side by side, distinguished by dot colour; single-condition node dot removed.
* Interpretation (AI) step removed.
* Unified high-resolution download box for all figures and tables (default 600 DPI).

### Changes in 2.3
* Step 1: new **minimum-value rule** (default 1,000, editable): an isotopologue with a value below it in any labeled
  sample is removed; if M+0 is below it, the whole peak is removed (12C reference / blank not checked; empty cells
  do not count). Runs before the one-peak-per-metabolite selection.
* Label arrival is a static side-by-side figure (no Play / animation anywhere in the app).
* Heatmap: annotation legend (group, time …) on the right, vertically centred; small colour scale at the bottom-right.
* Statistics: no p-value threshold or "significant" call — p and FDR listed for every channel, sorted by p.
* Fixed: isotope-incorporation stacked bars were misaligned for metabolites with missing isotopologues
  (e.g. glucose M+0, M+1, M+6); stacks now follow M+0 → M+n and sum to 100 %.

### Changes in 2.2
* Combined table now matches `sil_flux_pipeline.py` exactly (v2.1 differed: extra library-synonym merging joined
  *Glutathione* with *L-Glutathione (reduced)*, and a generic file-name parser renamed samples `shControl_…` instead
  of `shCtrl_…`). New rule: metabolites with only M+0 are excluded.
* Fixed a colour bug in pathway maps: values above the colour maximum were drawn almost white instead of darkest.
* Removed: box plots, multiple-testing selector, "All pathways (grouped)" and "Escher viewer (online)" tabs,
  labeling curves / slopes / first-order fits / speed ranking / static label-arrival figure, and the "Flow map · isotope clock ·
  fingerprint glyphs" step.

## 5 · Validation

`python -m pytest -q tests` — 62 checks: Welch / ANOVA / Tukey / two-way ANOVA / p-value adjustment against scipy and
statsmodels, Games–Howell and Welch ANOVA against pingouin (if installed), the zero rule, t-tests / Kruskal–Wallis against scipy, robust IQR normalization, the no-summing rule for generic tables,
natural-abundance correction, labeling metrics, kinetic fits, contrast builder, correction methods,
the overlap-free layout of every curated pathway map, KEGG binding of map nodes, the M+0-only rule, the reference sample naming, and the full pipeline on the example
data (identical to `sil_flux_pipeline.py`). Tested on Python 3.11 (pandas 1.5, numpy 1.24,
Streamlit 1.35, without scikit-learn) and Python 3.13 (pandas 3, numpy 2.5, Streamlit 1.65).

## 6 · Notes

* **Labeling is not flux.** Fractions, FE and the label-arrival figure describe where and how fast label appears,
  not absolute pathway flux.
* **Annotation library** (`data/metabolite_library.csv`) was curated by hand — verified online automatically in Step 7 (Pathway mapping)
  (annotation panel below the map; queries KEGG REST / PubChem from your computer).
* **Pathway maps** (`core/pathway_templates.py`) are curated drawings; metabolites are matched by KEGG ID (name as
  fallback). Add a pathway by adding a template there (same builder calls: met, stub, edge, arc, circle, box).
* **Example data** (`example_data/generic/`) are fully SYNTHETIC: shControl vs GeneKD, [U-13C6]-glucose, 2 / 8 / 24 h,
  3 replicates, one unlabeled reference and one blank; 63 metabolites (names and formulas from the public metabolite
  library), every value simulated with a seeded labeling model (`make_synthetic_example.py` regenerates them). They
  come from no real experiment. Load them in Step 1 → Generic isotopologue table + metadata → Load synthetic example.
* **Folder layout:** `core/` = analysis engine (no UI; importable from scripts or notebooks), `ui/` = one module per
  step, `app.py` = entry point.
