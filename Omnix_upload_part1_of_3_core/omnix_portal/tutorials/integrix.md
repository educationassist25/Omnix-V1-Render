## Overview

Integrix integrates transcriptomics, proteomics and metabolomics (or lipidomics) measured on the same samples. It follows the directional chain **gene → protein → enzymatic reaction → metabolite**: for every reaction it asks whether the change in its enzyme and the changes in its substrates and products point the same way, and calibrates each answer by permutation with false-discovery-rate (FDR) control. Results describe evidence *consistent with* a mechanism, not proof of causation.

The reaction, not the single molecule, is the unit of integration. Enzyme evidence comes from protein abundance through the gene–protein–reaction rules of the Human-GEM metabolic model (from mRNA, down-weighted, when proteomics is missing); the metabolite balance compares products with substrates. Their product, the **coherence**, is tested against sample-label and network-topology permutations.

Integrix works with three or two layers: **Transcriptomics + Proteomics + Metabolomics**, **Transcriptomics + Metabolomics** or **Proteomics + Metabolomics**. Paired samples give the strongest evidence (Tier 1); partly overlapping (Tier 2) and unpaired designs (Tier 3) are handled and labelled as such.

## Before you start: input files

Inside Omnix you normally need **no files**:

1. Analyse each dataset in its own platform (Transcriptomics, Proteomics, Metabolomics or Lipidomics) up to and including **Normalization**.
2. In that platform's side panel, press **Send to Integrix**. The normalized matrix, the sample metadata (QC injections removed) and the processing notes go to Integrix in this session.
3. Open Integrix and choose **Data → From this Omnix session**. Each dataset you sent is listed; tick the ones to integrate.

To keep the data for later, or to move it to another computer, press **Download for Integrix (.zip)** in each platform instead, and upload the zips in Integrix under **Data → Upload Omnix exports**. Each zip holds the matrix, the metadata and a small `*_omnix.json` file saying how the matrix was made (normalization, imputation, batch correction); Integrix uses the matrix as it is and cites those steps in its Methods.

You can also upload your own tables under **Upload my data**: one wide table per layer (identifier column(s) first, then one column per sample) and one sample metadata table, as **Raw intensities** (Integrix log2-transforms them) or **Normalized (log2)** (used as provided).

- **Identifiers:** gene symbols or Ensembl IDs (transcriptomics); UniProt accession plus a **Gene** column (proteomics); metabolite names plus **HMDB**, **KEGG** or **ChEBI** where possible (metabolomics).
- **Samples:** the same sample names in every matrix and in the metadata **Sample** column (capitals and surrounding spaces are ignored).

## Step 1: Data

The setup panel on the left has four steps; each collapses to a one-line summary when it is ready.

- **From this Omnix session**: the datasets sent from the platforms, with the number of features and samples and the time they were sent. Metabolomics and Lipidomics both fill the metabolite layer: tick the one to integrate. **Clear the sent data** empties the list.
- **Demo: Omnix normalized** and **Demo: Raw intensities**: a synthetic Condition × Treatment × Sex study of 40 samples with planted biology (glycolysis up, an ASS1 bottleneck, a glutaminase block, protein-only fatty-acid oxidation, mRNA-only cholesterol synthesis). **Download this example data** gives the tables, metadata, a README of comparisons to try and the ground truth.
- **Layers to integrate**: all three layers, or metabolomics with one enzyme layer.
- The step lists how each layer was processed (from its export), and what Integrix will still do.

## Step 2: Study design

- **Group column**: any metadata column (Condition, Treatment, Sex or a combined Group).
- **Case** and **Control (reference)**: the two groups compared. Results read case over control.
- **Batch column** and **Batch-correct** per layer: off by default, and off for layers the platform already corrected.

## Step 3: Statistics

- **Label permutations** (default 1,000) and **topology permutations** calibrate the reaction scores.
- **FDR (q)** threshold and **Pattern |z|** (the |z| used to call up / down / unchanged).
- Under **More options**: missing-value handling (k-nearest neighbours by default), the missing-value filter and the random seed.

## Step 4: Engines

- **Network, propagation and targets**: network propagation, choke points and target ranking.
- **Latent factors (MOFA+)**: shared and omics-specific variation.
- **Integrated biomarker signature**: multiblock sparse PLS-DA in nested cross-validation (needs at least 10 samples per group with all layers).

Press **Run analysis**. Finished analyses are kept on the server, so repeating a run with the same data and settings returns at once.

## Results

Four section cards sit above the results: **Results**, **Networks & targets**, **Advanced** and **Data & export**.

- **Overview**: metric cards, key findings, layer concordance across the top pathways and the reaction classification.
- **Reactions**: the pattern map (enzyme evidence E against metabolite balance M) and the reaction table with label, topology and balance q-values. Patterns include *Throughput-consistent*, *Candidate bottleneck*, *Constrained step*, *Non-abundance regulation* and *Compensated*.
- **Cascade**: enzyme → reaction → metabolite ribbons for a pathway, the top reactions or one enzyme.
- **Pathways**: IPA-style canonical pathways (Fisher's exact test, bars coloured by a directional z-score) and pathway activity per layer, classed as shared-concordant, shared-discordant or omics-specific.

## Networks & targets

- **Pathway maps & network**: IPA-style pathway maps (molecule shapes, measured red / green, predicted orange / blue, relationship-coloured edges, subcellular bands) and a network studio with 14 style presets (IPA classic, KEGG map, Cytoscape, journal minimal, dark, colour-blind safe, greyscale print, module map, hub focus, chord, hive plot, arc diagram and more). Exports: Cytoscape (.cyjs), GraphML, node and edge tables.
- **Targets**: enzymes ranked by robust rank aggregation of reaction coherence, mRNA–protein agreement, choke-point status and network propagation. Targets are hypotheses for experimental follow-up.

## Advanced

- **Concordance**: gene–protein agreement (coupled, opposing, protein-only, mRNA-only), when both transcriptomics and proteomics are integrated.
- **Factors**: MOFA+ variance explained per layer, factor scores, covariate associations and pathway enrichment of the loadings.
- **Signature**: cross-validated AUC with a bootstrap confidence interval, permutation p-value and the stability-selected panel. Validate any panel in an independent cohort before clinical use.

## Data & export

- **Data & mapping**: samples per layer, the evidence tier and identifier-mapping coverage, with the unmapped features to download.
- **Export**: the Excel workbook of every table, the run manifest (inputs, versions, settings, seeds) and an auto-written methods text that cites each platform's normalization. Every chart has **Download publication figure** under it. The file is the same chart as on screen (same data, marks and colours), typeset for a journal: a column width (single 3.5 in / 89 mm, 1.5 column 5.5 in / 140 mm or double 7.2 in / 183 mm; dense charts offer only the widths their labels fit), text at a real point size for the printed figure (8 pt by default; tick labels and legend 1 pt smaller, never below 6 pt), the height set by the content, tight margins and no title inside the image (it belongs in the figure legend). Formats: PNG, TIFF, JPEG, SVG, PDF, EPS or all six in a ZIP, at 300, 600 or 1200 DPI; SVG, PDF and EPS are vector files. **Build figure set (ZIP)** gives every chart at once in the same style, with the settings currently chosen on each page. The prioritised targets list, the cascade and the pathway maps are drawn for print and keep their own publication layout.
- **Methods**: the methods paragraph and references.

## Tips & troubleshooting

- **"The analysis could not run. Engine A needs a metabolomics layer"**: Integrix always needs a metabolite layer (metabolomics or lipidomics) and at least one of transcriptomics or proteomics.
- **Few observable reactions**: map more metabolites with HMDB, KEGG or ChEBI IDs, and keep a Gene column in proteomics.
- **Tier 2 or 3 instead of Tier 1**: the layers do not share all samples. Check that sample names match exactly; Data & mapping lists which samples each layer has.
- **Small studies**: with few samples per group, permutations cannot reach very small p-values; the app warns you. Treat results as exploratory.
- **Two platforms sent the metabolite layer**: tick only one of Metabolomics and Lipidomics under From this Omnix session.
- **Send to Integrix is greyed out or missing**: finish Normalization in the platform first. Sending your own data needs a subscription that includes Integrix (Multi-Omics, Lab or Enterprise).

## Glossary

- **Coherence (C)**: enzyme evidence × metabolite balance for one reaction; positive when they point the same way.
- **Enzyme evidence (E)**: protein z-score combined through the reaction's gene rules (complex: weakest subunit; isozymes: abundance-weighted mean).
- **Metabolite balance (M)**: (mean z of products − mean z of substrates) / √2.
- **Evidence tier**: 1 = paired samples, 2 = partly overlapping, 3 = unpaired.
- **Choke point**: a reaction that is the only producer or consumer of a metabolite.
- **Omnix export**: the zip a platform writes with Download for Integrix (matrix, metadata and a `*_omnix.json` description).
