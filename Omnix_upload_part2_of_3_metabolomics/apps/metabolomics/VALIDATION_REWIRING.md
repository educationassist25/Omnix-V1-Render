# Correlation Rewiring Map — validation

This file records how the Rewiring Map (`metabolomics_modules/rewiring.py`, `metabolomics_modules/rewiring_reactions.py`,
`metabolomics_modules/rewiring_ui.py`, `rewiring_app.py`, the *Network Analysis* page in `app.py`, and the reaction
library in `metabolomics_modules/data/humangem_*`) was built and tested. Latest run: 2026-09-27.

## How to reproduce

```bash
pip install -r requirements.txt pytest tabulate
python -m pytest tests/ -q               # 65 unit, UI and full-pipeline tests, ~80 s
python validate_reaction_library.py      # reaction-library checks (gold standard, controls, IDs, coverage)
python validate_rewiring.py              # simulation study below, ~2-3 min
# rebuild the reaction library from Human-GEM (needs PyYAML + RDKit; the app itself does not):
python build_reaction_library.py [path/to/Human-GEM]      # downloads the model files if no path is given
# optional browser tests (need Node + Playwright):
streamlit run rewiring_app.py --server.headless true --server.port 8599 &
node tests/e2e_browser.js http://localhost:8599
streamlit run app.py --server.headless true --server.port 8598 &
node tests/e2e_browser.js http://localhost:8598 --main
# no-overlap sweep over every network setting, on any interactive reports (HTML downloads):
node tests/e2e_layout.js report1.html report2.html
```

## 1. Reaction library (Human-GEM v2.0.1 → metabolite reaction pairs)

**Source.** Human-GEM v2.0.1 (SysBioChalmers, 2026-03-26; CC BY 4.0): 12,877 reactions. Cite Robinson JL et al.,
*Sci. Signal.* 13, eaaz1482 (2020) and Luo J et al., *PNAS* 123:e2516511123 (2026).

**Derivation** (`build_reaction_library.py`): compartments collapsed; pure transport dropped (4,129 reactions);
always-currency species removed (H⁺, H₂O, Pi, PPi, O₂, CO₂, NH₃, CoA, thiamin-PP, lipoamide, ions, electron
carriers ...); NAD(P)(H), FAD and nucleotides removed only when their cofactor partner is on the other side (so
NAD synthetase keeps NAD⁺, adenylate cyclase keeps ATP → cAMP); generic species (pools, R-group formulas, proteins
with >150 C) excluded; substrates matched to products by conserved atoms (mean of exact and skeleton maximum-common-
substructure fractions, carbon-count tie-breaker; acyl-CoAs and nucleotide sugars compared without the carrier);
routes through enzyme-bound acyl intermediates contracted (PDH, OGDH, BCKDH); transketolase/transaldolase paired
by the carbon unit they move. Human-GEM HMDB/KEGG IDs were checked against MetaboAI Pro's crosswalk and HMDB's own
names: 19 HMDB accessions that Human-GEM gives to two metabolites were resolved by name (e.g. HMDB0000618 belongs
to ribulose-5-phosphate, not ribose-5-phosphate); missing IDs are filled from the crosswalk only when unambiguous
and name-confirmed. 390 peptides are flagged and never used as two-step intermediates, and a two-step link
is rejected when both steps come from the same reaction (co-products of one hydrolysis, e.g. glutamate and glycine
from a glutathione conjugate, are not a route).

**Result:** 4,972 unordered metabolite pairs from 5,354 reactions, 3,178 metabolites. Rebuilding from the raw
Human-GEM files reproduces the bundled parquet files exactly.

### Automated checks (`validate_reaction_library.py`)

Metabolite names are resolved exactly as user data is (MSEA ID standardization → HMDB names → Human-GEM).

### 1. Gold-standard conversions: 89/89 correct

| substrate                  | product                    |   max_steps |   steps | result   | via / reaction                                               |
|:---------------------------|:---------------------------|------------:|--------:|:---------|:-------------------------------------------------------------|
| Glucose                    | Glucose-6-phosphate        |           1 |       1 | ok       | ATP:D-glucose 6-phosphotransferase                           |
| Glucose-6-phosphate        | Fructose-6-phosphate       |           1 |       1 | ok       | D-glucose-6-phosphate aldose-ketose-isomerase                |
| Fructose-6-phosphate       | Fructose-1,6-bisphosphate  |           1 |       1 | ok       | UTP:D-fructose-6-phosphate 1-phosphotransferase              |
| Fructose-1,6-bisphosphate  | Dihydroxyacetone phosphate |           1 |       1 | ok       | fructose-bisphosphate aldolase (DHAP)                        |
| Fructose-1,6-bisphosphate  | Glyceraldehyde-3-phosphate |           1 |       1 | ok       | fructose-bisphosphate aldolase (DHAP)                        |
| Dihydroxyacetone phosphate | Glyceraldehyde-3-phosphate |           1 |       1 | ok       | D-glyceraldehyde-3-phosphate aldose-ketose-isomerase         |
| 3-Phosphoglycerate         | 2-Phosphoglycerate         |           1 |       1 | ok       | D-phosphoglycerate 2,3-phosphomutase                         |
| 2-Phosphoglycerate         | Phosphoenolpyruvate        |           1 |       1 | ok       | 2-phospho-D-glycerate hydro-lyase (phosphoenolpyruvate-formi |
| Phosphoenolpyruvate        | Pyruvate                   |           1 |       1 | ok       | GTP:pyruvate 2-O-phosphotransferase                          |
| Pyruvate                   | Lactate                    |           1 |       1 | ok       | (S)-Lactate:NAD+ oxidoreductase                              |
| Pyruvate                   | Alanine                    |           1 |       1 | ok       | (R)-3-Amino-2-methylpropanoate:pyruvate aminotransferase     |
| Pyruvate                   | Acetyl-CoA                 |           1 |       1 | ok       | pyruvate:thiamin diphosphate acetaldehydetransferase (decarb |
| Pyruvate                   | Oxaloacetate               |           1 |       1 | ok       | Pyruvate:carbon-dioxide ligase (ADP-forming)                 |
| Oxaloacetate               | Citrate                    |           1 |       1 | ok       | acetyl-CoA:oxaloacetate C-acetyltransferase (thioester-hydro |
| Acetyl-CoA                 | Citrate                    |           1 |       1 | ok       | acetyl-CoA:oxaloacetate C-acetyltransferase (thioester-hydro |
| Citrate                    | cis-Aconitate              |           1 |       1 | ok       | citrate hydro-lyase (cis-aconitate-forming)                  |
| cis-Aconitate              | Isocitrate                 |           1 |       1 | ok       | isocitrate hydro-lyase (cis-aconitate-forming)               |
| Citrate                    | Isocitrate                 |           2 |       1 | ok       | citrate hydroxymutase                                        |
| Isocitrate                 | alpha-Ketoglutarate        |           1 |       1 | ok       | Isocitrate:NADP+ oxidoreductase (decarboxylating)            |
| alpha-Ketoglutarate        | Succinyl-CoA               |           1 |       1 | ok       | oxoglutarate dehydrogenase (succinyl-transferring) (AKG) → o |
| Succinyl-CoA               | Succinate                  |           1 |       1 | ok       | succinyl-CoA:acetoacetate CoA-transferase                    |
| Succinate                  | Fumarate                   |           1 |       1 | ok       | succinate:quinone oxidoreductase                             |
| Fumarate                   | Malate                     |           1 |       1 | ok       | (S)-malate hydro-lyase (fumarate-forming)                    |
| Malate                     | Oxaloacetate               |           1 |       1 | ok       | (S)-malate:NAD+ oxidoreductase                               |
| Glutamine                  | Glutamate                  |           1 |       1 | ok       | L-Glutamate:ammonia ligase (ADP-forming)                     |
| Glutamate                  | alpha-Ketoglutarate        |           1 |       1 | ok       | Phenylalanine Transaminase, Mitochondrial                    |
| Glutamine                  | alpha-Ketoglutarate        |           2 |       2 | ok       | glutamate                                                    |
| Aspartate                  | Oxaloacetate               |           1 |       1 | ok       | L-Aspartate:2-oxoglutarate aminotransferase                  |
| Aspartate                  | Asparagine                 |           1 |       1 | ok       | L-asparagine amidohydrolase                                  |
| Glutamate                  | GABA                       |           1 |       1 | ok       | L-glutamate 1-carboxy-lyase (4-aminobutanoate-forming)       |
| Glutamate                  | N-Acetylglutamate          |           1 |       1 | ok       | acetyl-CoA:L-glutamate N-acetyltransferase                   |
| Proline                    | Glutamate                  |           2 |       2 | ok       | 1-pyrroline-5-carboxylate                                    |
| Phenylalanine              | Tyrosine                   |           1 |       1 | ok       | L-phenylalanine,tetrahydrobiopterin:oxygen oxidoreductase(4- |
| Tyrosine                   | L-DOPA                     |           1 |       1 | ok       | L-Tyrosine,tetrahydrobiopterin:oxygen oxidoreductase (3-hydr |
| L-DOPA                     | Dopamine                   |           1 |       1 | ok       | 3,4-Dihydroxy-L-phenylalanine carboxy-lyase                  |
| Histidine                  | Urocanic acid              |           1 |       1 | ok       | L-histidine ammonia-lyase (urocanate-forming)                |
| Histidine                  | Histamine                  |           1 |       1 | ok       | L-histidine carboxy-lyase (histamine-forming)                |
| Leucine                    | 4-Methyl-2-oxopentanoate   |           1 |       1 | ok       | L-Leucine:2-oxoglutarate aminotransferase                    |
| Valine                     | 3-Methyl-2-oxobutanoate    |           1 |       1 | ok       | L-Valine:2-oxoglutarate aminotransferase                     |
| Arginine                   | Ornithine                  |           1 |       1 | ok       | L-Arginine amidinohydrolase                                  |
| Ornithine                  | Citrulline                 |           1 |       1 | ok       | Carbamoyl-phosphate:L-ornithine carbamoyltransferase         |
| Citrulline                 | Argininosuccinate          |           1 |       1 | ok       | L-Citrulline:L-aspartate ligase (AMP-forming)                |
| Argininosuccinate          | Arginine                   |           1 |       1 | ok       | 2-(Nomega-L-arginino)succinate arginine-lyase (fumarate-form |
| Arginine                   | Citrulline                 |           1 |       1 | ok       | L-Arginine, NADPH:Oxygen Oxidoreductase (Nitric-Oxide-Formin |
| Ornithine                  | Putrescine                 |           1 |       1 | ok       | L-ornithine carboxy-lyase (putrescine-forming)               |
| Putrescine                 | Spermidine                 |           1 |       1 | ok       | S-adenosylmethioninamine:putrescine 3-aminopropyltransferase |
| Spermidine                 | Spermine                   |           1 |       1 | ok       | S-adenosylmethioninamine:spermidine 3-aminopropyltransferase |
| Guanidinoacetate           | Creatine                   |           1 |       1 | ok       | S-Adenosyl-L-methionine:guanidinoacetate N-methyltransferase |
| Glycine                    | Guanidinoacetate           |           1 |       1 | ok       | L-Arginine:glycine amidinotransferase                        |
| Aspartate                  | Asparagine                 |           1 |       1 | ok       | L-asparagine amidohydrolase                                  |
| Glutamine                  | Glutamate                  |           1 |       1 | ok       | L-Glutamate:ammonia ligase (ADP-forming)                     |
| Methionine                 | S-Adenosylmethionine       |           1 |       1 | ok       | ATP:L-methionine S-adenosyltransferase                       |
| S-Adenosylmethionine       | S-Adenosylhomocysteine     |           1 |       1 | ok       | Methyltransferase Coq3                                       |
| S-Adenosylhomocysteine     | Homocysteine               |           1 |       1 | ok       | S-Adenosyl-L-homocysteine hydrolase                          |
| Homocysteine               | Methionine                 |           1 |       1 | ok       | 5-methyltetrahydrofolate:L-homocysteine S-methyltransferase  |
| Homocysteine               | Cystathionine              |           1 |       1 | ok       | L-serine hydro-lyase (adding homocysteine                    |
| Cystathionine              | Cysteine                   |           1 |       1 | ok       | L-cystathionine cysteine-lyase (deaminating                  |
| Glutathione                | Glutathione disulfide      |           1 |       1 | ok       | Glutathione: 5-HPETE oxidoreductase                          |
| Choline                    | Phosphocholine             |           1 |       1 | ok       | ATP:choline phosphotransferase                               |
| Choline                    | Betaine                    |           2 |       2 | ok       | betaine aldehyde                                             |
| Betaine                    | Dimethylglycine            |           1 |       1 | ok       | Trimethylaminoacetate:L-homosysteine S-methyltransferase     |
| Hypotaurine                | Taurine                    |           1 |       1 | ok       | hypotaurine:NAD+ oxidoreductase                              |
| Tryptophan                 | N-Formylkynurenine         |           1 |       1 | ok       | L-tryptophan:oxygen 2,3-oxidoreductase (decyclizing)         |
| N-Formylkynurenine         | Kynurenine                 |           1 |       1 | ok       | arylformamidase (L-formylkynurenine)                         |
| Tryptophan                 | Kynurenine                 |           2 |       2 | ok       | L-formylkynurenine                                           |
| Kynurenine                 | Kynurenic acid             |           1 |       1 | ok       | Kynurenine-Oxoglutarate Transaminase                         |
| Kynurenine                 | 3-Hydroxykynurenine        |           1 |       1 | ok       | L-Kynurenine,NADPH:oxygen oxidoreductase (3-hydroxylating)   |
| Kynurenine                 | Anthranilic acid           |           1 |       1 | ok       | L-Kynurenine hydrolase                                       |
| 3-Hydroxykynurenine        | 3-Hydroxyanthranilic acid  |           1 |       1 | ok       | 3-Hydroxy-L-kynurenine hydrolase                             |
| 3-Hydroxykynurenine        | Xanthurenic acid           |           1 |       1 | ok       | Kynurenine-Oxoglutarate Transaminase                         |
| Tryptophan                 | 5-Hydroxytryptophan        |           1 |       1 | ok       | L-Tryptophan, Tetrahydrobiopterin:Oxygen Oxidoreductase (5-H |
| 5-Hydroxytryptophan        | Serotonin                  |           1 |       1 | ok       | 5-Hydroxy-L-tryptophan decarboxy-lyase                       |
| Tryptophan                 | Serotonin                  |           2 |       2 | ok       | 5-hydroxy-L-tryptophan                                       |
| Serotonin                  | 5-Hydroxyindoleacetic acid |           2 |       2 | ok       | 5-hydroxyindoleacetaldehyde                                  |
| Tryptophan                 | Tryptamine                 |           1 |       1 | ok       | L-tryptophan decarboxy-lyase                                 |
| Adenosine                  | Inosine                    |           1 |       1 | ok       | Adenosine aminohydrolase                                     |
| Inosine                    | Hypoxanthine               |           1 |       1 | ok       | inosine:phosphate alpha-D-ribosyltransferase                 |
| Hypoxanthine               | Xanthine                   |           1 |       1 | ok       | hypoxanthine:NAD+ oxidoreductase                             |
| Xanthine                   | Uric acid                  |           1 |       1 | ok       | xanthine:NAD+ oxidoreductase                                 |
| Guanosine                  | Guanine                    |           1 |       1 | ok       | guanosine:phosphate alpha-D-ribosyltransferase               |
| Guanine                    | Xanthine                   |           1 |       1 | ok       | Guanine aminohydrolase                                       |
| Uridine                    | Uracil                     |           1 |       1 | ok       | uridine:phosphate alpha-D-ribosyltransferase                 |
| Cytidine                   | Uridine                    |           1 |       1 | ok       | Cytidine aminohydrolase                                      |
| Carnitine                  | Acetylcarnitine            |           1 |       1 | ok       | Acetyl-CoA:carnitine O-acetyltransferase                     |
| Glycerol                   | Glycerol-3-phosphate       |           1 |       1 | ok       | ATP:glycerol 3-phosphotransferase                            |
| Glucose-6-phosphate        | 6-Phosphogluconolactone    |           1 |       1 | ok       | D-glucose-6-phosphate:NADP+ 1-oxidoreductase                 |
| Ribulose-5-phosphate       | Ribose-5-phosphate         |           1 |       1 | ok       | D-ribose-5-phosphate aldose-ketose-isomerase                 |
| Cholesterol                | 7alpha-Hydroxycholesterol  |           1 |       1 | ok       | cytochrome P450 (cholesterol)                                |
| Acetoacetate               | 3-Hydroxybutyrate          |           1 |       1 | ok       | (R)-3-Hydroxybutanoate:NAD+ oxidoreductase                   |

### 2. Negative controls: 24/24 correctly not linked

| a             | b             | scope    | result   |
|:--------------|:--------------|:---------|:---------|
| Aspartate     | Glutamate     | direct   | ok       |
| Alanine       | Glutamate     | direct   | ok       |
| Leucine       | Glutamate     | direct   | ok       |
| Pyruvate      | Glutamate     | direct   | ok       |
| Glucose       | Taurine       | direct   | ok       |
| Citrate       | Tryptophan    | direct   | ok       |
| Lactate       | Kynurenine    | direct   | ok       |
| Glutamine     | Serotonin     | direct   | ok       |
| Hypoxanthine  | Taurine       | direct   | ok       |
| Palmitic acid | Hypoxanthine  | direct   | ok       |
| Tryptophan    | Leucine       | direct   | ok       |
| Alanine       | Glycine       | direct   | ok       |
| Serine        | Leucine       | direct   | ok       |
| Glutamine     | Citrate       | direct   | ok       |
| Glutamine     | Asparagine    | direct   | ok       |
| Tyrosine      | Glutamate     | direct   | ok       |
| Phenylalanine | Glutamate     | direct   | ok       |
| Tryptophan    | Leucine       | ≤2 steps | ok       |
| Tryptophan    | Phenylalanine | ≤2 steps | ok       |
| Valine        | Histidine     | ≤2 steps | ok       |
| Glucose       | Taurine       | ≤2 steps | ok       |
| Hypoxanthine  | Taurine       | ≤2 steps | ok       |
| Glutamine     | Serotonin     | ≤2 steps | ok       |
| Lactate       | Kynurenine    | ≤2 steps | ok       |

### 3. Always-currency / generic species in pairs: 0

Context cofactors kept where they are real substrates/products (name: pairs, route degree; two-step intermediates need route degree ≤ 45):

AMP: 15, 15, ATP: 15, 15, NAD+: 8, 8, GMP: 7, 7, GTP: 6, 6, UMP: 4, 4, NADP+: 3, 3, FAD: 3, 3, FMN: 2, 2, dUTP: 2, 2, ADP: 1, 1, dUDP: 1, 1, UTP: 1, 1, dCTP: 1, 1, dTDP: 1, 1, IDP: 1, 1, dCDP: 1, 1, dGTP: 1, 1, CMP: 1, 1, dGDP: 1, 1, CTP: 1, 1, ITP: 1, 1

### 4. HMDB↔KEGG agreement with MetaboAI Pro's reference: 173 agree, 3 disagree, 13 reference HMDB IDs not in Human-GEM
- D-Fructose (HMDB0000660): reference KEGG C00095, Human-GEM C02336 [fructose]
- D-Galactose (HMDB0000143): reference KEGG C00124, Human-GEM C00984 [galactose]
- L-Carnitine (HMDB0000062): reference KEGG C00318, Human-GEM C15025 [L-carnitine]

### 5. Coverage of bundled datasets

| dataset            |   features |   ID-mapped |   in Human-GEM |   in a KEGG set |   direct reaction pairs |   2-step pairs |
|:-------------------|-----------:|------------:|---------------:|----------------:|------------------------:|---------------:|
| Targeted demo      |         60 |          59 |             45 |              50 |                       5 |             10 |
| Untargeted demo    |        300 |         300 |            223 |             248 |                     152 |            140 |
| Rich Clinical demo |        200 |         200 |            196 |             174 |                     181 |             89 |
| Rewiring demo      |         43 |          43 |             43 |              43 |                      27 |             23 |

**PASS**

The 3 HMDB↔KEGG "disagreements" are anomer/stereo variants of the same compound (β-D-fructose, α-D-galactose,
unspecified carnitine); matching uses HMDB first and name-confirmed KEGG aliases, so they resolve correctly.

### Hand audit of randomly sampled reactions

Automated gold standards were used while tuning the pairing rules, so random reactions were also reviewed by hand
(equation, derived pairs, judged against the textbook atom mapping):

| Round | Sample | Result | Systematic errors found → fix |
|---|---|---|---|
| 1 | 30 multi-substrate reactions | 27/30 correct | CoA-transferases paired CoA ester ↔ CoA ester; leftover species joined unrelated partners → compare acyl groups without CoA; leftovers need ≥ 50% shared atoms |
| 2 | 40 new reactions | issues in 7 | serine:pyruvate-type transaminations (ties), NAD⁺ removed where it is the product, UDP-sugar transferases → skeleton MCS term, context-dependent cofactors, nucleotide-sugar carrier stripping |
| 3 | 50 new reactions (25 multi-substrate, 25 single-side) | 48/50 correct (≈ 96% of pairs) | transketolase paired look-alike sugar phosphates → EC 2.2.1.1/2.2.1.2 carbon-unit rule. Remaining: γ-glutamyl transferase with a generic acceptor |

## 2. Simulation study (`validate_rewiring.py`)

Simulated panel: 43 metabolites, Control vs Tumor, planted rewiring (21 + 5 lost, 4 flipped, 5 gained pairs), two
balanced batches with a batch mean shift, and three outlier samples that create a decoy correlation. Each dataset
is a fresh random draw. Sections 1–4 use the annotation-column grouping (unchanged from the first version);
section 5 uses the KEGG library from the MSEA page, the reaction scopes and the pathway-level test.

## Rewiring Map validation (30 planted datasets, 30 null datasets, 30 small-n datasets per scenario; 129s)

### 1. Recovery of planted rewiring (n = 30 per group)

| Metric | Mean | Min |
|---|---|---|
| Recall, lost pairs | 0.84 | 0.35 |
| Recall, gained pairs | 0.85 | 0.20 |
| Recall, flipped pairs | 1.00 | 1.00 |
| Class accuracy of recovered pairs | 0.999 | 0.971 |
| False discovery proportion (target ≤ 0.05) | 0.015 | max 0.128 |

### 2. Global null, n = 30 per group

- Datasets with any rewired pair called: 0.00 (target ≤ 0.05)
- Global permutation p < 0.05: 0.07 (target ≈ 0.05); median p = 0.53

### 3. Small groups (n = 9 per group), share of null datasets with any call

| Scenario | Benjamini–Hochberg | Permutation FDR |
|---|---|---|
| Gaussian, 147 correlated features | 0.27 (mean 0.8 pairs) | 0.10 (mean 0.1 pairs) |
| Rich Clinical demo, labels shuffled | 0.50 (mean 0.9 pairs) | 0.10 (mean 0.1 pairs) |

### 4. Hypothesis-card checks

- Outlier decoy called as a pair in 27/30 datasets; card flagged (mixed/fragile) in 27/27
- Modules made only of true rewired pairs that were flagged non-robust: 0

### 5. Sample size: pair-level vs reaction-scoped vs pathway-level (20 planted + 20 null datasets per size)

Planted rewiring as in section 1 (no outlier decoy), KEGG library from the MSEA page. Reaction-linked recall is over the planted pairs that are linked by a Human-GEM reaction (direct or 2 steps).

| n per group | Label splits | Recall, all pairs | Recall of reaction-linked pairs: all-pairs test | …: reaction-scoped test | FDP | Tryptophan metabolism detected | Ala/Asp/Glu metabolism detected |
|---|---|---|---|---|---|---|---|
| 3 | 20 | 0.00 | 0.00 | 0.00 | 0.000 | 0.00 | 0.00 |
| 4 | 70 | 0.00 | 0.00 | 0.00 | 0.000 | 0.00 | 0.00 |
| 6 | > 200 (random) | 0.01 | 0.01 | 0.06 | 0.000 | 0.05 | 0.00 |
| 10 | > 200 (random) | 0.12 | 0.14 | 0.24 | 0.056 | 0.40 | 0.25 |
| 20 | > 200 (random) | 0.60 | 0.66 | 0.78 | 0.051 | 0.90 | 0.85 |

Under the null (no rewiring):

| n per group | Datasets with any pair called (all / reaction-scoped) | Datasets with any pathway q < 0.05 | Pathway p < 0.05 (target ≈ 0.05) |
|---|---|---|---|
| 3 | 0.00 / 0.00 | 0.00 | 0.000 |
| 4 | 0.10 / 0.05 | 0.00 | 0.033 |
| 6 | 0.05 / 0.00 | 0.00 | 0.043 |
| 10 | 0.05 / 0.10 | 0.05 | 0.048 |
| 20 | 0.05 / 0.05 | 0.00 | 0.072 |

**Reading the numbers**

* With 30 samples per group, pairs are recovered with the right class almost every time, and the share of false
  calls stays below the 5% target on average.
* **Reaction-scoped testing helps.** Restricting the test to reaction-linked pairs raises recall of the planted
  reaction-linked pairs from 0.14 to 0.24 at 10 per group and from 0.66 to 0.78 at 20 per group, with the null
  call rate unchanged.
* **The pathway-level test** detects the planted pathways in 85–90% of datasets at 20 per group and 25–40% at 10,
  and is calibrated under the null (share of pathway p < 0.05 between 0.03 and 0.07).
* **3–6 samples per group.** The engine runs and every test stays valid (no false pathway calls at 3 per group),
  but nothing is detectable: with 3 + 3 there are 20 label splits and, because the statistic is unchanged when the
  labels are swapped, the smallest possible p-value is 0.10; with 4 + 4 it is 0.029, and multiple testing then
  prevents any FDR < 0.05. The app says this and shows a ranked, clearly labelled exploratory list instead.
  Differential correlation needs roughly 10–20 samples per group for effects of this size.
* A regression found during this work: with 3 samples, correlations sit at ±0.99999 and Fisher z amplified the
  rounding difference between the observed statistic (raw data) and its own permutation (group-centred data),
  giving p = 0 and false pathway calls in 25% of null datasets. The observed statistic is now computed from the same
  matrix as the null; a test guards this.

## 3. Automated tests (`tests/`)

| File | What it covers |
|---|---|
| `test_rewiring.py` (25 tests) | Correlations match numpy/scipy (Pearson, Spearman); Fisher z, p and BH q match a hand computation and statsmodels; Welch t matches scipy; permutation q-values match a brute-force implementation of the definition and are monotone; classification rules; invariance to feature/sample order and to group mean shifts; planted-truth recovery and FDR across seeds; null calibration at n = 30 and n = 9; decoy, imputation-artifact and loading-artifact scenarios; missing values, zero variance, QC exclusion, feature cap, errors; pathway matching; strict-JSON payload; HTML escaping of hostile names. |
| `test_rewiring_reactions.py` (27 tests) | Reaction library passes the gold standard and negative controls; no currency hubs; provenance; reaction direction and 2-step routes (tryptophan → N-formylkynurenine → kynurenine); no amino-acid shortcuts through peptides; name-resolution layers and ID carry-over into MSEA libraries; KEGG anomer groups (sugar phosphates in glycolysis); multi-dataset prefixes; reaction and direct scopes restrict testing and BH; ratio shift equals a manual computation; pathway test finds the planted pathways and not the intact TCA cycle; pathway p-values equal a brute-force permutation; null pathway test calibrated; 3 per group runs with exact enumeration and correct p floor; 3 per group needs permutations; 2 per group refused; exact global p equals brute-force enumeration over all 70 splits at 4 + 4; small-group cards marked not assessable; null at 3 per group never calls a pathway (regression); co-products of one reaction are not a two-step route (regression); 27 names with special characters (Greek letters, Unicode dashes/primes/spaces/superscripts, adduct and retention-time tags, footnote marks, isomer groups, locant conventions) match the correct HMDB ID, with negative controls (3- vs 2-hydroxybutyrate stay distinct, lipid slashes are not split, unknowns stay unmapped) and annotation-column IDs taking precedence; network figure geometry never overlaps (metabolite names vs each other, every dot and the arcs; pathway labels outside the ring and apart; leader lines outside the ring), for the demo at four settings and a dense 150-metabolite ring with long names; p-value cutoff uses the permutation p and a cutoff of 1 calls every pair passing \|Δr\|; top 5/10/15/20 pathway selection shows exactly the selected pathways' metabolites; network and correlation figures render as PNG, JPEG, TIFF, SVG, PDF and EPS with 'Relative Abundance (…)' axes. |
| `test_rewiring_apps.py` (12 tests) | Streamlit AppTest: standalone app runs the demo end to end with all seven downloads and pathway-level summary; settings change prompts a rerun; small-group note on the Rich Clinical demo; upload instructions; reaction-linked scope from the UI; FDR q-value pair cutoff up to 1 (no cutoff-type choice), network default top 10 by p-value ≤ 0.05, network controls and figure downloads in all six formats from the UI; main app still loads; Network Analysis offers the demo without data; pipeline-style normalized log2 data in single- and multi-dataset mode with KEGG grouping; only the four pathway libraries offered; the page writes p-value and FDR q-value in full (options, summary, table and CSV columns incl. every pathway-table column, embedded view); other pages' state untouched. |
| `e2e_layout.js` + `overlap_check.js` | Opens interactive reports offline and runs the geometric no-overlap check for every "pathways shown" setting (top 5/10/15/20/all × p-value/FDR q-value × cutoff), with and without a selected (bold) metabolite. Passed on the demo (12–43 metabolites) and on a dense 154-metabolite report with long metabolite and pathway names. A negative control (a label placed on a dot and across edges) is flagged. |
| `test_rewiring_pipeline.py` (1 test) | Drives the real main-app UI with no hand-built state: Multiple datasets → Demo Data → Load → Auto-Process → Generate Combined Normalized Data → Network Analysis → Run. |
| `e2e_browser.js` | Real Streamlit server + Chromium, standalone app and `--main`: 43 nodes, rewired edges = table rows, reaction layer drawn, pathway table listed and row click highlights its metabolites, a reaction edge opens the reaction/ratio panel, top-5 / all / no-pathway network filters, default network is top 10 by p-value ≤ 0.05, geometric no-overlap check (every text box vs every other text, dot, edge, reaction edge, arc and leader line) at the default and with all pathways, full (untruncated) metabolite and wrapped pathway names with no text clipped in the network or the correlation plot, SVG/PNG/JPEG/TIFF/PDF buttons directly under the network and the correlation plot, and six in-view downloads from inside the Streamlit iframe checked as valid files, no bare 'q' in headers or labels, 'Relative Abundance' axes, the legend inside the network image (edge classes, reactions, log2 FC scale), no separate figure-format section, 6 cards including a flagged one, iframe auto-sizing, card highlight, row selection, class-chip filtering, edge tooltip, no JavaScript errors. |

Latest run: `pytest tests/` 65 passed; `e2e_layout.js` passed (80 settings); `validate_reaction_library.py` PASS; `validate_rewiring.py` as above;
browser test passed on both apps. The interactive report was also rendered at n = 30 and n = 3 without JavaScript
errors.

## 4. Known limitations

* **Reaction library.** Human-GEM is a generic human model: tissue-specific or microbial reactions are absent, and
  about a quarter of untargeted features (e.g. drug metabolites, many lipid species written as pools) have no match
  in the network. Pairing is an approximation of atom mapping (≈ 96% of pairs correct in the hand audit); group-
  transfer reactions with generic acceptors (γ-glutamyl transferase) and acyl/alkyl exchange in ether-lipid
  synthesis can be paired wrongly. Two-step links follow Human-GEM reaction directions and skip hubs (> 45
  partners) and peptides, so some real two-step routes through hubs are not shown.
* **KEGG library (also affects the MSEA page).** The bundled KEGG library uses KEGG's anomer-specific IDs, so
  plain D-glucose-6-phosphate, D-fructose-6-phosphate and fructose-1,6-bisphosphate are not in the glycolysis set
  there. The Rewiring Map treats anomers as one compound; the MSEA page does not yet.
* **Small groups.** Below about 6 samples per group nothing can reach significance (see section 2); results are an
  exploratory ranking. Robustness checks that resample samples are not assessable below 5 per group.
* **Pathway assignment for the ring** picks, per metabolite, the best-covered pathway with ≥ 3 measured members;
  metabolites in several pathways are listed with all of them in the ID/pathway mapping download.
* Interpretations on the cards are rule-based text keyed on pathway names. They suggest what to test; they are not
  findings.
* Correlation is computed on the matrix the user supplies. In the main app that is the normalized log2 data; in
  the standalone app, uploads are log2-transformed and (by default) median-centred per sample.
* The existing `test_pipeline.py` fails before reaching any new code (`normalization.log2_transform` now returns
  one value, the test unpacks two). This was already the case in the uploaded version and was left unchanged.
