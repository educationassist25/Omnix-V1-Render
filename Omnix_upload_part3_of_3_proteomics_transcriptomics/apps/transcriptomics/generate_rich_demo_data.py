"""
generate_rich_demo_data.py - Rich Clinical Demo dataset generator for BulkRNAAI Pro.

Produces three files (written next to this script):
  1. sample_counts_matrix_richdemo.csv   - 2,200 genes x 46 samples, RAW integer counts
  2. sample_metadata_richdemo.csv        - 46 samples x clinical/demographic variables
  3. gene_row_annotations_richdemo.csv   - Gene -> Pathway annotation (for Heatmap row tracks)

Design
------
- 4 Diagnosis groups (Healthy Control, Prediabetic, Type 2 Diabetes, Metabolic
  Syndrome) x 10 biological samples each = 40, + 6 technical QC replicates
  (pooled-sample resequencing) = 46 total. "Diagnosis" doubles as the "Group"
  column the app's PCA/Statistics/Heatmap/Boxplot grouping-variable logic
  reads by default; Age/Gender/Treatment/Ethnicity/Body Weight/Batch ride
  along as extra clinical covariates.
- 300 genes carry REAL HGNC gene symbols drawn from 15 curated pathways
  (insulin signaling, glycolysis, oxidative phosphorylation, adipogenesis,
  inflammatory response, TNF/NF-kB, interferon-gamma, adaptive immunity,
  apoptosis, cell cycle, fatty acid metabolism, complement, xenobiotic
  metabolism, mTORC1 signaling, hypoxia) -- the SAME 15 gene sets bundled in
  gene_sets_demo.gmt, so the built-in demo GSEA returns genuinely enriched
  pathways rather than noise. The remaining ~1,900 genes are realistic
  background (GENE00001-style placeholders) at a range of expression levels,
  so gene filtering/QC/normalization all have real low-expression genes to
  filter out.
- Counts are simulated from a negative-binomial model (the standard RNA-seq
  generative model DESeq2/edgeR themselves assume): a lognormal per-gene mean
  count, a mean-dependent dispersion (low-count genes are noisier -- the
  well-known RNA-seq mean-variance trend), a per-sample library-size factor,
  a small per-batch multiplicative effect, and disease-driven log2 fold
  changes injected on a biologically themed subset of genes per diagnosis
  group (meaningful for DESeq2/PCA/heatmap/volcano/GSEA). A continuous
  Age-linked effect is injected on a handful of aging-associated genes for a
  genuine (if modest) relationship to a continuous covariate. QC replicates
  are simulated around the cohort-wide mean profile with only technical
  (sequencing-depth) noise, no biological variation -- what real pooled
  technical replicates should look like.
"""
import os
import numpy as np
import pandas as pd

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
SEED = 2024

DIAGNOSES = ["Healthy Control", "Prediabetic", "Type 2 Diabetes", "Metabolic Syndrome"]
N_PER_DIAGNOSIS = 10
N_QC = 6
N_BACKGROUND_GENES = 1900

GENDERS = ["Female", "Male"]
ETHNICITIES = ["Caucasian", "Hispanic/Latino", "African American", "Asian", "Other"]
TREATMENTS_BY_DIAGNOSIS = {
    "Healthy Control": ["Untreated"],
    "Prediabetic": ["Untreated", "Diet & Exercise"],
    "Type 2 Diabetes": ["Metformin", "Insulin", "Diet & Exercise"],
    "Metabolic Syndrome": ["Diet & Exercise", "Metformin"],
}

# Curated pathway gene sets — identical to gene_sets_demo.gmt, kept inline so this
# script has no runtime dependency on parsing the .gmt file.
PATHWAY_GENES = {
    "Insulin Signaling": ["INSR", "IRS1", "IRS2", "PIK3CA", "PIK3R1", "AKT1", "AKT2", "PDPK1",
                           "FOXO1", "FOXO3", "GSK3B", "PPP1R3A", "SLC2A4", "SLC2A1", "PRKCZ",
                           "SOS1", "GRB2", "RAF1", "MAPK1", "MAPK3"],
    "Glycolysis": ["HK1", "HK2", "GPI", "PFKL", "PFKM", "ALDOA", "TPI1", "GAPDH", "PGK1",
                   "PGAM1", "ENO1", "PKM", "LDHA", "LDHB", "PDK1", "PDK4", "PCK1", "PCK2",
                   "G6PC", "FBP1"],
    "Oxidative Phosphorylation": ["NDUFA1", "NDUFA9", "NDUFB3", "NDUFS1", "SDHA", "SDHB",
                                   "UQCRC1", "UQCRC2", "CYC1", "COX4I1", "COX5A", "COX6A1",
                                   "ATP5F1A", "ATP5F1B", "ATP5MC1", "PPARGC1A", "NRF1", "TFAM",
                                   "CPT1A", "CPT2"],
    "Adipogenesis": ["PPARG", "CEBPA", "CEBPB", "FABP4", "ADIPOQ", "LEP", "LEPR", "PLIN1",
                      "LIPE", "FASN", "SCD", "SREBF1", "ACACA", "GLUT4", "UCP1", "PCK1",
                      "RETN", "NAMPT", "CD36", "LPL"],
    "Inflammatory Response": ["IL6", "IL1B", "TNF", "NFKB1", "NFKB2", "RELA", "NLRP3", "CCL2",
                               "CCL5", "CXCL8", "CXCL10", "TLR4", "TLR2", "MYD88", "IRAK1",
                               "PTGS2", "ICAM1", "VCAM1", "SELE", "CRP"],
    "TNFa Signaling via NFkB": ["TNFAIP3", "NFKBIA", "NFKB1", "RELB", "TRAF1", "BIRC3", "CXCL1",
                                 "CXCL2", "CXCL3", "CCL20", "IL1A", "IL1B", "SOD2", "IER3",
                                 "PTGS2", "JUNB", "FOS", "EGR1", "SOCS3", "BCL3"],
    "Interferon-gamma Response": ["STAT1", "IRF1", "IRF7", "IRF9", "JAK2", "IFNGR1", "IFNGR2",
                                   "GBP1", "GBP2", "CXCL9", "CXCL10", "CXCL11", "HLA-A", "HLA-B",
                                   "HLA-DRA", "B2M", "PSMB9", "TAP1", "ICAM1", "SOCS1"],
    "Adaptive Immune Response": ["CD3D", "CD3E", "CD4", "CD8A", "CD19", "CD28", "CD40", "CD40LG",
                                  "IL2", "IL2RA", "IL4", "IL10", "FOXP3", "IKZF1", "LCK", "ZAP70",
                                  "PTPRC", "MS4A1", "CD79A", "CD79B"],
    "Apoptosis": ["CASP3", "CASP7", "CASP8", "CASP9", "BAX", "BAK1", "BCL2", "BCL2L1", "BID",
                  "FAS", "FASLG", "TNFRSF10A", "TNFRSF10B", "CYCS", "APAF1", "TP53", "BBC3",
                  "PMAIP1", "XIAP", "DIABLO"],
    "Cell Cycle": ["CCND1", "CCNE1", "CCNB1", "CCNA2", "CDK1", "CDK2", "CDK4", "CDK6", "CDKN1A",
                   "CDKN1B", "CDKN2A", "RB1", "E2F1", "MKI67", "PCNA", "MCM2", "MCM5", "TOP2A",
                   "AURKA", "PLK1"],
    "Fatty Acid Metabolism": ["ACOX1", "ACADVL", "ACADM", "HADHA", "HADHB", "CPT1A", "CPT1B",
                               "CPT2", "SCD", "FASN", "ACACA", "ACACB", "ELOVL6", "PPARA",
                               "PPARD", "FABP1", "FABP3", "SLC27A1", "SLC27A2", "ACSL1"],
    "Complement": ["C1QA", "C1QB", "C1QC", "C3", "C3AR1", "C5AR1", "CFB", "CFD", "CFH", "CFI",
                   "SERPING1", "CD55", "CD46", "CD59", "MASP1", "MASP2", "MBL2", "ITGAM",
                   "ITGB2", "VSIG4"],
    "Xenobiotic Metabolism": ["CYP1A1", "CYP1A2", "CYP2E1", "CYP3A4", "CYP2C9", "GSTA1", "GSTM1",
                               "GSTP1", "UGT1A1", "UGT2B7", "NQO1", "ABCB1", "ABCC2", "SLC22A1",
                               "ALDH1A1", "ALDH3A1", "AOX1", "SULT1A1", "EPHX1", "NAT2"],
    "mTORC1 Signaling": ["MTOR", "RPTOR", "RPS6KB1", "EIF4EBP1", "RPS6", "AKT1", "TSC1", "TSC2",
                          "RHEB", "PRKAA1", "PRKAA2", "DDIT4", "SLC7A5", "SLC3A2", "HIF1A",
                          "VEGFA", "SREBF1", "INSR", "IGF1", "IGF1R"],
    "Hypoxia": ["HIF1A", "EPAS1", "VEGFA", "VEGFB", "SLC2A1", "PDK1", "LDHA", "CA9", "BNIP3",
                "NDRG1", "PGK1", "ENO1", "ADM", "ANGPTL4", "PFKFB3", "EGLN1", "EGLN3", "VHL",
                "ARNT", "HK2"],
}

# Which pathways are up/down-regulated in which diagnosis, relative to Healthy Control
# (biologically themed, not meant as a literal clinical claim): Prediabetic and Type 2
# Diabetes show impaired insulin signaling / OXPHOS and rising inflammation; Metabolic
# Syndrome adds a stronger adipogenesis/lipid and complement signature.
DIAGNOSIS_PATHWAY_EFFECTS = {
    # NOTE: Prediabetic vs Healthy Control is the FIRST comparison a new user
    # sees by default (both selectboxes default to the first two groups in
    # discovery order). A previous version of this generator gave Prediabetic
    # only very subtle effects (|log2FC| ~0.5-0.9), which -- combined with
    # n=10/group and BH correction across ~2,100 genes -- produced ZERO genes
    # at FDR<0.05 for that default comparison (verified empirically through
    # the real pyDESeq2 pipeline). That is a bad first impression for a demo
    # dataset even though it's individually "realistic" (prediabetes is a
    # subtle transcriptomic state) -- so effect sizes below are calibrated,
    # via the real pipeline, to land a meaningful (double-digit) number of
    # FDR<0.05 genes for EVERY disease-vs-Healthy-Control pair, scaled so
    # Prediabetic < Type 2 Diabetes < Metabolic Syndrome in overall magnitude
    # (a defensible clinical severity gradient), rather than only for the
    # strongest pair.
    "Prediabetic": {"Insulin Signaling": -1.6, "Oxidative Phosphorylation": -1.0,
                    "Inflammatory Response": 1.3, "Glycolysis": 0.9},
    "Type 2 Diabetes": {"Insulin Signaling": -2.6, "Oxidative Phosphorylation": -1.7,
                        "Inflammatory Response": 2.1, "TNFa Signaling via NFkB": 1.6,
                        "Glycolysis": 1.4, "Hypoxia": 1.0, "mTORC1 Signaling": 0.7},
    "Metabolic Syndrome": {"Insulin Signaling": -2.0, "Adipogenesis": 2.4,
                           "Fatty Acid Metabolism": -1.4, "Inflammatory Response": 2.4,
                           "TNFa Signaling via NFkB": 1.8, "Complement": 1.5,
                           "Interferon-gamma Response": 1.2},
}

AGING_GENES = ["CDKN2A", "TP53", "IL6", "NFKB1", "SOD2", "TERT", "SIRT1"]

# Fraction of raw-count cells (biological samples only) turned into genuinely
# missing (NaN) values, so the app's Data Cleaning & Imputation tab has real
# missing-value handling to demonstrate on the built-in demo dataset.
MISSING_FRACTION = 0.05


def _inject_missing_values(counts_df: pd.DataFrame, meta: pd.DataFrame, rng) -> pd.DataFrame:
    """
    Realistically drop a fraction of raw-count cells to NaN (dropout / failed
    quantification), biased toward lower-expression genes -- in real RNA-seq,
    missing/unreliable calls are far more common for weakly-expressed
    transcripts than highly-expressed ones. Pooled QC replicates are left
    untouched (a technical replicate is expected to be complete), so their
    CV/correlation QC metrics stay meaningful.
    """
    out = counts_df.astype(float).copy()
    qc_samples = set(meta.loc[meta["IsQC"], "Sample"].values)
    bio_col_positions = np.array([j for j, c in enumerate(out.columns) if c not in qc_samples])

    gene_mean = out.iloc[:, bio_col_positions].mean(axis=1).values
    rank = pd.Series(gene_mean).rank(method="average").values
    weight = (len(rank) - rank + 1.0)
    weight = weight / weight.sum()

    n_total_cells = len(gene_mean) * len(bio_col_positions)
    n_missing = int(round(MISSING_FRACTION * n_total_cells))

    row_idx = rng.choice(len(gene_mean), size=n_missing, p=weight)
    col_idx = rng.choice(bio_col_positions, size=n_missing)

    arr = out.values
    arr[row_idx, col_idx] = np.nan
    return pd.DataFrame(arr, index=out.index, columns=out.columns)


def _build_sample_metadata(rng):
    sample_names, diagnosis, is_qc, batch = [], [], [], []
    age, gender, treatment, ethnicity, body_weight = [], [], [], [], []

    for dx in DIAGNOSES:
        age_center = 38 if dx == "Healthy Control" else rng.uniform(48, 58)
        weight_center = 68 if dx == "Healthy Control" else rng.uniform(82, 95)
        for i in range(N_PER_DIAGNOSIS):
            sample_names.append(f"{dx.replace(' ', '')}_{i+1}")
            diagnosis.append(dx)
            is_qc.append(False)
            batch.append("1" if i < N_PER_DIAGNOSIS // 2 else "2")
            age.append(int(np.clip(rng.normal(age_center, 8), 20, 85)))
            gender.append(rng.choice(GENDERS))
            treatment.append(rng.choice(TREATMENTS_BY_DIAGNOSIS[dx]))
            ethnicity.append(rng.choice(ETHNICITIES, p=[0.45, 0.20, 0.15, 0.15, 0.05]))
            body_weight.append(round(float(np.clip(rng.normal(weight_center, 12), 45, 140)), 1))

    mean_age = int(np.mean(age))
    mean_weight = round(float(np.mean(body_weight)), 1)
    for i in range(N_QC):
        sample_names.append(f"QC_{i+1}")
        diagnosis.append("QC")
        is_qc.append(True)
        batch.append("1" if i < N_QC // 2 else "2")
        age.append(mean_age)
        gender.append("Pooled")
        treatment.append("N/A (pooled QC)")
        ethnicity.append("Pooled")
        body_weight.append(mean_weight)

    meta = pd.DataFrame({
        "Sample": sample_names, "Group": diagnosis, "IsQC": is_qc, "Batch": batch,
        "Diagnosis": diagnosis, "Age": age, "Gender": gender, "Treatment": treatment,
        "Ethnicity": ethnicity, "Body Weight": body_weight,
    })
    return meta


def _build_gene_universe(rng):
    """Returns (gene_names, gene_pathway_map, base_mean_log). Curated pathway genes
    get realistic mid/high expression; background genes span a wide range including
    plenty of low-expression genes (for gene filtering to have something to do)."""
    pathway_genes, gene_pathway_map = [], {}
    for pathway, genes in PATHWAY_GENES.items():
        for g in genes:
            if g not in gene_pathway_map:  # a few genes legitimately appear in >1 pathway; keep first
                gene_pathway_map[g] = pathway
                pathway_genes.append(g)

    background_genes = [f"GENE{str(i+1).zfill(5)}" for i in range(N_BACKGROUND_GENES)]
    all_genes = pathway_genes + background_genes

    # log-mean baseline count per gene: pathway genes centered higher (they're
    # curated "real" genes, biologically expressed in this tissue context);
    # background spans a broad realistic range including a long low-expression tail.
    base_mean_log = np.concatenate([
        rng.normal(6.0, 1.0, size=len(pathway_genes)),      # ~e^6 ≈ 400 base counts
        rng.normal(3.5, 2.2, size=len(background_genes)),   # wide range incl. near-zero
    ])
    return all_genes, gene_pathway_map, base_mean_log


def generate():
    rng = np.random.default_rng(SEED)

    meta = _build_sample_metadata(rng)
    n_samples = len(meta)
    genes, gene_pathway_map, base_mean_log = _build_gene_universe(rng)
    n_genes = len(genes)
    gene_index = {g: i for i, g in enumerate(genes)}

    # --- per-sample library size (sequencing depth), realistic 12-35M reads, with a
    # mild batch effect on top -----------------------------------------------------
    lib_size = rng.uniform(1.2e7, 3.5e7, size=n_samples)
    batch_factor = np.where(meta["Batch"].values == "1", 1.0, rng.uniform(0.93, 1.07))
    lib_size = lib_size * batch_factor
    # Scale factor applied per-gene below: base_mean_log is calibrated so that
    # base_mean * this factor lands total simulated counts per sample in a
    # realistic ~10-30 million read range for a ~2,200-gene demo panel.
    lib_size_norm = lib_size / 1e6

    # --- per-gene log2 fold-change matrix (n_genes x n_samples), diagnosis- and
    # age-driven, built once and applied multiplicatively to the base mean ---------
    log2fc = np.zeros((n_genes, n_samples))
    sample_dx = meta["Diagnosis"].values
    sample_age = meta["Age"].values.astype(float)
    mean_age_all = sample_age[meta["IsQC"].values == False].mean()

    for dx, effects in DIAGNOSIS_PATHWAY_EFFECTS.items():
        dx_mask = sample_dx == dx
        n_dx = dx_mask.sum()
        for pathway, mean_effect in effects.items():
            pw_genes = [g for g in PATHWAY_GENES[pathway] if g in gene_index]
            pw_idx = [gene_index[g] for g in pw_genes]
            # per-gene effect scatter around the pathway's mean effect (not every
            # gene in a pathway moves identically -- realistic biological noise)
            per_gene_effect = rng.normal(mean_effect, 0.15, size=len(pw_idx))
            per_sample_noise = rng.normal(0, 0.15, size=(len(pw_idx), n_dx))
            log2fc[np.ix_(pw_idx, dx_mask)] += per_gene_effect[:, None] + per_sample_noise

    # continuous Age effect on a small handful of aging-associated genes
    for g in AGING_GENES:
        if g in gene_index:
            gi = gene_index[g]
            age_z = (sample_age - mean_age_all) / 15.0
            log2fc[gi, :] += 0.35 * age_z

    # QC replicates get NO biological (diagnosis/age) fold change -- pooled sample.
    qc_mask = meta["IsQC"].values
    log2fc[:, qc_mask] = 0.0

    # --- simulate negative-binomial counts -----------------------------------------
    # mu_ij = base_mean * 2^log2fc_ij * library_size_norm_j
    # dispersion: RNA-seq's well-known mean-variance trend -- low-count genes are
    # proportionally noisier. size (NB "r" parameter) shrinks the variance as
    # base expression rises, roughly matching real bulk RNA-seq behavior.
    base_mean = np.exp(base_mean_log)[:, None]  # genes x 1
    mu = base_mean * (2.0 ** log2fc) * lib_size_norm[None, :]
    mu = np.clip(mu, 1e-3, None)

    dispersion = 0.15 + 4.0 / (1.0 + base_mean[:, 0] / 50.0)  # higher for low-expression genes
    nb_size = np.clip(1.0 / dispersion, 0.3, 200.0)[:, None]  # genes x 1, broadcast over samples
    nb_size_b = np.broadcast_to(nb_size, mu.shape)
    prob = nb_size_b / (nb_size_b + mu)
    counts = rng.negative_binomial(nb_size_b, prob)

    counts_df = pd.DataFrame(counts, index=genes, columns=meta["Sample"].values)
    counts_df.index.name = "Gene"
    counts_df = _inject_missing_values(counts_df, meta, rng)

    # --- gene row annotations (Pathway) --------------------------------------------
    row_annot = pd.DataFrame({
        "Gene": genes,
        "Pathway": [gene_pathway_map.get(g, "Background / unannotated") for g in genes],
    })

    return counts_df, meta, row_annot


def main():
    counts_df, meta, row_annot = generate()
    counts_path = os.path.join(BASE_DIR, "sample_counts_matrix_richdemo.csv")
    meta_path = os.path.join(BASE_DIR, "sample_metadata_richdemo.csv")
    annot_path = os.path.join(BASE_DIR, "gene_row_annotations_richdemo.csv")

    # float_format keeps non-missing cells looking like plain integer counts
    # ("1234" not "1234.0"); NaN cells are still written as empty, regardless.
    counts_df.to_csv(counts_path, float_format="%.0f")
    meta.to_csv(meta_path, index=False)
    row_annot.to_csv(annot_path, index=False)

    n_missing = int(counts_df.isna().values.sum())
    n_cells = counts_df.size
    print(f"Wrote {counts_df.shape[0]} genes x {counts_df.shape[1]} samples -> {counts_path}")
    print(f"  ({n_missing} missing cells injected, {100 * n_missing / n_cells:.2f}% of {n_cells} total)")
    print(f"Wrote {meta.shape[0]} samples x {meta.shape[1]} fields -> {meta_path}")
    print(f"Wrote {row_annot.shape[0]} gene annotations -> {annot_path}")


if __name__ == "__main__":
    main()
