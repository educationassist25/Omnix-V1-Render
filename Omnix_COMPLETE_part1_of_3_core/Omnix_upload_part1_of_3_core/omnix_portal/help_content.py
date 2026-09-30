"""
Structured help for each platform: a quick start, the input files (what each column means, with a
live preview of the platform's own example file and a download of it), and the detailed guide
(omnix_portal/tutorials/<platform>.md, shown step by step).

Example files are the demo datasets that ship with each platform, so what a user downloads here is
exactly what the platform loads as its demo.
"""

QUICK_START = {
    "metabolomics": [
        ("Setup", "Data Upload",
         "Choose single or multiple datasets and the **Analysis type** (untargeted or targeted metabolomics or "
         "lipidomics). Load a demo dataset or your files, then click **Load & Validate Data**."),
        ("Preprocessing", "Data Cleaning & Imputation, QC Validation",
         "Filter features with too many missing values, impute the rest, and check reproducibility on the QC "
         "injections (CV), then click **Confirm QC & Proceed**."),
        ("Normalization", "Normalization",
         "Targeted assays: ISTD normalization, then strict log2. Untargeted assays: log2, then Median-IQR "
         "normalization."),
        ("Statistical Analysis", "Statistics",
         "Pick the grouping variable and compare two groups (Welch's t-test or Wilcoxon) or three or more "
         "(ANOVA with post-hoc tests)."),
        ("Explore and interpret", "PCA, Volcano, Heatmap, Boxplot, Biomarkers, Pathways, Network",
         "Visualize the results, list candidate biomarkers, run MSEA / MGPA pathway analysis and the "
         "Correlation Rewiring Map. Every figure and table can be downloaded."),
    ],
    "proteomics": [
        ("Setup", "Data Upload",
         "Choose the **Quantification method** (Label-Free or TMT). Load a demo dataset or your files, then "
         "click **Load & Validate Data**."),
        ("Preprocessing", "Data Cleaning & Imputation, QC Validation",
         "Filter proteins with too many missing values, impute the rest, check QC reproducibility, then click "
         "**Confirm QC & Proceed** (this carries the imputed data forward)."),
        ("Normalization", "Normalization",
         "Apply the log2 transformation, then normalize: reference-channel for TMT, median centering (default) "
         "or IQR for label-free."),
        ("Statistical Analysis", "Statistics",
         "Pick the grouping variable and compare two groups or run ANOVA across three or more."),
        ("Explore and interpret", "PCA, Volcano, Heatmap, Boxplot, Biomarkers, GSEA, P-P Interaction",
         "Visualize the results, list candidate biomarkers, run enrichment analysis and build the "
         "protein–protein interaction network. Every figure and table can be downloaded."),
    ],
    "transcriptomics": [
        ("Setup", "Data Upload",
         "Choose the **Input data type** (Raw Counts or Pre-normalized logCPM). Load the demo dataset or your "
         "files, then click **Load Dataset**."),
        ("Preprocessing", "Data Cleaning, QC Validation",
         "Filter lowly expressed genes (minimum CPM) and review library sizes, sample correlation and "
         "outlier flags."),
        ("Normalization", "Normalization",
         "Raw counts: choose the design factor (and optional covariates such as Batch) and run DESeq2 "
         "normalization. logCPM: use the matrix as it is."),
        ("Statistical Analysis", "Statistics",
         "Compare two groups, run DESeq2 contrasts, or ANOVA across three or more groups."),
        ("Explore and interpret", "PCA, Volcano, Heatmap, Boxplot, Biomarkers, GSEA",
         "Visualize the results, list candidate biomarker genes and run enrichment analysis. Every figure and "
         "table can be downloaded."),
    ],
}

# Each input file: title, required?, one-line purpose, column table, example file(s) (relative to the
# platform folder), checklist.
FILES = {
    "metabolomics": [
        {
            "title": "Peak area matrix", "required": True,
            "purpose": "One row per metabolite or lipid, one column per injection (biological samples and QCs).",
            "columns": [
                ("First column (e.g. Metabolite)", "Yes", "Metabolite or lipid name. Names must be unique.", "Glyceric acid"),
                ("One column per sample", "Yes", "Peak areas. Headers must match the Sample column of the metadata. "
                 "Non-numeric cells are read as missing.", "26207.5"),
                ("Internal standard row", "Targeted only", "The dataset's internal standard (ISTD) as an ordinary "
                 "row. Start the name with ISTD_ so it is selected automatically, or name it in the ISTD column "
                 "of the row annotations.", "ISTD_D4-Alanine"),
            ],
            "examples": [("Untargeted example", "sample_peak_area_matrix_untargeted.csv"),
                         ("Targeted example (with ISTD row)", "sample_peak_area_matrix_targeted.csv")],
        },
        {
            "title": "Sample metadata", "required": True,
            "purpose": "One row per injection, describing its group, whether it is a QC, and its batch.",
            "columns": [
                ("Sample", "Yes", "Exactly the column headers of the peak area matrix. Every sample must be listed.", "Control_1"),
                ("Group", "Yes", "Experimental group.", "Control"),
                ("IsQC", "No", "True / 1 / yes / qc marks a QC injection. Missing column = all biological.", "False"),
                ("Batch", "No", "Batch label, used by batch-wise Median-IQR normalization. Missing = batch 1.", "1"),
                ("Any other column", "No", "Kept. Text columns and numeric columns with at most 15 values become "
                 "grouping variables in Statistics, PCA, Heatmap and Boxplot.", "Diagnosis, Gender"),
            ],
            "examples": [("Example", "sample_metadata_untargeted.csv")],
        },
        {
            "title": "Metabolite row annotations", "required": False,
            "purpose": "Optional information per metabolite, used by the Heatmap, Pathway Analysis and, for targeted "
                       "assays, to link each metabolite to the dataset's internal standard.",
            "columns": [
                ("Metabolite", "Yes", "Same names as the first column of the peak area matrix.", "Glyceric acid"),
                ("HMDB_ID", "No", "HMDB identifier; improves pathway matching (names are used otherwise).", "HMDB0000139"),
                ("Class", "No", "Metabolite or lipid class, shown as a Heatmap row track.", "Metabolite"),
                ("ISTD", "Targeted", "The internal standard used to normalize the metabolite: the name of the ISTD "
                 "row in the matrix, the same for every metabolite of the dataset. The ISTD row names itself.",
                 "ISTD_D4-Alanine"),
                ("Other columns", "No", "e.g. Method, Pathway; shown as Heatmap row tracks.", "Untargeted Metabolomics"),
            ],
            "examples": [("Example", "metabolite_row_annotations_untargeted.csv"),
                         ("Targeted example (with ISTD column)", "metabolite_row_annotations_targeted.csv")],
        },
    ],
    "proteomics": [
        {
            "title": "Protein intensity matrix", "required": True,
            "purpose": "One row per protein, one column per sample (including QC replicates and, for TMT, the reference channel).",
            "columns": [
                ("First column (e.g. Protein)", "Yes", "Protein identifier. Must be unique.", "ABCA1_HUMAN"),
                ("Gene (second column)", "Recommended", "Gene symbol; header Gene, Gene Name, Gene_Name, GeneSymbol, "
                 "Gene Symbol, Gene Names or Genes. Used by GSEA and P-P Interaction.", "ABCA1"),
                ("One column per sample", "Yes", "Raw intensities. Headers must match the Sample column of the "
                 "metadata. Blank cells and zeros count as not detected.", "43504.97"),
            ],
            "examples": [("Label-Free example", "sample_protein_intensity_matrix_labelfree.csv"),
                         ("TMT example (with Reference_Pool channel)", "sample_protein_intensity_matrix_tmt.csv")],
        },
        {
            "title": "Sample metadata", "required": True,
            "purpose": "One row per sample, describing its group, whether it is a QC, and its batch.",
            "columns": [
                ("Sample", "Yes", "Exactly the sample column headers of the matrix. Every sample must be listed.", "Control_1"),
                ("Group", "Yes", "Experimental group or condition.", "Control"),
                ("IsQC", "No", "True / 1 / yes / QC marks a QC replicate. Missing column = all biological.", "False"),
                ("Batch", "No", "Batch label, used by batch-based IQR normalization. Missing = batch 1.", "1"),
                ("Any other column", "No", "Kept. Text columns and numeric columns with at most 15 values become "
                 "grouping and coloring variables.", "Diagnosis, Gender"),
            ],
            "examples": [("Example", "sample_metadata_labelfree.csv")],
        },
        {
            "title": "Protein row annotations", "required": False,
            "purpose": "Optional information per protein, shown as Heatmap row tracks.",
            "columns": [
                ("Protein", "Yes", "Same identifiers as the first column of the matrix.", "ABCC8_HUMAN"),
                ("Other columns", "No", "e.g. Gene, Method, Pathway.", "Amino Acid Metabolism"),
            ],
            "examples": [("Example", "protein_row_annotations_richdemo.csv")],
        },
    ],
    "transcriptomics": [
        {
            "title": "Expression matrix", "required": True,
            "purpose": "One row per gene, one column per sample: raw counts, or pre-normalized logCPM.",
            "columns": [
                ("First column (e.g. Gene)", "Yes", "Gene identifier, used as written. Must be unique. HGNC symbols "
                 "work best for pathway analysis.", "INSR"),
                ("One column per sample", "Yes", "Raw counts (no negative values; fractional counts are rounded) or "
                 "logCPM. Headers must match the Sample column of the metadata.", "4767"),
            ],
            "examples": [("Raw counts example", "sample_counts_matrix_richdemo.csv")],
        },
        {
            "title": "Clinical metadata", "required": True,
            "purpose": "One row per sample, with its group and any clinical variables.",
            "columns": [
                ("Sample", "Yes", "Exactly the column headers of the expression matrix. Every sample must be listed.", "HealthyControl_1"),
                ("Group", "Yes", "Main experimental or clinical group.", "Healthy Control"),
                ("IsQC", "No", "True / 1 / yes / qc marks a technical QC replicate (used only on QC Validation).", "False"),
                ("Batch", "No", "Batch label; can be used as a covariate. Missing = batch 1.", "1"),
                ("Any other column", "No", "Clinical variables. Text columns and numeric columns with at most 15 "
                 "values can be used for grouping, coloring or as covariates.", "Diagnosis, Age, Gender"),
            ],
            "examples": [("Example", "sample_metadata_richdemo.csv")],
        },
        {
            "title": "Gene row annotations", "required": False,
            "purpose": "Optional information per gene: Heatmap row tracks, and extra columns in Biomarker Discovery.",
            "columns": [
                ("First column (e.g. Gene)", "Yes", "Same identifiers as the expression matrix.", "INSR"),
                ("Pathway and other columns", "No", "e.g. Pathway.", "Insulin Signaling"),
            ],
            "examples": [("Example", "gene_row_annotations_richdemo.csv")],
        },
    ],
}

CHECKLIST = [
    "Sample names are identical in the data file's column headers and the metadata's **Sample** column "
    "(spelling, spaces and capitalization).",
    "Every sample in the data file appears in the metadata.",
    "Feature names in the first column are unique.",
    "Save as CSV (UTF-8) or Excel (.xlsx); keep one table per file, with the header in the first row. In Excel, "
    "put the table on the first sheet: Omnix reads only the first sheet.",
    "Not sure about the layout? Start from the Excel example files on the **Data preparation** page.",
    "Mark QC injections with **IsQC** = True so they are used for QC and kept out of the statistics.",
]
