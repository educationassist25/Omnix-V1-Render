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
         "Choose single or multiple datasets and the **Analysis type** (Untargeted or Targeted Metabolomics). "
         "Load a demo dataset or your files, then click **Load & Validate Data**."),
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
    "fluxomics": [
        ("Setup", "Upload & inspect, Experimental design",
         "Choose the **Input format** (Compound Discoverer SIL exports or a generic isotopologue table + metadata). "
         "Load a demo dataset or your files, inspect the isotopologue channels, then check groups, time points, "
         "tracer and natural-abundance correction and click **Confirm design**."),
        ("Normalization", "Normalization",
         "log2, then robust IQR normalization of each labeled sample. The sidebar switch **Intensity dataset for "
         "downstream steps** chooses log2 + IQR (default) or log2 for every later intensity analysis."),
        ("Isotope Labeling", "Isotope incorporation, Label arrival",
         "Isotopologue fractions, incorporation, fractional enrichment and mean label number per sample and "
         "group × time, and the side-by-side label-arrival figure."),
        ("Statistical Analysis", "Statistics",
         "Two-group tests or ANOVA (including two-way group × time) on every isotopologue channel, with p value "
         "and FDR for each channel."),
        ("Explore and interpret", "Dynamic bar plots, Heatmaps, Pathway mapping",
         "Bar plots per metabolite, Z-score heatmaps and 25 curated pathway maps. Every figure and table can be "
         "downloaded where it is shown."),
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
    "fluxomics": [
        {
            "title": "Isotopologue table (generic format)", "required": True,
            "purpose": "One row per metabolite isotopologue (M+0, M+1 … M+n), one column per sample. Isotopologues "
                       "are never summed. Alternatively, upload Compound Discoverer SIL exports (below).",
            "columns": [
                ("Metabolite", "Yes", "Metabolite name; one block of rows per metabolite.", "D-Glucose"),
                ("Formula", "Recommended", "Molecular formula; sets n, the number of carbons (maximum M+n).", "C6H12O6"),
                ("Method", "No", "LC-MS method.", "HILIC_Neg"),
                ("KEGG, HMDB, PubChemCID, ChEBI, BiGG", "No", "Identifiers; KEGG IDs improve pathway mapping.", "C00031"),
                ("Isotopologue", "Yes", "M+0, M+1 … (or 0, 1 …). Each once per metabolite; a missing row = not "
                 "detected.", "M+0"),
                ("One column per sample", "Yes", "Peak area of that isotopologue. Headers must match the Sample column "
                 "of the metadata. Empty cell = missing value (not imputed).", "10839055"),
            ],
            "examples": [("Synthetic example (with 12C reference and blank)",
                          "example_data/generic/synthetic_isotopologue_table.csv"),
                         ("Template without 12C reference",
                          "example_data/templates/generic/Generic_Isotopologue_Table_without_12C_Ref.csv")],
        },
        {
            "title": "Sample metadata", "required": False,
            "purpose": "One row per sample column of the table: its type, group, time point and replicate. Optional "
                       "but recommended (without it, Step 2 starts from the sample names).",
            "columns": [
                ("Sample", "Yes", "Exactly the sample column names of the table.", "Ctrl_2h_R1"),
                ("Type", "Yes", "Labeled, Unlabeled reference, Blank or QC.", "Labeled"),
                ("Group", "Yes", "Experimental group.", "shControl"),
                ("Treatment", "No", "Drug, dose …; each group × treatment becomes its own group.", ""),
                ("Time / TimeUnit", "Yes (labeled)", "Labeling time as a number, and min, h or d. Empty for the "
                 "reference and blank.", "2 / h"),
                ("Replicate", "Recommended", "Biological replicate number (pairs samples for paired tests).", "1"),
                ("Include", "No", "TRUE / FALSE; FALSE leaves the sample out of every analysis.", "True"),
            ],
            "examples": [("Synthetic example", "example_data/generic/synthetic_sample_metadata.csv")],
        },
        {
            "title": "Compound Discoverer SIL export (alternative input)", "required": False,
            "purpose": "Instead of the generic table: one Compound Discoverer SIL export per LC-MS method (Excel or "
                       "CSV, data on the first sheet; the workbook name contains the method, e.g. …_RP_Pos.xlsx). "
                       "Groups, time points and sample numbers are read from the raw-file names.",
            "columns": [
                ("Name, Formula, Molecular Weight", "Yes", "Compound as annotated in Compound Discoverer.", "Creatine"),
                ("RT in min, Area, Status, Study File ID", "Yes", "Peak retention time, total area, CD status and file "
                 "ID; used to match and rank peaks.", "1.304"),
                ("Filename", "Yes", "<date>_<method>_<CellLine>_<Group>_<SampleNo>_<Time>.raw; time 12C = unlabeled "
                 "reference; a name with Blank = blank.", "08302026_RP_Pos_CellLine_GeneKD_1001_12Hr.raw"),
                ("Relative Area 0 … Relative Area 25", "Yes", "Area of each isotopologue M+0 … M+25 (unused = 0).",
                 "4037138058"),
            ],
            "examples": [],
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

# Platforms whose metadata works differently get their own checklist.
CHECKLIST_BY = {
    "fluxomics": [
        "Sample names are identical in the table's column headers and the metadata's **Sample** column "
        "(spelling, spaces and capitalization).",
        "Each isotopologue (M+0, M+1 …) appears **once** per metabolite and sample; never sum isotopologues.",
        "Give each metabolite its **Formula** (sets the maximum M+n) and, where possible, its **KEGG** ID (pathway maps).",
        "Mark the unlabeled reference and blanks with **Type** = Unlabeled reference / Blank, and give every labeled "
        "sample its **Time** and **Replicate**.",
        "Leave cells empty for compounds not reported in a sample: empty is missing, while values below the minimum "
        "value (default 1,000) in a labeled sample remove that isotopologue.",
        "Compound Discoverer: one export per LC-MS method, data on the first sheet, method in the workbook name.",
        "Save as CSV (UTF-8) or Excel (.xlsx), header in the first row, table on the first sheet.",
        "Not sure about the layout? Start from the templates on the **Data preparation** page.",
    ],
}


def checklist(platform):
    """The 'Before you upload' checklist for one platform."""
    return CHECKLIST_BY.get(platform, CHECKLIST)


# Lipidomics runs the Metabolomics pipeline with only the two lipidomics analysis types, so its quick start
# and input files follow Metabolomics, with lipid wording and the lipidomics demo files.
QUICK_START["lipidomics"] = [
    ("Setup", "Data Upload",
     "Choose single or multiple datasets and the **Analysis type** (Untargeted or Targeted Lipidomics). Load a "
     "demo dataset or your files, then click **Load & Validate Data**."),
    *QUICK_START["metabolomics"][1:4],
    ("Explore and interpret", "PCA, Volcano, Heatmap, Boxplot, Biomarkers, Pathways, Network",
     "Visualize the results, list candidate lipid biomarkers, run MSEA / MGPA pathway analysis (e.g. with the "
     "LIPID MAPS library) and the Correlation Rewiring Map. Every figure and table can be downloaded."),
]

FILES["lipidomics"] = [
    {
        "title": "Peak area matrix", "required": True,
        "purpose": "One row per lipid, one column per injection (biological samples and QCs).",
        "columns": [
            ("First column (e.g. Metabolite)", "Yes", "Lipid name. Names must be unique.", "PE(16:0/20:3(8Z,11Z,14Z))"),
            ("One column per sample", "Yes", "Peak areas. Headers must match the Sample column of the metadata. "
             "Non-numeric cells are read as missing.", "13374.43"),
            ("Internal standard row", "Targeted only", "The dataset's lipid internal standard (ISTD) as an ordinary "
             "row. Start the name with ISTD_ so it is selected automatically, or name it in the ISTD column "
             "of the row annotations.", "ISTD_PC(15:0/18:1-d7)"),
        ],
        "examples": [("Untargeted example", "sample_peak_area_matrix_multimethod_untargeted_lipidomics.csv"),
                     ("Targeted example (with ISTD row)", "sample_peak_area_matrix_multimethod_targeted_lipidomics.csv")],
    },
    FILES["metabolomics"][1],
    {
        "title": "Lipid row annotations", "required": False,
        "purpose": "Optional information per lipid, used by the Heatmap, Pathway Analysis and, for targeted "
                   "assays, to link each lipid to the dataset's internal standard.",
        "columns": [
            ("Metabolite", "Yes", "Same names as the first column of the peak area matrix.", "PE(16:0/20:3(8Z,11Z,14Z))"),
            ("HMDB_ID", "No", "HMDB identifier; improves pathway matching (names are used otherwise).", "HMDB0008936"),
            ("Class", "No", "Lipid class, shown as a Heatmap row track.", "Lipid"),
            ("ISTD", "Targeted", "The internal standard used to normalize the lipid: the name of the ISTD row in "
             "the matrix, the same for every lipid of the dataset. The ISTD row names itself.",
             "ISTD_PC(15:0/18:1-d7)"),
            ("Other columns", "No", "e.g. Method, Pathway; shown as Heatmap row tracks.", "Untargeted Lipidomics"),
        ],
        "examples": [("Example (lipidomics and metabolomics methods)", "metabolite_row_annotations_multimethod.csv")],
    },
]
