"""
Omnix portal settings: edit the values in this file to fill in your own details.

Everything in [SQUARE BRACKETS] is a placeholder that is shown on the site until you replace it.
"""

BRAND = "Omnix"
TAGLINE = "Integrated Multi-Omics Data Analytics"

# ---- Contact (Contact menu and Contact page) ------------------------------------------
CONTACT = {
    "email": "education.assist25@gmail.com",
    "phone": "",                                   # leave empty to hide
    "organization": "",                            # leave empty to hide
    "address": "3606 Park Vista Dr, Missouri City, TX 77459",
}

# ---- Subscription page -----------------------------------------------------------------
# Prices are list prices in US dollars per user; edit them freely. "platforms" is how many analysis
# platforms the plan unlocks ("one", "two" or "all"); "keys" is the number of access keys issued (users).
# To issue keys for a plan: python tools/make_access_key.py --plan "<name>" --platforms ... --count ...
SUBSCRIPTION_INTRO = ("Every visitor can run all Omnix analyses on the built-in demo datasets, free. A subscription "
                      "unlocks the same workflows for your own data. Choose one platform, two, or all five (with "
                      "Integrix tri-omics integration), for yourself or for your whole lab, billed monthly or "
                      "annually.")
CURRENCY = "$"
SUBSCRIPTION_PLANS = [
    {
        "id": "single", "name": "Single-Omics", "audience": "One researcher, one data type",
        "monthly": 25, "annual": 250, "platforms": "one", "keys": 1,
        "features": [
            "Your own data in one platform of your choice",
            "Metabolomics, Lipidomics, Proteomics, Transcriptomics or Fluxomics",
            "All analyses, figures and tables in that platform",
            "Demo access to the other four platforms and Integrix",
            "Email support",
        ],
    },
    {
        "id": "duo", "name": "Duo-Omics", "audience": "One researcher, two data types",
        "monthly": 45, "annual": 450, "platforms": "two", "keys": 1,
        "features": [
            "Your own data in two platforms of your choice",
            "e.g. Metabolomics + Lipidomics, or Proteomics + Transcriptomics",
            "All analyses, figures and tables in both platforms",
            "Demo access to the other three platforms and Integrix",
            "Email support",
        ],
    },
    {
        "id": "multi", "name": "Multi-Omics", "audience": "One researcher, all five data types",
        "monthly": 89, "annual": 890, "platforms": "all", "keys": 1, "badge": "Most popular",
        "features": [
            "Your own data in all five platforms",
            "Metabolomics, Lipidomics, Proteomics, Transcriptomics and Fluxomics",
            "Integrix: integrate transcriptomics, proteomics and metabolomics",
            "All analyses, figures and tables",
            "Priority email support",
        ],
    },
    {
        "id": "lab", "name": "Lab", "audience": "Research groups and core facilities",
        "monthly": None, "annual": 3900, "platforms": "all", "keys": 10,
        "features": [
            "Up to 10 users, each with a personal access key",
            "All five platforms + Integrix tri-omics integration",
            "One invoice for the group",
            "Onboarding session for your team",
            "Priority support",
        ],
    },
    {
        "id": "enterprise", "name": "Enterprise", "audience": "Institutions and industry",
        "monthly": None, "annual": None, "platforms": "all", "keys": None,
        "features": [
            "Any number of users",
            "All five platforms + Integrix tri-omics integration",
            "Purchase order and invoice billing",
            "Training and dedicated support",
            "Custom terms",
        ],
    },
]
ANNUAL_SAVING = "2 months free"
PRICING_NOTES = [
    "Prices are per year (or per month) and exclude applicable taxes.",
    "Academic and non-profit pricing is shown. Commercial organizations: please contact us for a quote.",
    "Annual plans are billed once and include two months free compared with monthly billing.",
    "Need more than 10 users? Add users to the Lab plan at $390 per user per year, or choose Enterprise.",
    "Single-Omics, Duo-Omics and Multi-Omics are personal plans for one named user. Research groups and core "
    "facilities sharing Omnix should choose the Lab plan, which gives every member their own key.",
]
PRICING_FAQ = [
    ("What can I do without a subscription?",
     "Everything except uploading your own files: every analysis in all five platforms runs on the built-in demo "
     "datasets, and all figures and tables from the demo can be downloaded."),
    ("What is Integrix, and which plans include it?",
     "Integrix integrates transcriptomics, proteomics and metabolomics (or lipidomics) of the same samples along "
     "gene → protein → enzyme → metabolite, with pathway maps, targets and an integrated signature. Normalize each "
     "dataset in its own platform, press Send to Integrix, and integrate them in one place. Integrix with your own "
     "data is included in Multi-Omics, Lab and Enterprise; every visitor can run it on its demo datasets."),
    ("How do I get access after subscribing?",
     "You receive a personal access key (OMX-XXXX-XXXX-XXXX-XXXX). Enter it below under Activate your subscription; "
     "uploads unlock immediately in the platforms your plan covers."),
    ("Can I change platforms or upgrade later?",
     "Yes. Single-Omics can be upgraded to Duo-Omics, Multi-Omics or Lab at any time; you pay only the difference for the rest "
     "of your billing period."),
    ("Can a key be shared?",
     "Each key is for one person. The Lab plan gives each member of your group their own key."),
    ("What happens when a subscription ends?",
     "Your key stops unlocking uploads after its end date, and Omnix returns to demo access. Renewing reactivates "
     "the same key."),
]

# ---- The five analysis platforms --------------------------------------------------------
# key -> (folder under apps/, menu name, one-line description, chips, highlights).
# "marker": the platform's own code folder inside apps/<folder>/ (default "<folder>_modules"); the
# portal uses it to recognise the platform's app.py in the repository.
# "app_family": Metabolomics and Lipidomics run the same apps/metabolomics/app.py; app_family limits each tab to
# its own two analysis types (Untargeted / Targeted), exactly as metabolomics_app.py and lipidomics_app.py do.
PLATFORMS = {
    "metabolomics": {
        "folder": "metabolomics",
        "app_family": "Metabolomics",
        "name": "Metabolomics",
        "subtitle": "LC-MS metabolomics",
        "chips": ["Untargeted", "Targeted"],
        "highlights": [
            "QC, imputation and normalization",
            "Statistics, PCA, volcano, heatmap, biomarkers",
            "Pathway analysis (MSEA, MGPA)",
            "Correlation rewiring network",
        ],
    },
    "lipidomics": {
        "folder": "metabolomics",
        "marker": "metabolomics_modules",
        "app_family": "Lipidomics",
        "name": "Lipidomics",
        "subtitle": "LC-MS lipidomics",
        "chips": ["Untargeted", "Targeted"],
        "highlights": [
            "QC, imputation and normalization",
            "Statistics, PCA, volcano, heatmap, biomarkers",
            "Pathway analysis (MSEA, MGPA)",
            "Correlation rewiring network",
        ],
    },
    "proteomics": {
        "folder": "proteomics",
        "name": "Proteomics",
        "subtitle": "LC-MS/MS protein quantification",
        "chips": ["Label-Free", "TMT"],
        "highlights": [
            "QC, imputation and normalization",
            "Statistics, PCA, volcano, heatmap, biomarkers",
            "Gene set enrichment (GSEA)",
            "Protein–protein interaction network",
        ],
    },
    "transcriptomics": {
        "folder": "transcriptomics",
        "name": "Transcriptomics",
        "subtitle": "Bulk RNA-seq differential expression",
        "chips": ["Raw counts", "logCPM"],
        "highlights": [
            "Gene filtering, QC, and outlier checks",
            "Normalization and contrasts",
            "PCA, volcano, heatmap, biomarkers",
            "Gene set enrichment (GSEA)",
        ],
    },
    "fluxomics": {
        "folder": "fluxomics",
        "marker": "core",
        "name": "Fluxomics",
        "subtitle": "Stable-isotope tracing LC-MS metabolomics",
        "chips": ["13C tracing", "Compound Discoverer", "Generic table"],
        "highlights": [
            "Natural-abundance correction and IQR normalization",
            "Isotopologue fractions, incorporation, enrichment",
            "Statistics, bar plots, heatmaps, label arrival",
            "25 curated pathway maps with Z-score cards",
        ],
    },
    # Integration engine, not a data type: "integration" keeps it out of the platforms a Single- or Duo-Omics plan
    # can choose, and out of the "five platforms" counts; plans that cover all platforms include it.
    "integrix": {
        "folder": "integrix",
        "marker": "trace_core",
        "integration": True,
        "name": "Integrix",
        "subtitle": "Tri-omics integration: transcriptomics · proteomics · metabolomics",
        "chips": ["Gene → protein → enzyme → metabolite", "Paired or partial samples"],
        "highlights": [
            "Reaction-level coherence of enzyme and metabolite changes",
            "Shared and omics-specific pathways, IPA-style pathway maps",
            "Bottlenecks, regulators and prioritized targets",
            "MOFA+ factors and a cross-validated integrated signature",
        ],
    },
}

# Platforms that analyse one data type (the five platforms); integration engines are listed separately.
DATA_PLATFORMS = [k for k, p in PLATFORMS.items() if not p.get("integration")]
INTEGRATION_PLATFORMS = [k for k, p in PLATFORMS.items() if p.get("integration")]
