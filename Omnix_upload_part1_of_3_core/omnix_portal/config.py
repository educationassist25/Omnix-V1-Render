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
# platforms the plan unlocks ("one" or "all"); "keys" is the number of access keys issued (users).
# To issue keys for a plan: python tools/make_access_key.py --plan "<name>" --platforms ... --count ...
SUBSCRIPTION_INTRO = ("Every visitor can run all Omnix analyses on the built-in demo datasets, free. A subscription "
                      "unlocks the same workflows for your own data. Choose one platform or all three, for yourself "
                      "or for your whole lab, billed monthly or annually.")
CURRENCY = "$"
SUBSCRIPTION_PLANS = [
    {
        "id": "single", "name": "Single-Omics", "audience": "One researcher, one data type",
        "monthly": 49, "annual": 490, "platforms": "one", "keys": 1,
        "features": [
            "Your own data in one platform of your choice",
            "Metabolomics, Proteomics or Transcriptomics",
            "All analyses, figures and tables in that platform",
            "Demo access to the other two platforms",
            "Email support",
        ],
    },
    {
        "id": "multi", "name": "Multi-Omics", "audience": "One researcher, all three data types",
        "monthly": 99, "annual": 990, "platforms": "all", "keys": 1, "badge": "Most popular",
        "features": [
            "Your own data in all three platforms",
            "Metabolomics, Proteomics and Transcriptomics",
            "One consistent workflow across data types",
            "All analyses, figures and tables",
            "Priority email support",
        ],
    },
    {
        "id": "lab", "name": "Lab", "audience": "Research groups and core facilities",
        "monthly": None, "annual": 3900, "platforms": "all", "keys": 10,
        "features": [
            "Up to 10 users, each with a personal access key",
            "All three platforms",
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
            "All three platforms",
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
]
PRICING_FAQ = [
    ("What can I do without a subscription?",
     "Everything except uploading your own files: every analysis in all three platforms runs on the built-in demo "
     "datasets, and all figures and tables from the demo can be downloaded."),
    ("How do I get access after subscribing?",
     "You receive a personal access key (OMX-XXXX-XXXX-XXXX-XXXX). Enter it below under Activate your subscription; "
     "uploads unlock immediately in the platforms your plan covers."),
    ("Can I change platforms or upgrade later?",
     "Yes. Single-Omics can be upgraded to Multi-Omics or Lab at any time; you pay only the difference for the rest "
     "of your billing period."),
    ("Can a key be shared?",
     "Each key is for one person. The Lab plan gives each member of your group their own key."),
    ("What happens when a subscription ends?",
     "Your key stops unlocking uploads after its end date, and Omnix returns to demo access. Renewing reactivates "
     "the same key."),
]

# ---- The three analysis platforms -------------------------------------------------------
# key -> (folder under apps/, menu name, one-line description, chips, highlights)
PLATFORMS = {
    "metabolomics": {
        "folder": "metabolomics",
        "name": "Metabolomics",
        "subtitle": "LC-MS metabolomics and lipidomics",
        "chips": ["Untargeted", "Targeted", "Lipidomics"],
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
}
