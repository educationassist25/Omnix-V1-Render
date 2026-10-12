"""
MetaboAI Pro Metabolomics — the Metabolomics app (Untargeted Metabolomics and Targeted Metabolomics).

    streamlit run metabolomics_app.py

Runs the shared MetaboAI Pro pipeline (app.py) with only the Metabolomics analysis types offered;
every page and analysis is identical to app.py.
"""
import os
import runpy

runpy.run_path(os.path.join(os.path.dirname(os.path.abspath(__file__)), "app.py"),
               init_globals={"APP_FAMILY": "Metabolomics"}, run_name="__main__")
