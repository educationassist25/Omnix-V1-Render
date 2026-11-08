"""
MetaboAI Pro Lipidomics — the Lipidomics app (Untargeted Lipidomics and Targeted Lipidomics).

    streamlit run lipidomics_app.py

Runs the shared MetaboAI Pro pipeline (app.py) with only the Lipidomics analysis types offered;
every page and analysis is identical to app.py.
"""
import os
import runpy

runpy.run_path(os.path.join(os.path.dirname(os.path.abspath(__file__)), "app.py"),
               init_globals={"APP_FAMILY": "Lipidomics"}, run_name="__main__")
