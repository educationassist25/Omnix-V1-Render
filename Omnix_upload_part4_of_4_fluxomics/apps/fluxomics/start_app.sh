#!/usr/bin/env bash
# Start Omnix Fluxomics (macOS / Linux). Usage: bash start_app.sh
cd "$(dirname "$0")"
if command -v conda >/dev/null 2>&1; then
  eval "$(conda shell.bash hook)"; conda activate fluxomics 2>/dev/null || conda activate base
fi
python -c "import streamlit" 2>/dev/null || pip install -r requirements.txt
streamlit run app.py
