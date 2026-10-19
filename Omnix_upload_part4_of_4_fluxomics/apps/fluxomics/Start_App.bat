@echo off
REM Double-click to start Omnix Fluxomics (Windows, Anaconda).
cd /d "%~dp0"
call conda activate fluxomics 2>nul || call conda activate base
python -c "import streamlit" 2>nul || pip install -r requirements.txt
streamlit run app.py
pause
