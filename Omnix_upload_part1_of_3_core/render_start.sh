#!/usr/bin/env bash
# Start Omnix on Render.
# Render keeps the Secret File "secrets.toml" at /etc/secrets/secrets.toml; Streamlit reads
# .streamlit/secrets.toml, so it is copied there on every start (never committed to GitHub).
set -euo pipefail
mkdir -p .streamlit
if [ -f /etc/secrets/secrets.toml ]; then
  cp /etc/secrets/secrets.toml .streamlit/secrets.toml
  echo "Omnix: secrets.toml loaded from Render Secret Files."
else
  echo "Omnix: no Secret File 'secrets.toml' found - the site runs in demo mode without sign-in or payments."
fi
# Render gives the service 2 CPUs, but the libraries see every CPU of the host machine and start that many threads
# (NumPy/BLAS, scikit-learn, PyDESeq2, MOFA+). Too many threads on 2 CPUs slow every platform down, so cap them.
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMEXPR_MAX_THREADS=2 VECLIB_MAXIMUM_THREADS=2
export LOKY_MAX_CPU_COUNT=2 NUMBA_NUM_THREADS=2
export MALLOC_ARENA_MAX=2        # less memory lost to fragmentation in a long-running multi-user server
# Integrix "download this figure as shown": Plotly's image engine (Kaleido) draws the chart with a headless Chrome.
# Fetch it once in the background (about 100 MB, a few seconds; skipped when already there) so start-up is not
# delayed; until it is ready, Integrix offers the publication-style figure instead.
( python -c "import kaleido; kaleido.get_chrome_sync()" > /tmp/omnix_chrome.log 2>&1 \
    && echo "Omnix: Chrome for figure export ready." || echo "Omnix: Chrome for figure export not available (see /tmp/omnix_chrome.log)." ) &
exec streamlit run omnix_app.py \
  --server.port "${PORT:-10000}" \
  --server.address 0.0.0.0 \
  --server.headless true \
  --browser.gatherUsageStats false \
  --server.fileWatcherType none \
  --server.runOnSave false
