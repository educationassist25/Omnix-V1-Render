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
exec streamlit run omnix_app.py \
  --server.port "${PORT:-10000}" \
  --server.address 0.0.0.0 \
  --server.headless true \
  --browser.gatherUsageStats false
