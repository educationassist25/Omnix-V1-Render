# Omnix: speed and "Connection failed with status 429" fix (2026-11-09)

Only two files change, both in `Omnix_upload_part1_of_3_core/`.

## requirements.txt
`streamlit[auth]>=1.52` let every Render build install the newest Streamlit (1.67). From 1.58 Streamlit runs on a
new web server (uvicorn) and splits its page code into many more files: opening Integrix made 115 requests to the
server with 1.67 against 36 with 1.52.2 (home page 73 against 14). On Render each first visit sends that burst
through Render's request protection, which answers "429 Too Many Requests" and makes pages slow.
Now pinned: `streamlit[auth]==1.52.2` (Tornado server, the version Omnix was built on).

## render_start.sh
- Caps NumPy/BLAS, scikit-learn, PyDESeq2 and MOFA+ threads at 2 (Render gives the service 2 CPUs, but the libraries
  see all of the host's CPUs and start that many threads).
- `MALLOC_ARENA_MAX=2`: less memory fragmentation in the long-running server.
- `--server.fileWatcherType none --server.runOnSave false`: no source-file watching on the live site.

## Checks (Streamlit 1.52.2)
- Browser: home, all six platforms, Subscription (billing toggle), Account, Help, Data preparation, Fluxomics demo,
  Proteomics demo load: no errors. Integrix demo run: 27 s, server stayed responsive, peak memory about 0.9 GB.
- Automated tests: 53 of 74 pass. The 21 that fail are limits of the older Streamlit test tool (it cannot operate the
  Annual/Monthly switch or a hidden radio button), not app errors; the same screens work in the browser.

## After uploading
Render redeploys; in the build log check "Successfully installed ... streamlit-1.52.2".
If you start the service with your own Start Command in the Render dashboard, set it to `bash render_start.sh`.
