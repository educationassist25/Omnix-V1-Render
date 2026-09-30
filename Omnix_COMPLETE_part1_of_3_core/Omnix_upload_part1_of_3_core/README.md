# Omnix™ — Integrated Multi-Omics Data Analytics

One site holding three analysis platforms behind a shared Omnix header:

| Menu | Platform | What it analyses |
|---|---|---|
| **Metabolomics** | `apps/metabolomics` | LC-MS metabolomics and lipidomics (untargeted and targeted, single or combined datasets) |
| **Proteomics** | `apps/proteomics` | LC-MS/MS proteomics (Label-Free and TMT) |
| **Transcriptomics** | `apps/transcriptomics` | Bulk RNA-seq (raw counts with DESeq2, or logCPM) |

plus **Subscription**, **Contact** (menu and page) and **Help** (a detailed tutorial for each platform).

## Run

```bash
pip install -r requirements.txt
streamlit run omnix_app.py
```

Pages can be linked directly: `?page=metabolomics`, `?page=proteomics`, `?page=transcriptomics`,
`?page=subscription`, `?page=contact`, `?page=help-metabolomics` (and `help-proteomics`, `help-transcriptomics`).
To host it on **Render** (recommended), follow `docs/DEPLOY_ON_RENDER.md`: `render.yaml` and `render_start.sh` are ready. On Streamlit Community Cloud, set the main file to `omnix_app.py`.

Each platform still runs on its own too: `streamlit run apps/metabolomics/app.py` (likewise for the others).

## Your details

Edit **`omnix_portal/config.py`**: contact email, phone, organization and address, the Subscription page's
text and plans, and the description shown on each platform's card. Everything in `[SQUARE BRACKETS]` is a
placeholder shown on the site until you replace it. The tutorials are Markdown files in
`omnix_portal/tutorials/` and can be edited directly.

## Subscriptions (demo access vs. your own data)

Without a subscription, visitors can run every analysis in all three platforms on the built-in demo
datasets (and download all results); every file upload is locked. Subscribers unlock uploads with a
personal access key.

1. Create a key for a subscriber: `python tools/make_access_key.py --name "Jane Doe" --plan "Standard" --expires 2027-09-30`
2. Send them the printed key (`OMX-....`), and paste the printed line into the app's secrets
   (Streamlit Community Cloud: *App settings → Secrets*; locally `.streamlit/secrets.toml`), under `[omnix.access_keys]`.
   See `.streamlit/secrets.toml.example` for the full format. Only a hash of the key is stored.
3. The subscriber enters the key on the **Subscription** page. To end access, delete the line or let `expires` pass.

Plans (prices and texts in `omnix_portal/config.py`):

| Plan | Your own data in | Users | Academic list price |
|---|---|---|---|
| Single-Omics | one platform | 1 | $49 / month or $490 / year |
| Multi-Omics | all three | 1 | $99 / month or $990 / year |
| Lab | all three | up to 10 | $3,900 / year |
| Enterprise | all three | as agreed | custom quote |

Issue keys to match the plan: Single-Omics `--platforms proteomics` (the key then unlocks uploads only in
that platform; the others stay in demo mode with an upgrade prompt), Lab `--count 10` (one key per member).

### Accounts and online billing (Login / Account page)

Customers sign in with Google or Microsoft, subscribe and pay through Stripe, download invoices,
update their card, cancel or resume, and add Lab team members. Stripe holds all billing data; Omnix
reads each person's subscription when they sign in. Setup, step by step: `docs/ACCOUNT_AND_BILLING_SETUP.md`.
To try it without Stripe, set `billing_preview = true` under `[omnix]` (simulated, never on the live site).

### Your owner key and the Owner console

1. On your own computer: `python tools/make_access_key.py --owner --name "Your Name"`.
2. Paste the printed line into **Settings → Secrets** under `[omnix.access_keys]` and save.
3. Enter the key on the Subscription page. You get full access plus **Owner console** in the side panel:
   subscriber list with expiry status (CSV export), key creation for any plan, renew/revoke instructions,
   key lookup, and a monitor of wrong-key attempts.

The app cannot edit its own secrets, so keys created in the console take effect once you paste their lines
into Secrets. Keys are shown once and never stored; only their SHA-256 hashes are kept. After 5 wrong keys
within 15 minutes, activation pauses for that visitor for 15 minutes. `.gitignore` keeps a local
`.streamlit/secrets.toml` out of GitHub.

Optional: with Streamlit authentication configured (`[auth]` in secrets), signed-in users listed under
`[omnix.subscriber_emails]` are subscribers too. `open_access = true` under `[omnix]` unlocks everything
(internal use or testing).

## How it is put together

- `omnix_app.py` — entry point; `omnix_portal/portal.py` — header, menus, portal pages, and runs the
  selected platform's `app.py` inside the portal.
- Each platform keeps its own code package (`metabolomics_modules`, `proteomics_modules`,
  `transcriptomics_modules`) so they never load each other's code.
- `omnix_portal/isolation.py` keeps each platform's data separate within one browser session: moving from
  Metabolomics to Proteomics sets Metabolomics' data aside and brings it back, untouched, when you return.
  Widget selections on a platform's page reset when you come back; loaded data and results do not.

## Tests

```bash
python -m pytest tests -q                      # Omnix: pages, menus, deep links, isolation, full workflows
(cd apps/metabolomics && python -m pytest tests -q)
(cd apps/proteomics && python test_pipeline.py)
```
