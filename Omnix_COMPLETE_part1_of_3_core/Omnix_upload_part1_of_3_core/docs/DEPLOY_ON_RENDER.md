# Deploying Omnix on Render

The project is ready for Render: `render.yaml` defines the web service and `render_start.sh` starts it.
Your settings (sign-in, Stripe, owner email) go into one **Secret File** on Render, never into GitHub.

## What you need before you start

| # | Item | Where it comes from | Needed for |
|---|------|--------------------|-----------|
| 1 | GitHub account and a repository with the Omnix files (`omnix_app.py`, `render.yaml`, `render_start.sh`, `requirements.txt`, `omnix_portal/`, `apps/`) at its top level | github.com | Required |
| 2 | Render account, signed up with GitHub, with a payment card | render.com | Required |
| 3 | The site address: `https://omnix.onrender.com` (free) or your own domain | You decide; your own domain needs access to its DNS settings | Required |
| 4 | Your owner email(s) | You | Owner console |
| 5 | A cookie secret (a long random string) | `python -c "import secrets; print(secrets.token_hex(32))"` | Sign-in |
| 6 | Google OAuth client ID and client secret | console.cloud.google.com → APIs & Services → Credentials | Google sign-in |
| 7 | Microsoft client ID and secret (optional) | portal.azure.com → App registrations | Microsoft sign-in |
| 8 | Stripe restricted key and the five price IDs (optional at launch) | dashboard.stripe.com | Online payments |
| 9 | Gmail App Password (optional) | myaccount.google.com/apppasswords | Contact form sending email |
| 10 | Access-key hashes (optional) | `python tools/make_access_key.py ...` | Keys for purchase-order customers |

Items 1 to 6 are enough to launch; add 7 to 10 whenever you are ready.

## Step by step

### 1. Put the code on GitHub
Create a repository (private is fine) and upload the three Omnix parts so that `omnix_app.py` and `render.yaml` are at the top
level. Do not upload a `secrets.toml` (it is in `.gitignore`).

### 2. Create the service from the Blueprint
1. On render.com, go to **New → Blueprint** and connect your GitHub account.
2. Pick the Omnix repository. Render reads `render.yaml` and shows one web service, **omnix**, on the **Standard** plan (2 GB).
3. Click **Apply**. The first build installs the packages and takes about 5 to 10 minutes.

The Standard plan is needed because one full analysis uses about 450 MB of memory. The Free and Starter plans have 512 MB, and
the Free plan also sleeps after 15 minutes without visitors.

### 3. Add your settings as a Secret File
1. Open the **omnix** service → **Environment** → **Secret Files** → **Add Secret File**.
2. Filename: `secrets.toml`. Contents: start from `.streamlit/secrets.toml.example` and fill in:
   - `[omnix]` `owner_emails` (item 4) and `app_url = "https://omnix.onrender.com"` (item 3)
   - `[auth]` `redirect_uri = "https://omnix.onrender.com/oauth2callback"` and `cookie_secret` (item 5)
   - `[auth.google]` client ID and secret (item 6)
   - `[stripe]`, `[omnix.smtp]`, `[omnix.access_keys]` when you have them (items 8 to 10)
3. Click **Save Changes**. Render redeploys, and the log shows `Omnix: secrets.toml loaded from Render Secret Files.`

### 4. Tell Google about the new address
In Google Cloud Console → **APIs & Services → Credentials** → your OAuth client, add this **Authorized redirect URI**:
`https://omnix.onrender.com/oauth2callback`, then save. The address must match `redirect_uri` exactly. If you use Stripe, also
check that `app_url` matches the address.

### 5. Test the live site
1. Open `https://omnix.onrender.com`: the Home page appears.
2. Run each platform's demo.
3. Click **Login**, sign in with Google using your owner email, and check the **Owner console**.
4. With a Stripe test key, buy a plan with card `4242 4242 4242 4242`, then switch to the live key.

### 6. Your own domain (optional)
1. In the service, open **Settings → Custom Domains → Add** and enter, for example, `www.omnix-analytics.com`.
2. At your domain registrar, add the DNS record Render shows (a CNAME to `omnix.onrender.com`). Render adds the HTTPS
   certificate itself.
3. Change `app_url` and `redirect_uri` in the Secret File to the new address, and add the new redirect URI in Google.

### 7. Updating the site
Push changes to the `main` branch: Render rebuilds and redeploys automatically. **Events** lists every deploy and can roll back
to an earlier one; **Logs** shows errors.

## Notes
- Omnix stores nothing on the server: subscriptions live in Stripe and uploads stay in the visitor's session. The server's disk
  being reset on each deploy is therefore fine.
- Health check: Render checks `/_stcore/health` and restarts the service if it stops answering.
- Secrets shown in chats or emails (Google client secret, owner access key, cookie secret) should be replaced with new ones before
  going live.
