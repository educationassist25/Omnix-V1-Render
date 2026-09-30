# Omnix accounts and online billing: setup guide

The **Login / Account** page lets customers sign in with Google or Microsoft, subscribe and pay
online, download invoices, update their card, cancel or resume, and (Lab plan) add team members.

* **Sign-in** uses Streamlit's built-in login. Omnix never sees or stores passwords.
* **Payments, invoices, renewals, receipts and failed-payment retries** are handled by **Stripe**.
  Card numbers are entered on Stripe's pages and never reach Omnix.
* **Nothing is stored by Omnix.** Each time someone signs in, their subscription is read from Stripe,
  and that decides which platforms accept their own data. No database or webhook is needed.

Access keys and the Owner console keep working alongside this, e.g. for institutions that pay by
purchase order.

Everything below goes into **Streamlit Cloud → your app → Settings → Secrets**. The template is
`.streamlit/secrets.toml.example`. Leave out any section you are not using.

---

## 1. Try it first (preview mode, 1 minute)

```toml
[omnix]
billing_preview = true
```

The Login page then shows **Preview with a demo account**. Subscriptions, invoices and the card are
simulated and nothing is charged. **Remove this line before going live.**

## 2. Your own owner sign-in

```toml
[omnix]
owner_emails = ["your.name@bcm.edu"]
app_url = "https://YOUR-APP.streamlit.app"
```

When you sign in with that email, you get full access and the **Owner console**. Its **Online
subscriptions** tab lists every customer, plan, status, renewal date and annual recurring revenue.

## 3. Sign-in with Google and/or Microsoft

Make a random cookie secret: `python -c "import secrets; print(secrets.token_hex(32))"`

```toml
[auth]
redirect_uri = "https://YOUR-APP.streamlit.app/oauth2callback"
cookie_secret = "the random string"
```

**Google** (console.cloud.google.com)
1. Go to *APIs & Services → OAuth consent screen*: choose External, app name *Omnix*, add your support email, publish.
2. Go to *Credentials → Create credentials → OAuth client ID → Web application*.
   Set the Authorized redirect URI to `https://YOUR-APP.streamlit.app/oauth2callback`.
3. Copy the client ID and secret:

```toml
[auth.google]
client_id = "....apps.googleusercontent.com"
client_secret = "..."
server_metadata_url = "https://accounts.google.com/.well-known/openid-configuration"
```

**Microsoft** (portal.azure.com → Microsoft Entra ID → App registrations)
1. Click *New registration*. Name it *Omnix*. For account types, choose *Accounts in any organizational
   directory and personal Microsoft accounts*. Set the Web redirect URI to
   `https://YOUR-APP.streamlit.app/oauth2callback`.
2. Go to *Certificates & secrets → New client secret*.

```toml
[auth.microsoft]
client_id = "Application (client) ID"
client_secret = "..."
server_metadata_url = "https://login.microsoftonline.com/common/v2.0/.well-known/openid-configuration"
```

If Microsoft sign-in reports an issuer error, replace `common` with a specific tenant ID. Note that
this limits sign-in to that one organization. Google alone is enough to start.

## 4. Stripe payments

1. **Create and activate an account** at stripe.com: business details and a bank account for payouts.
   Stripe will ask for a public website with your **contact details, terms, and refund/cancellation
   policy**, so put these on the Contact page or link them.
2. **Build and test everything in Test mode first** (toggle at the top of the dashboard).
3. **Create the products and prices** (*Product catalog → Add product*, recurring prices, USD):

   | Product | Price 1 | Price 2 |
   |---|---|---|
   | Omnix Single-Omics | $49 / month | $490 / year |
   | Omnix Multi-Omics | $99 / month | $990 / year |
   | Omnix Lab | $3,900 / year | — |

   Copy each price ID (`price_...`).
4. **Set up the customer portal** (*Settings → Billing → Customer portal*). Turn on:
   * invoice history
   * update payment methods
   * update billing address and tax ID
   * cancel subscriptions at the end of the billing period
   * switch plans (add the three products)
5. **Turn on the emails** (*Settings → Billing → Subscriptions and emails*):
   * send receipts for successful payments
   * send upcoming-renewal reminders
   * send failed-payment emails with Smart Retries
6. **Create a restricted API key** (*Developers → API keys → Create restricted key*) with these
   permissions:
   * **Write:** Customers, Subscriptions, Checkout Sessions, Customer portal
   * **Read:** Invoices, Payment methods, Prices, Products

   A restricted key cannot issue refunds or move money out, which limits the damage if it ever leaks.

```toml
[stripe]
secret_key = "rk_test_..."          # rk_live_... when you go live
automatic_tax = false               # true after you set up Stripe Tax

[stripe.prices]
single_monthly = "price_..."
single_annual  = "price_..."
multi_monthly  = "price_..."
multi_annual   = "price_..."
lab_annual     = "price_..."
```

7. **Run a test purchase:**
   1. Sign in on Omnix and choose a plan, then pay with card `4242 4242 4242 4242` (any future date, any CVC).
   2. You come back to the Account page, and uploads unlock.
   3. Check the invoice, the card, cancel and resume, and adding a Lab member.
8. **Go live:** create the same products in Live mode, then replace the key and price IDs with the live ones.

## Going-live checklist

- [ ] `billing_preview` removed (or `false`)
- [ ] `app_url` and `redirect_uri` use the real app address
- [ ] your email is in `owner_emails`
- [ ] live Stripe restricted key and live price IDs
- [ ] customer portal and emails configured in Live mode
- [ ] terms and refund policy published; Contact page details filled in (`omnix_portal/config.py`)
- [ ] one real purchase and refund made to confirm everything

## What happens when…

| Situation | Result |
|---|---|
| A card payment fails at renewal | Stripe retries and emails the customer. The status shows *Payment due* and access continues during the retries. Access stops if the subscription is finally cancelled. |
| A customer cancels | Access continues to the end of the paid period, and they can resume until then. |
| A customer switches plan or billing period | They do it in the Stripe portal, with prorated charges. Omnix picks up the change within two minutes. |
| An institution wants to pay by purchase order | Send a Stripe invoice from the dashboard, or give them an access key or a `subscriber_emails` entry. |
| A refund is needed | Refund in the Stripe dashboard, and cancel the subscription there if access should stop. |
| A Lab member leaves | The Lab owner removes them on Account → Team. |
