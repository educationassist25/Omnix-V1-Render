"""
The Account page (menu: Login / Account): sign in or register, then manage the subscription,
invoices, payment method and Lab team. Billing itself lives in Stripe (omnix_portal/billing.py).
"""

import datetime as dt
import html
import re

import pandas as pd
import streamlit as st

from . import access, account, billing, config

TEAL, TEAL_DARK, ORANGE, INK, INK2, LINE = "#00695C", "#00574D", "#E65100", "#111816", "#4A5754", "#DCE3E0"
EMAIL_RE = re.compile(r"^[^@\s]+@[^@\s]+\.[^@\s]+$")

_ICON = {
    "renew": '<polyline points="23 4 23 10 17 10"/><polyline points="1 20 1 14 7 14"/>'
             '<path d="M3.51 9a9 9 0 0 1 14.85-3.36L23 10M1 14l4.64 4.36A9 9 0 0 0 20.49 15"/>',
    "invoice": '<path d="M14 2H6a2 2 0 0 0-2 2v16a2 2 0 0 0 2 2h12a2 2 0 0 0 2-2V8z"/><polyline points="14 2 14 8 20 8"/>'
               '<line x1="16" y1="13" x2="8" y2="13"/><line x1="16" y1="17" x2="8" y2="17"/>',
    "card": '<rect x="1" y="4" width="22" height="16" rx="2" ry="2"/><line x1="1" y1="10" x2="23" y2="10"/>',
    "team": '<path d="M17 21v-2a4 4 0 0 0-4-4H5a4 4 0 0 0-4 4v2"/><circle cx="9" cy="7" r="4"/>'
            '<path d="M23 21v-2a4 4 0 0 0-3-3.87"/><path d="M16 3.13a4 4 0 0 1 0 7.75"/>',
    "support": '<path d="M21 11.5a8.38 8.38 0 0 1-.9 3.8 8.5 8.5 0 0 1-7.6 4.7 8.38 8.38 0 0 1-3.8-.9L3 21l1.9-5.7a8.38 '
               '8.38 0 0 1-.9-3.8 8.5 8.5 0 0 1 4.7-7.6 8.38 8.38 0 0 1 3.8-.9h.5a8.48 8.48 0 0 1 8 8v.5z"/>',
    "lock": '<rect x="3" y="11" width="18" height="11" rx="2" ry="2"/><path d="M7 11V7a5 5 0 0 1 10 0v4"/>',
}


def _icon(name, size=40, color=TEAL):
    return (f'<svg width="{size}" height="{size}" viewBox="0 0 24 24" fill="none" stroke="{color}" stroke-width="1.6" '
            f'stroke-linecap="round" stroke-linejoin="round">{_ICON[name]}</svg>')


def css():
    return f"""
[class*="st-key-omnix_nav_account"] [data-testid="stButton"] button {{ background: {TEAL} !important; border: 0 !important;
  border-radius: 999px !important; padding: 3px 16px 3px 12px !important; min-height: 0 !important; margin-bottom: 9px; }}
[class*="st-key-omnix_nav_account"] [data-testid="stButton"] button p {{ font-size: 16.5px !important; }}
[class*="st-key-omnix_header"] [class*="st-key-omnix_nav"] {{ column-gap: 28px !important; }}
[class*="st-key-omnix_nav_account"] [data-testid="stButton"] button p,
[class*="st-key-omnix_nav_account"] [data-testid="stButton"] button span {{ color: #FFFFFF !important; font-weight: 600 !important; }}
[class*="st-key-omnix_nav_account"] [data-testid="stButton"] button:hover {{ background: {TEAL_DARK} !important; }}
.omx-hero {{ border-radius: 18px; padding: 54px 28px; text-align: center; color: #FFFFFF; margin: 8px 0 34px 0;
  background: linear-gradient(115deg, #0B4F6C 0%, {TEAL} 48%, #3FA36B 100%); position: relative; overflow: hidden; }}
.omx-hero:after {{ content: ""; position: absolute; inset: 0; opacity: .18;
  background: radial-gradient(1200px 300px at 20% 120%, #FFFFFF 0%, transparent 60%),
              radial-gradient(900px 260px at 90% -30%, #FFFFFF 0%, transparent 55%); }}
.omx-hero .t {{ font-family: 'Sora', system-ui, sans-serif; font-weight: 300; font-size: 40px; letter-spacing: 12px;
  text-transform: uppercase; position: relative; z-index: 1; }}
.omx-hero .s {{ font-size: 17px; opacity: .92; margin-top: 12px; letter-spacing: .3px; position: relative; z-index: 1; }}
.omx-hero.small {{ padding: 30px 28px; text-align: left; }}
.omx-hero.small .t {{ font-size: 28px; letter-spacing: 1px; text-transform: none; font-weight: 600; }}
.omx-hero.small .s {{ margin-top: 6px; font-size: 15.5px; }}
.omx-tile {{ text-align: center; padding: 8px 14px 0 14px; }}
.omx-tile h4 {{ white-space: nowrap; font-family: 'Sora', system-ui, sans-serif; font-size: 19px !important; margin: 14px 0 8px 0 !important;
  padding: 0 !important; color: {INK}; }}
.omx-tile p {{ color: {INK2}; font-size: 15px; letter-spacing: .4px; line-height: 1.6; margin: 0; }}
[class*="st-key-omnix_quiet_haskey"] {{ max-width: 460px; margin: 10px auto 0 auto; }}
[class*="st-key-omnix_signin_card"] {{ max-width: 460px; margin: 18px auto 0 auto; background: #FFFFFF;
  border-radius: 16px !important; padding: 10px 8px 18px 8px; box-shadow: 0 10px 30px rgba(17,24,22,.07); }}
.omx-signin-h {{ font-family: 'Sora', system-ui, sans-serif; font-size: 21px; font-weight: 700; text-align: center;
  margin: 6px 0 4px 0; color: {INK}; }}
.omx-signin-s {{ text-align: center; color: {INK2}; font-size: 14.5px; margin-bottom: 8px; }}
.omx-fine {{ text-align: center; color: {INK2}; font-size: 13px; margin-top: 6px; line-height: 1.5; }}
[class*="st-key-omnix_signin_"] [data-testid="stButton"] button {{ min-height: 48px; border-radius: 10px !important;
  border: 1px solid {LINE} !important; background: #FFFFFF !important; }}
[class*="st-key-omnix_signin_"] [data-testid="stButton"] button p {{ font-weight: 600 !important; color: {INK} !important;
  font-size: 15.5px !important; }}
[class*="st-key-omnix_signin_"] [data-testid="stButton"] button:hover {{ border-color: {TEAL} !important; }}
.omx-sub h3 {{ margin: 0 !important; padding: 0 !important; font-size: 22px !important; }}
.omx-pill {{ display: inline-block; border-radius: 999px; padding: 3px 11px; font-size: 12.5px; font-weight: 700;
  letter-spacing: .3px; margin-left: 10px; vertical-align: middle; }}
.omx-pill.ok {{ background: #E3F4EC; color: #11734B; }} .omx-pill.warn {{ background: #FFF1E0; color: #A64B00; }}
.omx-pill.bad {{ background: #FDE7E7; color: #B42318; }} .omx-pill.off {{ background: #EEF1F0; color: {INK2}; }}
.omx-facts {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(170px, 1fr)); gap: 14px; margin: 16px 0 6px 0; }}
.omx-facts .k {{ font-size: 12.5px; text-transform: uppercase; letter-spacing: .6px; color: {INK2}; font-weight: 700; }}
.omx-facts .v {{ font-size: 16px; color: {INK}; margin-top: 3px; }}
.omx-cardline {{ display: flex; align-items: center; gap: 16px; font-size: 17px; }}
.omx-cardline .n {{ font-family: 'IBM Plex Mono', ui-monospace, monospace; letter-spacing: 1px; }}
.omx-sum ul {{ min-height: 0 !important; }}
.omx-preview {{ background: #FFF6E5; border: 1px solid #F3D19C; color: #7A4A00; border-radius: 10px;
  padding: 10px 14px; font-size: 14.5px; margin-bottom: 14px; }}
"""


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _money(amount, currency="USD"):
    return f"${amount:,.2f}" if currency.upper() == "USD" else f"{amount:,.2f} {currency.upper()}"


def _d(day):
    return day.strftime("%b %d, %Y") if day else "—"


def _pill(sub):
    st_ = sub["status"]
    if st_ in ("active", "trialing") and sub["cancel_at_period_end"]:
        return f"<span class='omx-pill warn'>Ends {_d(sub['renews_on'])}</span>"
    return {"active": "<span class='omx-pill ok'>Active</span>", "trialing": "<span class='omx-pill ok'>Trial</span>",
            "past_due": "<span class='omx-pill bad'>Payment due</span>",
            "unpaid": "<span class='omx-pill bad'>Unpaid</span>",
            "paused": "<span class='omx-pill off'>Paused</span>"}.get(st_, "<span class='omx-pill off'>Ended</span>")


def _return_url():
    base = ""
    try:
        base = str((st.secrets.get("omnix", {}) or {}).get("app_url", "") or "")
    except Exception:
        pass
    if not base:
        try:
            base = str(st.context.url or "")
        except Exception:
            base = ""
    base = base.split("?")[0].rstrip("/") or "http://localhost:8501"
    return base + "/?page=account"


def _platforms_text(keys):
    from . import portal
    return "All three platforms" if len(keys) == len(billing.ALL) else portal._platform_names(keys)


def _refresh():
    billing.forget_cache()


def _cached(name, user, fn):
    """Per-session cache for Stripe reads (cleared with the billing cache after any change)."""
    import time
    key = f"omnix_bill_{name}"
    hit = st.session_state.get(key)
    if hit and time.time() - hit[0] < billing.CACHE_SECONDS:
        return hit[1]
    try:
        val = fn(user["email"])
    except Exception:
        val = None
    st.session_state[key] = (time.time(), val)
    return val


# ---------------------------------------------------------------------------
# signed out
# ---------------------------------------------------------------------------
def _signed_out():
    st.markdown("<div class='omx-hero'><div class='t'>Omnix Account</div>"
                "<div class='s'>Your subscription, invoices and team, in one place.</div></div>", unsafe_allow_html=True)
    tiles = [("renew", "Your subscription", "Subscribe, change plan, renew or cancel at any time."),
             ("invoice", "Invoices & payments", "Download invoices and receipts, pay online and update your card."),
             ("team", "Your lab team", "Give each member of your group their own access on the Lab plan."),
             ("support", "Priority support", "Reach the Omnix team directly about your data and analyses.")]
    cols = st.columns(4, gap="medium")
    for col, (icon, title, text) in zip(cols, tiles):
        col.markdown(f"<div class='omx-tile'>{_icon(icon, 46)}<h4>{title}</h4><p>{text}</p></div>",
                     unsafe_allow_html=True)

    provs = account.providers()
    with st.container(key="omnix_signin_card", border=True):
        st.markdown(f"<div style='text-align:center;margin-top:6px'>{_icon('lock', 30)}</div>"
                    "<div class='omx-signin-h'>Sign in or register</div>"
                    "<div class='omx-signin-s'>Use your university or work account.</div>", unsafe_allow_html=True)
        errs = st.session_state.pop("omnix_signin_error", None)
        if errs:
            st.error("Sign-in is not configured correctly:\n\n" + "\n".join(f"- {e}" for e in errs))
        if provs:
            for pid, label in provs.items():
                with st.container(key=f"omnix_signin_{pid or 'default'}"):
                    st.button(f"Continue with {label}", key=f"omnix_login_{pid or 'default'}",
                              on_click=account.sign_in, args=(pid,), use_container_width=True)
        if account.preview_enabled():
            with st.container(key="omnix_signin_preview"):
                st.button("Preview with a demo account", key="omnix_login_preview", on_click=account.sign_in,
                          args=("preview",), use_container_width=True)
        if not provs and not account.preview_enabled():
            st.info("Online accounts are being set up. Subscribers can use their access key on the Subscription page.")
        st.markdown("<div class='omx-fine'>New to Omnix? Signing in for the first time creates your account.<br>"
                    "Omnix never sees or stores your password.</div>", unsafe_allow_html=True)
    from . import portal
    with st.container(key="omnix_quiet_haskey"):
        st.button("I have an access key", key="omnix_acct_haskey",
                  on_click=portal.go, args=("subscription",), use_container_width=True)


# ---------------------------------------------------------------------------
# signed in
# ---------------------------------------------------------------------------
def _plan_chooser(be, user, heading, key):
    st.markdown(f"#### {heading}")
    options = [p for p in config.SUBSCRIPTION_PLANS if p["id"] in billing.SELF_SERVICE]
    pre = st.session_state.pop("omnix_acct_plan", None)
    if pre in [p["id"] for p in options]:
        st.session_state[f"omnix_acct_pick_{key}"] = pre
    c1, c2 = st.columns([3, 2], gap="large")
    with c1:
        pid = st.radio("Plan", [p["id"] for p in options], key=f"omnix_acct_pick_{key}", horizontal=True,
                       format_func=lambda i: billing.plan(i)["name"])
        p = billing.plan(pid)
        periods = ["Annual"] + (["Monthly"] if p.get("monthly") is not None else [])
        period = st.segmented_control("Billing period", periods, default="Annual", key=f"omnix_acct_per_{key}") or "Annual"
        interval = "year" if period == "Annual" else "month"
        platform = ""
        if p.get("platforms") == "one":
            platform = st.selectbox("Platform", list(billing.ALL), key=f"omnix_acct_plat_{key}",
                                    format_func=lambda k: config.PLATFORMS[k]["name"])
    with c2:
        amount = billing.price_of(pid, interval)
        per = "year" if interval == "year" else "month"
        feats = "".join(f"<li>{html.escape(f)}</li>" for f in p["features"][:4])
        st.markdown(f"<div class='omx-plan omx-sum' style='background:#fff;border:1px solid {LINE};border-radius:14px;"
                    f"padding:18px 20px'><h3>{html.escape(p['name'])}</h3>"
                    f"<div class='price'>{_money(amount)}<span> / {per}</span></div>"
                    f"<div class='per'>Renews automatically. Cancel any time; access continues to the end of the "
                    f"paid period.</div><ul>{feats}</ul></div>", unsafe_allow_html=True)
    if be.live:
        if st.button("Continue to secure checkout", key=f"omnix_acct_checkout_{key}", type="primary"):
            try:
                url = be.checkout(user["email"], user["name"], pid, interval, platform, _return_url())
                st.session_state["omnix_acct_checkout_url"] = url
            except Exception as e:
                st.error(str(e))
        url = st.session_state.get("omnix_acct_checkout_url")
        if url:
            st.link_button("Pay securely with Stripe  →", url, type="primary")
            st.caption("You will be taken to Stripe's secure payment page and brought back here afterwards.")
    else:
        if st.button("Subscribe (preview, no charge)", key=f"omnix_acct_checkout_{key}", type="primary"):
            be.checkout(user["email"], user["name"], pid, interval, platform, _return_url())
            _refresh()
            st.session_state["omnix_acct_flash"] = f"Subscribed to {p['name']}. Your first invoice is in Invoices."
            st.rerun()
    st.caption("Institution paying by purchase order, or need more than 10 users? Contact us for a quote.")


def _sub_card(be, user, sub, own=True):
    plan_line = sub["plan_name"] + ("" if own else " · team member")
    billing_line = (f"{_money(sub['amount'], sub['currency'])} / {'year' if sub['interval'] == 'year' else 'month'}"
                    if own else f"Paid by {sub.get('owner') or 'your lab'}")
    live = sub["status"] in billing.LIVE
    end_label = ("Access until" if sub["cancel_at_period_end"] or not live else "Renews on")
    with st.container(border=True, key=f"omnix_sub_{sub['id']}"):
        st.markdown(f"<div class='omx-sub'><h3>{html.escape(plan_line)} {_pill(sub)}</h3>"
                    f"<div class='omx-facts'>"
                    f"<div><div class='k'>Your data in</div><div class='v'>{html.escape(_platforms_text(sub['platforms']))}</div></div>"
                    f"<div><div class='k'>Billing</div><div class='v'>{html.escape(billing_line)}</div></div>"
                    f"<div><div class='k'>{end_label}</div><div class='v'>{_d(sub['renews_on'] if live else sub.get('canceled_at') or sub['renews_on'])}</div></div>"
                    f"<div><div class='k'>Users</div><div class='v'>{1 + len(sub['members']) if sub['plan'] == 'lab' else 1}"
                    f"{' of ' + str(sub['seats']) if sub['plan'] == 'lab' else ''}</div></div>"
                    f"</div></div>", unsafe_allow_html=True)
        if not own or not live:
            return
        if sub["status"] == "past_due":
            st.warning("Your last payment did not go through. Update your payment method to keep access.")
        confirm_key = f"omnix_acct_confirm_{sub['id']}"
        c = st.container(horizontal=True, gap="small")
        with c:
            if sub["cancel_at_period_end"]:
                if st.button("Resume subscription", key=f"omnix_resume_{sub['id']}", type="primary"):
                    be.set_cancel(user["email"], sub["id"], False)
                    _refresh()
                    st.session_state["omnix_acct_flash"] = "Your subscription will renew automatically again."
                    st.rerun()
            elif not st.session_state.get(confirm_key):
                if st.button("Cancel subscription", key=f"omnix_cancel_{sub['id']}"):
                    st.session_state[confirm_key] = True
                    st.rerun()
            if be.live:
                url = _portal_url(be, user)
                if url:
                    st.link_button("Change plan or billing", url)
        if st.session_state.get(confirm_key) and not sub["cancel_at_period_end"]:
            st.warning(f"Your subscription will not renew. You keep full access until {_d(sub['renews_on'])}, and "
                       "you can resume it any time before then.")
            c = st.container(horizontal=True, gap="small")
            with c:
                if st.button("Confirm cancellation", key=f"omnix_cancel_ok_{sub['id']}", type="primary"):
                    be.set_cancel(user["email"], sub["id"], True)
                    st.session_state.pop(confirm_key, None)
                    _refresh()
                    st.session_state["omnix_acct_flash"] = (f"Cancellation confirmed. Access continues until "
                                                            f"{_d(sub['renews_on'])}.")
                    st.rerun()
                if st.button("Keep my subscription", key=f"omnix_cancel_no_{sub['id']}"):
                    st.session_state.pop(confirm_key, None)
                    st.rerun()


def _portal_url(be, user):
    key = "omnix_bill_portal_url"
    if key not in st.session_state:
        try:
            st.session_state[key] = be.portal(user["email"], _return_url())
        except Exception:
            st.session_state[key] = None
    return st.session_state[key]


def _tab_subscription(be, user, snap):
    own_live = [s for s in snap["subs"] if s["status"] in billing.LIVE]
    teams = snap["teams"]
    if snap.get("error"):
        st.error("We could not reach the billing service just now. Please try again in a moment.")
    for s in own_live:
        _sub_card(be, user, s, own=True)
    for s in teams:
        _sub_card(be, user, s, own=False)
    if not own_live:
        if not teams:
            st.markdown("You don't have an active subscription. Without one, every analysis runs on the built-in "
                        "demo datasets.")
        _plan_chooser(be, user, "Choose a plan" if not snap["subs"] else "Renew or choose a new plan", "new")
    elif not be.live:
        with st.expander("Change plan"):
            _plan_chooser(be, user, "Switch to another plan", "change")
    ended = [s for s in snap["subs"] if s["status"] not in billing.LIVE]
    if ended:
        st.markdown("#### Previous subscriptions")
        st.dataframe(pd.DataFrame([{"Plan": s["plan_name"], "Your data in": _platforms_text(s["platforms"]),
                                    "Ended": _d(s.get("canceled_at") or s["renews_on"])} for s in ended]),
                     hide_index=True, use_container_width=True)


def _invoice_html(inv, user):
    c = config.CONTACT
    return f"""<!doctype html><html><head><meta charset="utf-8"><title>Invoice {inv['number']}</title>
<style>body{{font-family:Helvetica,Arial,sans-serif;color:#111816;max-width:760px;margin:40px auto;padding:0 24px}}
h1{{font-size:30px;margin:0;color:{TEAL_DARK}}} .m{{color:#4A5754;font-size:14px}} table{{width:100%;border-collapse:collapse;margin-top:28px}}
th,td{{text-align:left;padding:12px 8px;border-bottom:1px solid #DCE3E0}} th:last-child,td:last-child{{text-align:right}}
.tot td{{font-weight:700;border-bottom:0}} .row{{display:flex;justify-content:space-between;margin-top:26px}}
.paid{{display:inline-block;border:2px solid #11734B;color:#11734B;border-radius:6px;padding:4px 12px;font-weight:700;letter-spacing:1px}}</style>
</head><body><div class="row"><div><h1>{config.BRAND}™</h1><div class="m">{config.TAGLINE}</div></div>
<div style="text-align:right"><div style="font-size:22px;font-weight:700">INVOICE</div><div class="m">{inv['number']}</div>
<div class="m">Date: {_d(inv['date'])}</div></div></div>
<div class="row"><div><div class="m">Billed to</div><b>{html.escape(user['name'])}</b><br>{html.escape(user['email'])}</div>
<div style="text-align:right"><div class="m">From</div>{html.escape(c['organization'])}<br>{html.escape(c['address'])}<br>{html.escape(c['email'])}</div></div>
<table><tr><th>Description</th><th>Amount</th></tr><tr><td>{html.escape(inv['description'])}</td><td>{_money(inv['amount'], inv['currency'])}</td></tr>
<tr class="tot"><td>Total</td><td>{_money(inv['amount'], inv['currency'])}</td></tr></table>
<p style="margin-top:28px"><span class="paid">{inv['status'].upper()}</span></p>
<p class="m">Preview invoice generated by Omnix. In live mode invoices are issued by Stripe.</p></body></html>"""


def _tab_invoices(be, user):
    invs = _cached("invoices", user, be.invoices) or []
    if not invs:
        st.info("No invoices yet. Your invoices and receipts appear here after your first payment.")
        return
    status = {"paid": "Paid", "open": "Due", "void": "Void", "uncollectible": "Uncollectible"}
    table = pd.DataFrame([{"Invoice": i["number"], "Date": _d(i["date"]), "Description": i["description"],
                           "Amount": _money(i["amount"], i["currency"]), "Status": status.get(i["status"], i["status"]),
                           "PDF": i["pdf"], "View or pay": i["url"]} for i in invs])
    if be.live:
        st.dataframe(table, hide_index=True, use_container_width=True, column_config={
            "PDF": st.column_config.LinkColumn("PDF", display_text="Download"),
            "View or pay": st.column_config.LinkColumn("View or pay", display_text="Open")})
        due = [i for i in invs if i["status"] == "open" and i["url"]]
        if due:
            st.warning(f"{len(due)} invoice{'s' if len(due) > 1 else ''} awaiting payment.")
            st.link_button(f"Pay invoice {due[0]['number']}  →", due[0]["url"], type="primary")
    else:
        st.dataframe(table.drop(columns=["PDF", "View or pay"]), hide_index=True, use_container_width=True)
        pick = st.selectbox("Download an invoice", [i["number"] for i in invs], key="omnix_acct_inv_pick")
        inv = next(i for i in invs if i["number"] == pick)
        st.download_button(f"Download invoice {pick}", _invoice_html(inv, user).encode(), f"Omnix_invoice_{pick}.html",
                           "text/html", key="omnix_acct_inv_dl")
    st.caption("Need a quote, a W-9, or to pay by purchase order or bank transfer? Contact us and we will issue an "
               "invoice to your institution.")


def _tab_payment(be, user):
    pm = _cached("card", user, be.payment_method)
    with st.container(border=True):
        if pm:
            st.markdown(f"<div class='omx-cardline'>{_icon('card', 34)}<div><b>{html.escape(pm['brand'])}</b> "
                        f"<span class='n'>•••• {html.escape(pm['last4'])}</span><div style='font-size:14px;color:{INK2}'>"
                        f"Expires {html.escape(pm['exp'])} · Receipts go to {html.escape(pm['email'])}</div></div></div>",
                        unsafe_allow_html=True)
        else:
            st.markdown("No payment method on file. One is added when you subscribe.")
    if be.live:
        url = _portal_url(be, user)
        if url:
            st.link_button("Update payment method", url)
        st.caption("Card details are entered on and stored by Stripe (PCI DSS Level 1). Omnix never sees your card "
                   "number.")
    else:
        st.caption("Preview: in live mode, “Update payment method” opens Stripe's secure billing portal.")


def _tab_team(be, user, snap):
    labs = [s for s in snap["subs"] if s["plan"] == "lab" and s["status"] in billing.LIVE]
    if not labs:
        if snap["teams"]:
            st.markdown(f"You are a member of {html.escape(snap['teams'][0].get('owner') or 'a')} Lab plan.")
        st.info("Team access is part of the Lab plan: up to 10 people, each signing in with their own account.")
        return
    lab = labs[0]
    members = list(lab["members"])
    free = lab["seats"] - 1 - len(members)
    st.markdown(f"**{1 + len(members)} of {lab['seats']} users.** Members sign in with the email address you add "
                "here and get the same access as you.")
    rows = [{"Member": user["email"], "Role": "Account owner"}] + [{"Member": m, "Role": "Member"} for m in members]
    st.dataframe(pd.DataFrame(rows), hide_index=True, use_container_width=True)
    if free > 0:
        with st.form("omnix_team_add", border=True, clear_on_submit=True):
            email = st.text_input("Add a member by email", placeholder="colleague@university.edu",
                                  key="omnix_team_email")
            add = st.form_submit_button("Add member")
        if add:
            email = (email or "").strip().lower()
            if not EMAIL_RE.match(email):
                st.error("Enter a valid email address.")
            elif email == user["email"] or email in members:
                st.error("This person already has access.")
            else:
                be.set_members(user["email"], lab["id"], members + [email])
                _refresh()
                st.session_state["omnix_acct_flash"] = f"{email} now has access. Ask them to sign in with that email."
                st.rerun()
    else:
        st.info("All seats are in use. Remove a member or contact us to add seats.")
    if members:
        rm = st.selectbox("Remove a member", members, index=None, placeholder="Choose a member",
                          key="omnix_team_rm_pick")
        if rm and st.button(f"Remove {rm}", key="omnix_team_rm"):
            be.set_members(user["email"], lab["id"], [m for m in members if m != rm])
            _refresh()
            st.session_state["omnix_acct_flash"] = f"{rm} no longer has access."
            st.rerun()


def _tab_profile(be, user):
    s = access.status()
    via = {"sign-in": "Google / Microsoft sign-in", "preview": "Preview account"}.get(user["via"], user["via"])
    st.markdown(f"<div class='omx-facts'>"
                f"<div><div class='k'>Name</div><div class='v'>{html.escape(user['name'])}</div></div>"
                f"<div><div class='k'>Email</div><div class='v'>{html.escape(user['email'])}</div></div>"
                f"<div><div class='k'>Signed in with</div><div class='v'>{html.escape(via)}</div></div>"
                f"<div><div class='k'>Access</div><div class='v'>{html.escape(s.get('plan') or 'Demo')}</div></div>"
                f"</div>", unsafe_allow_html=True)
    st.markdown("#### Support")
    email = config.CONTACT["email"]
    st.markdown(f"Questions about your account, billing or an analysis? Email **{html.escape(email)}** from "
                f"{html.escape(user['email'])} and include your invoice number if it is about a payment.")
    from . import portal
    c = st.container(horizontal=True, gap="small")
    with c:
        for k in billing.ALL:
            st.button(f"{config.PLATFORMS[k]['name']} guide", key=f"omnix_acct_help_{k}", on_click=portal.go,
                      args=(f"help-{k}",))


def page():
    user = account.current_user()
    be = billing.backend()
    if not user:
        _signed_out()
        return
    first = user["name"].split(" ")[0]
    st.markdown(f"<div class='omx-hero small'><div class='t'>Welcome, {html.escape(first)}</div>"
                f"<div class='s'>{html.escape(user['email'])}</div></div>", unsafe_allow_html=True)
    if be is not None and not be.live:
        st.markdown("<div class='omx-preview'><b>Preview mode.</b> Subscriptions, invoices and cards on this page are "
                    "simulated; nothing is charged. Turn preview off before going live.</div>", unsafe_allow_html=True)
    outcome = st.query_params.get("checkout")
    if outcome:
        del st.query_params["checkout"]
        _refresh()
        if outcome == "success":
            st.success("Thank you. Your payment was received and your subscription is active.")
        else:
            st.info("Checkout was cancelled. Nothing was charged.")
    flash = st.session_state.pop("omnix_acct_flash", None)
    if flash:
        st.success(flash)
    with st.container(key="omnix_quiet_signout"):
        st.button("Sign out", key="omnix_signout", on_click=account.sign_out)
    if be is None:
        s = access.status()
        st.info("Online billing is not set up yet. "
                + ("Your access: " + s["plan"] + "." if s["subscribed"] else
                   "To subscribe, contact us or enter an access key on the Subscription page."))
        return
    snap = billing.snapshot(user["email"])
    t1, t2, t3, t4, t5 = st.tabs(["Subscription", "Invoices", "Payment method", "Team", "Profile & support"])
    with t1:
        _tab_subscription(be, user, snap)
    with t2:
        _tab_invoices(be, user)
    with t3:
        _tab_payment(be, user)
    with t4:
        _tab_team(be, user, snap)
    with t5:
        _tab_profile(be, user)
