"""
Online subscriptions, invoices and payments through Stripe.

Stripe is the record of truth: Omnix stores nothing itself. When a signed-in person opens Omnix,
their subscriptions are read from Stripe (cached for two minutes) and decide what they can upload.
Payments happen on Stripe's hosted Checkout page; changing the card, switching plan or billing
period, and downloading receipts happen in Stripe's hosted billing portal. Cancelling, resuming
and Lab team members are handled on the Omnix Account page.

Secrets (see docs/ACCOUNT_AND_BILLING_SETUP.md):

    [stripe]
    secret_key = "rk_live_..."             # a restricted key is best
    [stripe.prices]
    single_monthly = "price_..."           # one Stripe Price per plan and billing period
    single_annual  = "price_..."
    multi_monthly  = "price_..."
    multi_annual   = "price_..."
    lab_annual     = "price_..."

Without a Stripe key, `billing_preview = true` under [omnix] turns on a simulated backend for
trying the pages (no real payment). Otherwise online billing is simply off.
"""

import datetime as _dt
import hashlib
import threading
import time

import streamlit as st

from . import config

ALL = ("metabolomics", "proteomics", "transcriptomics")
SELF_SERVICE = ("single", "multi", "lab")            # Enterprise goes through sales
LIVE = ("active", "trialing", "past_due")             # statuses that give access (past_due: grace period)
CACHE_SECONDS = 120
LAB_SEATS = 10


def plan(plan_id):
    return next((p for p in config.SUBSCRIPTION_PLANS if p["id"] == plan_id), None)


def price_of(plan_id, interval):
    p = plan(plan_id) or {}
    return p.get("annual" if interval == "year" else "monthly")


def member_tag(email: str) -> str:
    return "m_" + hashlib.sha256(email.strip().lower().encode()).hexdigest()[:16]


def _date(ts):
    return _dt.datetime.fromtimestamp(int(ts), _dt.timezone.utc).date() if ts else None


def _sub_record(plan_id, interval, status, platform, period_end, cancel_at_end, amount, currency, sid,
                members=(), owner="", canceled_at=None):
    p = plan(plan_id) or {"name": plan_id or "Subscription", "platforms": "all"}
    platforms = [platform] if p.get("platforms") == "one" and platform in ALL else list(ALL)
    return {"id": sid, "plan": plan_id, "plan_name": p["name"], "interval": interval, "status": status,
            "platforms": platforms, "renews_on": period_end, "cancel_at_period_end": bool(cancel_at_end),
            "amount": amount, "currency": (currency or "usd").upper(), "members": list(members),
            "seats": LAB_SEATS if plan_id == "lab" else 1, "owner": owner, "canceled_at": canceled_at}


# ---------------------------------------------------------------------------
# Stripe
# ---------------------------------------------------------------------------
class StripeBilling:
    live = True

    def __init__(self, cfg):
        import stripe
        self.stripe = stripe
        stripe.api_key = cfg["secret_key"]
        stripe.max_network_retries = 2
        self.prices = dict(cfg.get("prices", {}))
        self.price_to_plan = {}
        for k, v in self.prices.items():
            pid, _, per = k.partition("_")
            self.price_to_plan[v] = (pid, "year" if per == "annual" else "month")
        self.automatic_tax = bool(cfg.get("automatic_tax", False))

    # customers
    def customer_id(self, email, name="", create=False):
        found = self.stripe.Customer.list(email=email, limit=1).data
        if found:
            return found[0].id
        if create:
            return self.stripe.Customer.create(email=email, name=name or None,
                                               metadata={"source": "omnix"}).id
        return None

    def _record(self, s, owner=""):
        item = s["items"]["data"][0] if s.get("items") and s["items"]["data"] else {}
        price = item.get("price") or {}
        md = dict(s.get("metadata") or {})
        pid, per = self.price_to_plan.get(price.get("id"), (md.get("omnix_plan", ""), None))
        per = per or ((price.get("recurring") or {}).get("interval") or "year")
        end = s.get("current_period_end") or item.get("current_period_end")   # newer API: on the item
        members = [m for m in (md.get("omnix_members") or "").split(",") if m]
        return _sub_record(pid, per, s["status"], md.get("omnix_platform", ""), _date(end),
                           s.get("cancel_at_period_end"), (price.get("unit_amount") or 0) / 100,
                           price.get("currency"), s["id"], members, owner, _date(s.get("canceled_at")))

    def subscriptions(self, email):
        cid = self.customer_id(email)
        if not cid:
            return []
        subs = self.stripe.Subscription.list(customer=cid, status="all", limit=20).data
        return [self._record(s, email) for s in subs if s["status"] not in ("incomplete", "incomplete_expired")]

    def memberships(self, email):
        """Lab subscriptions (of other people) that list this email as a team member."""
        try:
            found = self.stripe.Subscription.search(query=f"metadata['{member_tag(email)}']:'1'", limit=10).data
        except Exception:
            return []
        return [self._record(s) for s in found if s["status"] in LIVE]

    def invoices(self, email):
        cid = self.customer_id(email)
        if not cid:
            return []
        out = []
        for inv in self.stripe.Invoice.list(customer=cid, limit=24).data:
            if inv["status"] == "draft":
                continue
            line = (inv.get("lines") or {}).get("data") or [{}]
            out.append({"number": inv.get("number") or inv["id"], "date": _date(inv["created"]),
                        "description": line[0].get("description") or "Omnix subscription",
                        "amount": (inv.get("total") or 0) / 100, "currency": (inv.get("currency") or "usd").upper(),
                        "status": inv["status"], "pdf": inv.get("invoice_pdf"), "url": inv.get("hosted_invoice_url")})
        return out

    def payment_method(self, email):
        cid = self.customer_id(email)
        if not cid:
            return None
        cust = self.stripe.Customer.retrieve(cid, expand=["invoice_settings.default_payment_method"])
        pm = (cust.get("invoice_settings") or {}).get("default_payment_method")
        if not pm:
            pms = self.stripe.PaymentMethod.list(customer=cid, type="card", limit=1).data
            pm = pms[0] if pms else None
        if not pm or not pm.get("card"):
            return None
        c = pm["card"]
        return {"brand": str(c.get("brand", "card")).title(), "last4": c.get("last4", ""),
                "exp": f"{int(c.get('exp_month', 0)):02d}/{c.get('exp_year', '')}", "email": cust.get("email", email)}

    def checkout(self, email, name, plan_id, interval, platform, return_url):
        price = self.prices.get(f"{plan_id}_{'annual' if interval == 'year' else 'monthly'}")
        if not price:
            raise ValueError("This plan is not available for online purchase yet. Please contact us.")
        cid = self.customer_id(email, name, create=True)
        md = {"omnix_plan": plan_id, "omnix_platform": platform or "", "omnix_owner": email}
        kw = dict(mode="subscription", customer=cid, line_items=[{"price": price, "quantity": 1}],
                  subscription_data={"metadata": md}, metadata=md, allow_promotion_codes=True,
                  billing_address_collection="required", tax_id_collection={"enabled": True},
                  customer_update={"address": "auto", "name": "auto"},
                  success_url=return_url + "&checkout=success", cancel_url=return_url + "&checkout=cancel")
        if self.automatic_tax:
            kw["automatic_tax"] = {"enabled": True}
        return self.stripe.checkout.Session.create(**kw).url

    def portal(self, email, return_url):
        cid = self.customer_id(email)
        if not cid:
            return None
        return self.stripe.billing_portal.Session.create(customer=cid, return_url=return_url).url

    def _owned(self, email, sub_id):
        s = self.stripe.Subscription.retrieve(sub_id)
        if s["customer"] != self.customer_id(email):
            raise PermissionError("This subscription belongs to another account.")
        return s

    def set_cancel(self, email, sub_id, cancel: bool):
        self._owned(email, sub_id)
        self.stripe.Subscription.modify(sub_id, cancel_at_period_end=bool(cancel))

    def set_members(self, email, sub_id, members):
        s = self._owned(email, sub_id)
        md = {k: "" for k in (s.get("metadata") or {}) if k.startswith("m_")}      # "" deletes a key
        for m in members:
            md[member_tag(m)] = "1"
        md["omnix_members"] = ",".join(members)
        self.stripe.Subscription.modify(sub_id, metadata=md)

    def all_subscriptions(self):
        subs = self.stripe.Subscription.list(status="all", limit=100, expand=["data.customer"]).data
        return [self._record(s, (s.get("customer") or {}).get("email", "")) for s in subs
                if s["status"] not in ("incomplete", "incomplete_expired")]


# ---------------------------------------------------------------------------
# Preview (simulated; nothing is charged)
# ---------------------------------------------------------------------------
_store = {}                    # email -> {"subs": [...], "invoices": [...], "card": {...}}
_store_lock = threading.Lock()


class PreviewBilling:
    live = False

    def _acct(self, email):
        return _store.setdefault(email, {"subs": [], "invoices": [], "card": None})

    def customer_id(self, email, name="", create=False):
        return "cus_preview_" + hashlib.sha256(email.encode()).hexdigest()[:10]

    def subscriptions(self, email):
        with _store_lock:
            return [dict(s) for s in self._acct(email)["subs"]]

    def memberships(self, email):
        with _store_lock:
            return [dict(s) for acct in _store.values() for s in acct["subs"]
                    if email in s["members"] and s["status"] in LIVE]

    def invoices(self, email):
        with _store_lock:
            return [dict(i) for i in reversed(self._acct(email)["invoices"])]

    def payment_method(self, email):
        with _store_lock:
            card = self._acct(email)["card"]
            return dict(card, email=email) if card else None

    def checkout(self, email, name, plan_id, interval, platform, return_url):
        """Preview: subscribe at once with a test card and issue the first invoice."""
        amount = price_of(plan_id, interval)
        if amount is None:
            raise ValueError("This plan is not available for online purchase.")
        today = _dt.date.today()
        end = today + (_dt.timedelta(days=365) if interval == "year" else _dt.timedelta(days=30))
        with _store_lock:
            acct = self._acct(email)
            for s in acct["subs"]:                                    # replaces an existing plan
                if s["status"] in LIVE:
                    s["status"], s["canceled_at"] = "canceled", today
            sid = f"sub_preview_{len(acct['subs']) + 1}"
            acct["subs"].insert(0, _sub_record(plan_id, interval, "active", platform, end, False, float(amount),
                                               "usd", sid, owner=email))
            acct["card"] = {"brand": "Visa", "last4": "4242", "exp": f"12/{today.year + 3}"}
            n = len(acct["invoices"]) + 1
            p = plan(plan_id)
            what = p["name"] + (f" ({platform.title()})" if p.get("platforms") == "one" and platform else "")
            acct["invoices"].append({"number": f"OMX-{today.year}-{n:04d}", "date": today,
                                     "description": f"Omnix {what}, {'annual' if interval == 'year' else 'monthly'}"
                                                    f" subscription, {today:%b %d, %Y} – {end:%b %d, %Y}",
                                     "amount": float(amount), "currency": "USD", "status": "paid",
                                     "pdf": None, "url": None, "name": name, "email": email})
        return None

    def portal(self, email, return_url):
        return None

    def set_cancel(self, email, sub_id, cancel):
        with _store_lock:
            for s in self._acct(email)["subs"]:
                if s["id"] == sub_id:
                    s["cancel_at_period_end"] = bool(cancel)
                    return
        raise PermissionError("This subscription belongs to another account.")

    def set_members(self, email, sub_id, members):
        with _store_lock:
            for s in self._acct(email)["subs"]:
                if s["id"] == sub_id:
                    s["members"] = list(members)
                    return
        raise PermissionError("This subscription belongs to another account.")

    def all_subscriptions(self):
        with _store_lock:
            return [dict(s, owner=e) for e, a in _store.items() for s in a["subs"]]


# ---------------------------------------------------------------------------
# entry points
# ---------------------------------------------------------------------------
def backend():
    try:
        cfg = st.secrets.get("stripe", {}) or {}
        omnix = st.secrets.get("omnix", {}) or {}
    except Exception:
        cfg, omnix = {}, {}
    if cfg.get("secret_key"):
        return StripeBilling(dict(cfg))
    if omnix.get("billing_preview"):
        return PreviewBilling()
    return None


def forget_cache():
    for k in [k for k in st.session_state.keys() if str(k).startswith("omnix_bill_")]:
        del st.session_state[k]


def snapshot(email, refresh=False):
    """Everything the Account page and the access check need, cached per session for two minutes."""
    be = backend()
    if be is None or not email:
        return None
    key = "omnix_bill_" + hashlib.sha256(email.encode()).hexdigest()[:12]
    hit = st.session_state.get(key)
    if hit and not refresh and time.time() - hit["at"] < CACHE_SECONDS:
        return hit
    try:
        data = {"at": time.time(), "subs": be.subscriptions(email), "teams": be.memberships(email), "error": ""}
    except Exception as e:                                   # Stripe unreachable: fail closed, say why
        data = {"at": time.time(), "subs": [], "teams": [], "error": str(e)[:200]}
    st.session_state[key] = data
    return data


def entitlement(email):
    """{'plan', 'platforms', 'expires'} from the person's live subscriptions and Lab memberships, or None."""
    snap = snapshot(email)
    if not snap:
        return None
    live = [s for s in snap["subs"] + snap["teams"] if s["status"] in LIVE]
    if not live:
        return None
    platforms = sorted({p for s in live for p in s["platforms"]}, key=ALL.index)
    best = max(live, key=lambda s: len(s["platforms"]))
    own = [s for s in live if s in snap["subs"]]
    name = best["plan_name"] + ("" if best in own else " (team member)")
    ending = [s["renews_on"] for s in live if s["renews_on"] and s["cancel_at_period_end"]]
    return {"plan": name, "platforms": platforms,
            "expires": str(max(ending)) if ending and len(ending) == len(live) else ""}
