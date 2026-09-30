"""
Subscription access for Omnix.

Without a subscription a visitor can run every analysis on the built-in demo datasets; uploading
their own files (data, metadata, annotations, gene sets) requires an active subscription.

Subscriptions are granted in the app's secrets (Streamlit Cloud: App settings -> Secrets;
locally: .streamlit/secrets.toml). Two ways, usable together:

1. Access keys. Generate one per subscriber with `python tools/make_access_key.py`, give the key to
   the subscriber, and paste the printed line under [omnix.access_keys]. Only a SHA-256 hash of the
   key is stored, never the key itself. The subscriber enters the key on the Subscription page.

       [omnix.access_keys]
       "3f9a...c1" = { name = "Jane Doe", plan = "Multi-Omics", expires = "2027-09-30" }
       "8b21...7e" = { name = "Sam Lee", plan = "Single-Omics", platforms = ["proteomics"], expires = "2027-03-31" }

   `platforms` limits a key to some platforms (metabolomics, proteomics, transcriptomics); without
   it the key covers all three.

2. Accounts (Account page). People sign in with Google or Microsoft; their Stripe subscriptions
   (omnix_portal/billing.py) decide what they can upload, and so do Lab team memberships. Emails
   listed here are subscribers without paying online (e.g. invoiced by purchase order):

       [omnix.subscriber_emails]
       "jane.doe@university.edu" = { plan = "Multi-Omics", expires = "2027-09-30" }

`open_access = true` under [omnix] gives everyone full access (internal use, testing).
An `expires` date is optional; access ends after that day.

Owner: a key with `role = "admin"`, or signing in with an email listed in `owner_emails = [...]`
under [omnix], also opens the Owner console (subscriber list, key
creation, lock-out monitor). Make yours with `python tools/make_access_key.py --owner --name "..."`.

Guessing protection: after MAX_FAILS wrong keys within LOCK_MINUTES, activation is paused for that
visitor (by IP address when available, and by browser session) for LOCK_MINUTES.
"""

import datetime as _dt
import hmac
import threading
import time

import streamlit as st

from . import keys as _keys

_SESSION_KEY = "omnix_access"          # portal-owned session key (kept when moving between apps)


# ---------------------------------------------------------------------------
# configuration
# ---------------------------------------------------------------------------
def _settings() -> dict:
    try:
        return dict(st.secrets.get("omnix", {}))
    except Exception:              # no secrets file at all
        return {}


hash_key = _keys.hash_key


def _active(info: dict) -> bool:
    exp = str(info.get("expires", "") or "").strip()
    if not exp:
        return True
    try:
        return _dt.date.today() <= _dt.date.fromisoformat(exp)
    except ValueError:
        return False                # malformed date: treat as expired rather than open-ended


ALL_PLATFORMS = ("metabolomics", "proteomics", "transcriptomics")


def _platforms(info: dict) -> list:
    """Platforms a subscription covers; all of them when `platforms` is missing, empty or "all"."""
    raw = info.get("platforms")
    if raw is None or raw == "" or raw == "all":
        return list(ALL_PLATFORMS)
    if isinstance(raw, str):
        raw = raw.replace(";", ",").split(",")
    wanted = {str(x).strip().lower() for x in raw}
    if "all" in wanted:
        return list(ALL_PLATFORMS)
    return [p for p in ALL_PLATFORMS if p in wanted]


def sign_in_available() -> bool:
    """True when Streamlit authentication is configured (st.login can be used)."""
    try:
        return hasattr(st, "login") and "auth" in st.secrets
    except Exception:
        return False


# ---------------------------------------------------------------------------
# status
# ---------------------------------------------------------------------------
def status() -> dict:
    """{'subscribed': bool, 'name', 'plan', 'expires', 'via', 'platforms'} for the current visitor."""
    cfg = _settings()
    if cfg.get("open_access"):
        return {"subscribed": True, "name": "", "plan": "Full access", "expires": "", "via": "open",
                "platforms": list(ALL_PLATFORMS)}
    held = st.session_state.get(_SESSION_KEY)
    if held and _active(held):
        return {"platforms": list(ALL_PLATFORMS), **held, "subscribed": True}
    from . import account, billing                           # late import: both read st.secrets
    user = account.current_user()
    if user:
        email = user["email"]
        owners = {str(e).strip().lower() for e in (cfg.get("owner_emails") or [])}
        role = "admin" if email in owners else ""
        if role:
            return {"subscribed": True, "name": email, "plan": "Owner", "expires": "", "via": "account",
                    "platforms": list(ALL_PLATFORMS), "role": role}
        emails = {str(k).strip().lower(): dict(v) if hasattr(v, "keys") else {} for k, v in
                  dict(cfg.get("subscriber_emails", {})).items()}
        info = emails.get(email)
        if info is not None and _active(info):
            return {"subscribed": True, "name": email, "plan": info.get("plan", "Subscription"),
                    "expires": str(info.get("expires", "") or ""), "via": "account", "platforms": _platforms(info),
                    "role": str(info.get("role", "") or "")}
        ent = billing.entitlement(email)
        if ent:
            return {"subscribed": True, "name": email, "plan": ent["plan"], "expires": ent["expires"],
                    "via": "account", "platforms": ent["platforms"], "role": ""}
        return {"subscribed": False, "name": email, "plan": "", "expires": "", "via": "account", "platforms": []}
    return {"subscribed": False, "name": "", "plan": "", "expires": "", "via": "", "platforms": []}


def is_subscriber(platform: str = None) -> bool:
    """True with an active subscription; with `platform`, only if the plan covers that platform."""
    s = status()
    return s["subscribed"] and (platform is None or platform in s["platforms"])


def is_admin() -> bool:
    """Owner: an active key or email with role = "admin" (open_access does not make anyone owner)."""
    s = status()
    return s["subscribed"] and s.get("via") in ("key", "account") and s.get("role") == "admin"


# ---------------------------------------------------------------------------
# guessing protection
# ---------------------------------------------------------------------------
MAX_FAILS = 5
LOCK_MINUTES = 15
_fails = {}                       # client id -> [failure times]; shared by all sessions in this process
_fails_lock = threading.Lock()


def _client_ids():
    ids = []
    try:
        from streamlit.runtime.scriptrunner import get_script_run_ctx
        ctx = get_script_run_ctx()
        if ctx is not None:
            ids.append("s:" + ctx.session_id)
    except Exception:
        pass
    ip = None
    try:
        fwd = st.context.headers.get("X-Forwarded-For") if st.context.headers else None
        ip = (fwd.split(",")[0].strip() if fwd else None) or st.context.ip_address
    except Exception:
        pass
    if isinstance(ip, str) and ip:
        ids.append("ip:" + ip)
    return ids or ["s:unknown"]


def _recent(times, now):
    return [t for t in times if now - t < LOCK_MINUTES * 60]


def locked_minutes() -> int:
    """Minutes left before this visitor may try again (0 = not locked)."""
    now = time.time()
    with _fails_lock:
        worst = 0.0
        for cid in _client_ids():
            times = _recent(_fails.get(cid, []), now)
            if len(times) >= MAX_FAILS:
                worst = max(worst, times[-MAX_FAILS] + LOCK_MINUTES * 60 - now)
    return int(worst // 60) + 1 if worst > 0 else 0


def _record_failure():
    now = time.time()
    with _fails_lock:
        for cid in _client_ids():
            _fails[cid] = _recent(_fails.get(cid, []), now) + [now]
        if len(_fails) > 5000:                                  # keep memory bounded
            for cid in [c for c, t in _fails.items() if not _recent(t, now)]:
                _fails.pop(cid, None)


def _clear_failures():
    with _fails_lock:
        for cid in _client_ids():
            _fails.pop(cid, None)


def failure_report():
    """[(client, failures in the window, locked)] for the Owner console; IPs only, sessions hidden."""
    now = time.time()
    with _fails_lock:
        rows = [(cid[3:], len(_recent(t, now)), len(_recent(t, now)) >= MAX_FAILS)
                for cid, t in _fails.items() if cid.startswith("ip:") and _recent(t, now)]
        sessions = sum(1 for cid, t in _fails.items() if cid.startswith("s:") and _recent(t, now))
    return sorted(rows, key=lambda r: -r[1]), sessions


# ---------------------------------------------------------------------------
# keys
# ---------------------------------------------------------------------------
def subscribers():
    """Every entry in [omnix.access_keys] (for the Owner console)."""
    rows = []
    today = _dt.date.today()
    for stored, info in dict(_settings().get("access_keys", {})).items():
        info = dict(info) if hasattr(info, "keys") else {}
        exp = str(info.get("expires", "") or "").strip()
        days = None
        if exp:
            try:
                days = (_dt.date.fromisoformat(exp) - today).days
            except ValueError:
                days = -1
        rows.append({"name": info.get("name", ""), "plan": info.get("plan", ""), "platforms": _platforms(info),
                     "expires": exp, "days_left": days, "active": _active(info),
                     "role": str(info.get("role", "") or ""), "key_id": str(stored)[:8]})
    return rows


def email_subscribers():
    """Every entry in [omnix.subscriber_emails] (people who sign in with Google/Microsoft and get access
    without paying online, e.g. invoiced institutions or trial users)."""
    rows = []
    today = _dt.date.today()
    for email, info in dict(_settings().get("subscriber_emails", {})).items():
        info = dict(info) if hasattr(info, "keys") else {}
        exp = str(info.get("expires", "") or "").strip()
        days = None
        if exp:
            try:
                days = (_dt.date.fromisoformat(exp) - today).days
            except ValueError:
                days = -1
        rows.append({"email": str(email).strip().lower(), "plan": info.get("plan", ""),
                     "platforms": _platforms(info), "expires": exp, "days_left": days, "active": _active(info)})
    return rows


def find_key(key: str):
    """Key ID (first 8 characters of its hash) if this key is in the secrets, else None."""
    digest = hash_key(key or "")
    for stored in dict(_settings().get("access_keys", {})):
        if hmac.compare_digest(str(stored).lower(), digest):
            return str(stored)[:8]
    return None


def redeem_key(key: str):
    """Check an access key; on success remember it for this session. Returns (ok, message)."""
    key = (key or "").strip()
    if not key:
        return False, "Enter your access key."
    wait = locked_minutes()
    if wait:
        return False, (f"Too many incorrect access keys. For your security, activation is paused for {wait} "
                       f"minute{'s' if wait != 1 else ''}. Please try again later or contact us.")
    keys = dict(_settings().get("access_keys", {}))
    digest = hash_key(key)
    for stored, info in keys.items():
        if hmac.compare_digest(str(stored).lower(), digest):
            info = dict(info) if hasattr(info, "keys") else {}
            if not _active(info):
                return False, f"This access key expired on {info.get('expires')}. Please renew your subscription."
            st.session_state[_SESSION_KEY] = {"name": info.get("name", ""), "plan": info.get("plan", "Subscription"),
                                              "expires": str(info.get("expires", "") or ""), "via": "key",
                                              "platforms": _platforms(info),
                                              "role": str(info.get("role", "") or "")}
            _clear_failures()
            return True, "Access key accepted. Your subscription is active."
    _record_failure()
    left = MAX_FAILS - max((len(_recent(_fails.get(c, []), time.time())) for c in _client_ids()), default=0)
    extra = f" {left} attempt{'s' if left != 1 else ''} left before a {LOCK_MINUTES}-minute pause." if 0 < left <= 2 else ""
    return False, "This access key was not recognized. Check it and try again, or contact us." + extra


def forget_key():
    st.session_state.pop(_SESSION_KEY, None)


# ---------------------------------------------------------------------------
# demo-only mode: every file upload in every platform needs a subscription
# ---------------------------------------------------------------------------
LOCKED_HELP = "Uploading your own files requires an Omnix subscription. Without one, use the built-in demo datasets."
_installed = False


def install_upload_gate():
    """Wrap Streamlit's file uploader once per process. The wrapper reads the CURRENT session's
    state on every call, so subscribed and demo visitors can use the site at the same time."""
    global _installed
    if _installed:
        return
    import streamlit as _st
    from streamlit.elements.widgets.file_uploader import FileUploaderMixin

    original = FileUploaderMixin.file_uploader

    def gated(self, label, *args, **kwargs):
        if _st.session_state.get("omnix_demo_only", False):
            kwargs["disabled"] = True
            kwargs["help"] = LOCKED_HELP
            original(self, f"{label}  ·  subscription required", *args, **kwargs)
            return None
        return original(self, label, *args, **kwargs)

    gated.__doc__ = original.__doc__
    FileUploaderMixin.file_uploader = gated
    main_dg = getattr(_st, "_main", None)
    if main_dg is not None:
        _st.file_uploader = lambda label, *a, **k: gated(main_dg, label, *a, **k)
    _installed = True


def set_demo_only(flag: bool):
    st.session_state["omnix_demo_only"] = bool(flag)
