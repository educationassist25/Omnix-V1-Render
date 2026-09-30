"""
Who is signed in.

Sign-in uses Streamlit's built-in OpenID Connect login (st.login) with Google and/or Microsoft, so
Omnix never handles passwords. Configure the providers in the app's secrets ([auth], [auth.google],
[auth.microsoft]; see docs/ACCOUNT_AND_BILLING_SETUP.md). Signing in for the first time is the
registration: the Stripe customer record is created when the person first subscribes.

Preview: with `billing_preview = true` under [omnix] (and no Stripe key), a "Preview with a demo
account" button signs in a demo user so the Account pages can be tried without real sign-in or
payment. Never enable preview on the public site.
"""

import streamlit as st

PROVIDERS = {"google": "Google", "microsoft": "Microsoft"}
_PREVIEW_KEY = "omnix_preview_user"
PREVIEW_USER = {"email": "demo.researcher@example.edu", "name": "Demo Researcher", "via": "preview"}


def _secrets(section):
    try:
        return st.secrets.get(section, {}) or {}
    except Exception:
        return {}


def preview_enabled() -> bool:
    return bool(_secrets("omnix").get("billing_preview")) and not _secrets("stripe").get("secret_key")


def auth_packages_ok() -> bool:
    """Sign-in needs Authlib and httpx (requirements: streamlit[auth]). Without them the login
    route fails with a bare "Internal server error", so buttons are hidden instead."""
    try:
        from authlib.integrations import starlette_client  # noqa: F401
        import httpx  # noqa: F401
        return True
    except Exception:
        return False


def providers() -> dict:
    """Configured sign-in providers {id: label}."""
    auth = _secrets("auth")
    if not hasattr(st, "login") or not auth or not auth_packages_ok():
        return {}
    found = {k: v for k, v in PROVIDERS.items() if k in auth}
    if not found and auth.get("client_id"):          # a single unnamed provider
        found = {"": "your account"}
    return found


def current_user():
    """{'email', 'name', 'via'} for the signed-in person, or None."""
    try:
        u = st.user
        if u.get("is_logged_in"):
            email = str(u.get("email", "") or "").strip().lower()
            if email and u.get("email_verified", True) is not False:
                name = str(u.get("name", "") or "") or email.split("@")[0]
                return {"email": email, "name": name, "via": "sign-in"}
    except Exception:
        pass
    if preview_enabled() and st.session_state.get(_PREVIEW_KEY):
        return dict(PREVIEW_USER)
    return None


def config_problems(provider: str = "google") -> list:
    """Plain-language list of what is missing or wrong in the sign-in secrets (never shows values)."""
    out = []
    auth = _secrets("auth")
    if not auth:
        return ["There is no [auth] section in Secrets."]
    uri = str(auth.get("redirect_uri", "") or "")
    if not uri:
        out.append("[auth] is missing redirect_uri.")
    elif not uri.startswith("https://") and "localhost" not in uri:
        out.append("redirect_uri must start with https://")
    elif not uri.rstrip("/").endswith("/oauth2callback"):
        out.append("redirect_uri must end with /oauth2callback")
    elif "CHANGE-ME" in uri or "YOUR-APP" in uri or "YOUR-SITE" in uri:
        out.append("redirect_uri still contains the placeholder app address.")
    if not auth.get("cookie_secret"):
        out.append("[auth] is missing cookie_secret.")
    sec = auth.get(provider)
    if sec is None:
        out.append(f"There is no [auth.{provider}] section.")
    else:
        for k in ("client_id", "client_secret", "server_metadata_url"):
            v = str(sec.get(k, "") or "")
            if not v:
                out.append(f"[auth.{provider}] is missing {k}.")
            elif "CHANGE-ME" in v or "PASTE" in v:
                out.append(f"[auth.{provider}] {k} still contains a placeholder.")
            elif v != v.strip():
                out.append(f"[auth.{provider}] {k} has a space at the start or end.")
    return out


def sign_in(provider: str):
    if provider == "preview":
        st.session_state[_PREVIEW_KEY] = True
        return
    try:
        st.login(provider) if provider else st.login()
    except Exception:
        probs = config_problems(provider or "google")
        st.session_state["omnix_signin_error"] = probs or [
            "Sign-in settings could not be used. Check the [auth] sections in Secrets."]


def sign_out():
    user = current_user()
    for k in [k for k in st.session_state.keys() if str(k).startswith("omnix_bill_")]:
        del st.session_state[k]
    if user and user["via"] == "preview":
        st.session_state.pop(_PREVIEW_KEY, None)
        return
    st.logout()
