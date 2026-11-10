"""
Contact form (Contact page, "Get a quote", plan buttons).

The visitor picks a reason, edits the pre-filled subject and writes a message. Two ways to deliver it:

1. Email sending from the server, when an SMTP account is configured in the secrets:

       [omnix.smtp]
       host = "smtp.gmail.com"
       port = 587
       username = "you@gmail.com"
       password = "app password"          # Gmail: an App Password, not your normal password
       to = "you@gmail.com"               # where enquiries arrive (defaults to the Contact email)

   The message is sent to you with Reply-To set to the visitor, so you simply press Reply.

2. Otherwise, an "Open in your email app" button opens the visitor's own email program with the
   recipient, subject and message already filled in, plus a copy option.
"""

import re
import smtplib
import ssl
import time
import urllib.parse
from email.message import EmailMessage

import streamlit as st

from . import config

REASONS = {
    "Request a quote": "Quote request: Omnix subscription",
    "Subscription and billing": "Subscription and billing question",
    "Purchase order or invoice": "Purchase order / invoice request",
    "Free trial or demo account": "Trial / demo account request",
    "Technical support": "Technical support request",
    "Question about an analysis": "Question about an Omnix analysis",
    "Partnership or collaboration": "Partnership / collaboration enquiry",
    "Other": "Omnix enquiry",
}
EMAIL_RE = re.compile(r"^[^@\s]+@[^@\s]+\.[^@\s]+$")
_SENT_KEY = "omnix_contact_sent_at"
MAX_PER_WINDOW, WINDOW_S = 3, 600          # at most 3 messages per 10 minutes per visitor session


def _smtp():
    try:
        cfg = dict((st.secrets.get("omnix", {}) or {}).get("smtp", {}) or {})
    except Exception:
        return None
    if any("PASTE" in str(cfg.get(k, "")) or "CHANGE-ME" in str(cfg.get(k, "")) for k in ("username", "password")):
        return None                      # template placeholder not filled in yet: use the email-app fallback
    return cfg if cfg.get("host") and cfg.get("username") and cfg.get("password") else None


def set_reason(reason: str, plan: str = ""):
    """Callback for 'Get a quote' / plan buttons: open the Contact page with the form pre-filled."""
    st.session_state["omnix_contact_reason"] = reason
    st.session_state["omnix_contact_subject"] = REASONS.get(reason, "Omnix enquiry") + (f": {plan} plan" if plan else "")
    st.session_state.omnix_page = "contact"


def _on_reason_change():
    st.session_state["omnix_contact_subject"] = REASONS.get(st.session_state.get("omnix_contact_reason"),
                                                            "Omnix enquiry")


def _body(name, email, org, reason, message):
    return (f"{message.strip()}\n\n---\nName: {name}\nEmail: {email}\n"
            + (f"Organization: {org}\n" if org else "") + f"Reason: {reason}\nSent from the Omnix website")


def _send(to_addr, name, email, org, reason, subject, message):
    cfg = _smtp()
    msg = EmailMessage()
    msg["Subject"] = f"[Omnix] {subject}"
    msg["From"] = cfg.get("from", cfg["username"])
    msg["To"] = cfg.get("to", to_addr)
    msg["Reply-To"] = f"{name} <{email}>"
    msg.set_content(_body(name, email, org, reason, message))
    port = int(cfg.get("port", 587))
    if port == 465:
        with smtplib.SMTP_SSL(cfg["host"], port, context=ssl.create_default_context(), timeout=20) as s:
            s.login(cfg["username"], cfg["password"])
            s.send_message(msg)
    else:
        with smtplib.SMTP(cfg["host"], port, timeout=20) as s:
            s.starttls(context=ssl.create_default_context())
            s.login(cfg["username"], cfg["password"])
            s.send_message(msg)


def form():
    """The contact form. Returns nothing; shows success / the email-app fallback itself."""
    to_addr = config.CONTACT.get("email", "")
    if "omnix_contact_reason" not in st.session_state:
        st.session_state["omnix_contact_reason"] = "Request a quote"
    if "omnix_contact_subject" not in st.session_state:
        _on_reason_change()
    st.markdown("### Send us a message")
    st.selectbox("Reason for contacting", list(REASONS), key="omnix_contact_reason", on_change=_on_reason_change)
    with st.form("omnix_contact_form", border=True):
        c1, c2 = st.columns(2)
        name = c1.text_input("Your name *", key="omnix_contact_name")
        email = c2.text_input("Your email *", key="omnix_contact_email", placeholder="name@university.edu")
        org = st.text_input("Organization (optional)", key="omnix_contact_org")
        subject = st.text_input("Subject *", key="omnix_contact_subject")
        message = st.text_area("Message *", key="omnix_contact_message", height=160,
                               placeholder="For a quote: the plan, number of users, platform(s) and whether you "
                                           "need a purchase order.")
        sent = st.form_submit_button("Send message")
    if not sent:
        return
    reason = st.session_state["omnix_contact_reason"]
    problems = []
    if not name.strip():
        problems.append("your name")
    if not EMAIL_RE.match(email.strip()):
        problems.append("a valid email address")
    if not subject.strip():
        problems.append("a subject")
    if len(message.strip()) < 5:
        problems.append("a message")
    if problems:
        st.error("Please add " + ", ".join(problems) + ".")
        return
    if _smtp():
        now = time.time()
        recent = [t for t in st.session_state.get(_SENT_KEY, []) if now - t < WINDOW_S]
        if len(recent) >= MAX_PER_WINDOW:
            st.error("You have sent several messages in a short time. Please wait a few minutes and try again.")
            return
        try:
            _send(to_addr, name.strip(), email.strip(), org.strip(), reason, subject.strip(), message)
            st.session_state[_SENT_KEY] = recent + [now]
            st.success("Thank you. Your message has been sent, and we will reply to "
                       f"{email.strip()} as soon as possible.")
            return
        except Exception:
            st.warning("The message could not be sent from the website just now. Please use the button below "
                       "to send it from your own email app.")
    body = _body(name.strip(), email.strip(), org.strip(), reason, message)
    link = (f"mailto:{to_addr}?subject={urllib.parse.quote('[Omnix] ' + subject.strip())}"
            f"&body={urllib.parse.quote(body)}")
    st.success("Your message is ready. Click below to open it in your email app and press Send.")
    st.link_button("Open in your email app", link, type="primary")
    with st.expander("Or copy the message"):
        st.code(f"To: {to_addr}\nSubject: [Omnix] {subject.strip()}\n\n{body}", language=None, wrap_lines=True)
