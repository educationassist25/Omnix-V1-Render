"""
Omnix portal: the header, menu and side panel shared by every page, the portal's own pages (Home,
Subscription, Contact, Help), subscription access, and the three analysis platforms run inside it.

Pages ("routes"): home | metabolomics | proteomics | transcriptomics | subscription | contact |
help-metabolomics | help-proteomics | help-transcriptomics. The current one is kept in
st.session_state["omnix_page"] and mirrored in the address bar (?page=...), so links can be shared.
"""

import base64
import datetime as dt
import functools
import html
import os
import re

import pandas as pd
import streamlit as st

from omnix_portal import access, account, account_page, billing, config, contact, demo_files, help_content, isolation, keys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TUTORIALS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tutorials")
APPS = list(config.PLATFORMS)
PAGES = ["home", *APPS, "subscription", "contact", "account", "owner", "data-preparation"] + [f"help-{k}" for k in APPS]

TEAL, TEAL_DARK, ORANGE, INK, INK2, LINE, GROUND = "#00695C", "#00574D", "#E65100", "#111816", "#4A5754", "#DCE3E0", "#F4F6F5"

ICONS = {
    "metabolomics": '<path d="M9 3h6"></path><path d="M10 3v6.5L4.5 19a1.5 1.5 0 0 0 1.3 2.2h12.4a1.5 1.5 0 0 0 1.3-2.2L14 9.5V3"></path><path d="M7.5 15h9"></path>',
    "proteomics": '<circle cx="6" cy="7" r="2.5"></circle><circle cx="18" cy="6" r="2.5"></circle><circle cx="12" cy="17" r="2.5"></circle><circle cx="19" cy="17" r="1.8"></circle><path d="M8.3 8.3l2.6 6.4"></path><path d="M15.6 7.1l-2.6 7.6"></path><path d="M8.5 7l7-1"></path><path d="M14.5 17h2.7"></path>',
    "transcriptomics": '<path d="M7 3c0 4 10 5 10 9s-10 5-10 9"></path><path d="M17 3c0 4-10 5-10 9s10 5 10 9"></path><path d="M8.5 6.5h7"></path><path d="M8.5 17.5h7"></path>',
    "workflow": '<path d="M4 6h10"></path><path d="M4 12h16"></path><path d="M4 18h7"></path><circle cx="17" cy="6" r="2"></circle><circle cx="14" cy="18" r="2"></circle>',
    "stats": '<path d="M4 20V10"></path><path d="M10 20V4"></path><path d="M16 20v-7"></path><path d="M3 20h18"></path>',
    "figure": '<rect x="3" y="4" width="18" height="14" rx="2"></rect><path d="M7 14l3-3 3 2 4-5"></path><path d="M8 21h8"></path>',
}


def _svg(paths, size=28, color=TEAL):
    return (f'<svg width="{size}" height="{size}" viewBox="0 0 24 24" fill="none" stroke="{color}" '
            f'stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">{paths}</svg>')


ASSETS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "assets")


@functools.lru_cache(maxsize=8)
def _img_uri(name):
    """An image in omnix_portal/assets as a data URI ('' if missing, so the site still works without it)."""
    path = os.path.join(ASSETS, name)
    if not os.path.isfile(path):
        return ""
    with open(path, "rb") as fh:
        return "data:image/png;base64," + base64.b64encode(fh.read()).decode()


def _logo(cls):
    uri = _img_uri("omnix_logo_web.png")
    return f'<img src="{uri}" class="omx-logo {cls}" alt="Omnix logo">' if uri else ""


def _slug(text):
    return re.sub(r"[^a-z0-9]+", "-", text.lower()).strip("-")


# ---------------------------------------------------------------------------
# routing
# ---------------------------------------------------------------------------
def go(page):
    """Button callback: open a page (runs before the next script run, so it renders at once)."""
    st.session_state.omnix_page = page


def _current_page():
    ss = st.session_state
    if "omnix_page" not in ss:
        wanted = st.query_params.get("page", "home")
        ss.omnix_page = wanted if wanted in PAGES else "home"
    page = ss.omnix_page if ss.omnix_page in PAGES else "home"
    if st.query_params.get("page") != page:
        st.query_params["page"] = page
    return page


# ---------------------------------------------------------------------------
# styles
# ---------------------------------------------------------------------------
def _styles(page, portal_page):
    active = "help" if page.startswith("help-") or page == "data-preparation" else page
    side_active = "help" if page.startswith("help-") else page
    portal_only = ""
    if portal_page:
        # portal pages use the Omnix typefaces throughout; inside a platform only the header does,
        # so every platform looks exactly as it does on its own
        portal_only = f"""
        .stMainBlockContainer, .stMainBlockContainer p, .stMainBlockContainer li, .stMainBlockContainer td {{ font-family: 'IBM Plex Sans', system-ui, sans-serif; }}
        .stMainBlockContainer h1, .stMainBlockContainer h2, .stMainBlockContainer h3, .stMainBlockContainer h4 {{ font-family: 'Sora', system-ui, sans-serif; color: {INK}; }}
        [data-testid="stAppViewContainer"] {{ background: {GROUND}; }}
        /* responsive: side-by-side cards wrap onto new rows instead of squeezing, and words never break
           in the middle */
        .stMainBlockContainer [data-testid="stHorizontalBlock"] {{ flex-wrap: wrap !important; row-gap: 16px; }}
        .stMainBlockContainer [data-testid="stHorizontalBlock"] > [data-testid="stColumn"] {{ min-width: min(220px, 100%) !important; }}
        [class*="st-key-omnix_footer"] [data-testid="stHorizontalBlock"] > [data-testid="stColumn"] {{ min-width: min(150px, 100%) !important; }}
        .stMainBlockContainer p, .stMainBlockContainer li, .stMainBlockContainer h1, .stMainBlockContainer h2,
        .stMainBlockContainer h3, .stMainBlockContainer h4, .stMainBlockContainer button p {{
          word-break: normal !important; overflow-wrap: break-word !important; hyphens: manual; }}
        .stMainBlockContainer h2 {{ font-size: clamp(26px, 3.2vw, 40px) !important; }}
        .omx-card h2 {{ font-size: clamp(20px, 2vw, 26px) !important; }}
        @media (max-width: 1250px) {{
          .omx-plan ul {{ min-height: 0 !important; }}
          .omx-feature-card {{ min-height: 0 !important; }}
          .omx-tile h4 {{ white-space: normal !important; }}
          .omx-hero .t {{ font-size: 30px !important; letter-spacing: 7px !important; }}
          .omx-about .hl {{ font-size: 22px !important; }}
          .omx-about p, .omx-lead {{ font-size: 16.5px !important; }}
          [class*="st-key-omnix_footer"] [data-testid="stHorizontalBlock"] > [data-testid="stColumn"]:first-child {{
            min-width: 100% !important; margin-bottom: 8px; }}
        }}
        @media (max-width: 700px) {{
          .omx-hero {{ padding: 34px 18px !important; }}
          .omx-hero .t {{ font-size: 22px !important; letter-spacing: 4px !important; }}
          [class*="st-key-omnix_footer"] {{ padding: 26px 20px 22px 20px !important; }}
        }}
        /* side panel on portal pages: same look as the platforms' own side panel */
        [data-testid="stSidebar"] {{ background: linear-gradient(180deg, #f8fafc 0%, #e6f4f2 100%) !important;
          border-right: 1px solid #e2e8f0; }}
        div[data-testid="stVerticalBlock"][class*="st-key-omnix_side_nav"] {{ background: #e6f4f2; border: 1px solid #bfe3dc;
          border-radius: 10px; padding: 8px !important; gap: 4px !important; }}
        [class*="st-key-omnix_side_nav"] [data-testid="stButton"] button {{ background: transparent !important;
          border: 1px solid transparent !important; border-left: 3px solid #bfe3dc !important; box-shadow: none !important;
          border-radius: 6px !important; padding: 9px 12px !important; width: 100%; }}
        [class*="st-key-omnix_side_nav"] [data-testid="stButton"] button p {{ color: #00695c !important; font-weight: 700 !important;
          font-size: 17px !important; white-space: nowrap !important; }}
        [class*="st-key-omnix_side_nav"] [data-testid="stButton"] button:hover {{ background: rgba(0, 105, 92, 0.10) !important;
          border-left-color: #00897b !important; }}
        [class*="st-key-omnix_side_nav"] [class*="st-key-omnix_sidebtn_{_slug(side_active)}"] [data-testid="stButton"] button {{
          background: #00695c !important; border-left-color: {ORANGE} !important; }}
        [class*="st-key-omnix_side_nav"] [class*="st-key-omnix_sidebtn_{_slug(side_active)}"] [data-testid="stButton"] button p {{
          color: #ffffff !important; }}
        """
    st.markdown(
        f"""
<style>
@import url('https://fonts.googleapis.com/css2?family=Sora:wght@500;600;700;800&family=IBM+Plex+Sans:wght@400;500;600&display=swap');
{portal_only}
/* the side panel is always shown: no collapse control on desktop (phones keep Streamlit's toggle) */
@media (min-width: 1100px) {{
  [data-testid="stSidebarCollapseButton"], [data-testid="stSidebarCollapsedControl"],
  [data-testid="stExpandSidebarButton"] {{ display: none !important; }}
}}
/* ---- responsive: header scales down on smaller windows ---- */
@media (max-width: 1250px) {{
  .stApp [class*="st-key-omnix_brand"]:not([class*="compact"]) [data-testid="stButton"] button p {{ font-size: 50px !important; }}
  .stApp [data-testid="stMarkdownContainer"] div.omx-tagline {{ font-size: 18px; }}
  .stApp [class*="st-key-omnix_header"] [class*="st-key-omnix_nav"] {{ column-gap: 24px !important; row-gap: 4px !important; }}
  .stApp [class*="st-key-omnix_nav"] [data-testid="stButton"] button p,
  .stApp [class*="st-key-omnix_nav"] [data-testid="stPopover"] button p {{ font-size: 17px !important; }}
}}
@media (max-width: 900px) {{
  .stApp [class*="st-key-omnix_brand"]:not([class*="compact"]) [data-testid="stButton"] button p {{ font-size: 40px !important; }}
  .stApp [class*="st-key-omnix_brand_compact"] [data-testid="stButton"] button p {{ font-size: 24px !important; }}
  .stApp [data-testid="stMarkdownContainer"] div.omx-tagline {{ font-size: 16px; }}
  .stApp [class*="st-key-omnix_header"] [class*="st-key-omnix_nav"] {{ column-gap: 18px !important; }}
  .stApp [class*="st-key-omnix_nav"] [data-testid="stButton"] button p,
  .stApp [class*="st-key-omnix_nav"] [data-testid="stPopover"] button p {{ font-size: 15.5px !important; }}
}}
@media (max-width: 640px) {{
  .stApp .stMainBlockContainer {{ padding-top: 2.2rem !important; padding-left: 1rem !important; padding-right: 1rem !important; }}
  .stApp [class*="st-key-omnix_header"] {{ padding: 6px 2px 0 2px !important; }}
  .stApp [class*="st-key-omnix_brand"]:not([class*="compact"]) [data-testid="stButton"] button p {{ font-size: 34px !important; }}
  .stApp [data-testid="stMarkdownContainer"] div.omx-tagline {{ font-size: 14px !important; margin: 2px 0 2px 0; }}
  .stApp [class*="st-key-omnix_header"] [class*="st-key-omnix_nav"] {{ column-gap: 14px !important; row-gap: 0 !important; }}
  .stApp [class*="st-key-omnix_nav"] [data-testid="stButton"] button,
  .stApp [class*="st-key-omnix_nav"] [data-testid="stPopover"] button {{ padding: 2px 0 6px 0 !important; }}
  .stApp [class*="st-key-omnix_nav"] [data-testid="stButton"] button p,
  .stApp [class*="st-key-omnix_nav"] [data-testid="stPopover"] button p {{ font-size: 14.5px !important; }}
  .stApp [class*="st-key-omnix_nav_account"] [data-testid="stButton"] button {{ padding: 2px 12px 2px 10px !important; }}
}}


/* logo */
.omx-logo.full {{ height: 92px; width: auto; display: block; }}
.omx-logo.compact {{ height: 44px; width: auto; display: block; }}
.omx-side-head {{ display: flex; align-items: center; gap: 10px; }}
.omx-logo.side {{ height: 46px; width: auto; }}
.omx-fhead-row {{ display: flex; align-items: center; gap: 14px; }}
.omx-fbadge {{ background: #FFFFFF; border-radius: 12px; padding: 6px 8px; display: inline-flex; }}
.omx-logo.foot {{ height: 44px; width: auto; }}
@media (max-width: 900px) {{ .omx-logo.full {{ height: 68px; }} }}
@media (max-width: 640px) {{ .omx-logo.full {{ height: 50px; }} .omx-logo.compact {{ height: 34px; }} }}

/* Omnix wordmark: a button that returns to the home page */
[class*="st-key-omnix_brand"] [data-testid="stButton"] button {{
  background: transparent !important; border: none !important; box-shadow: none !important; transform: none !important;
  padding: 0 !important; min-height: 0 !important; height: auto !important; }}
[class*="st-key-omnix_brand"] [data-testid="stButton"] button p {{
  font-family: 'Sora', system-ui, sans-serif !important; font-weight: 800 !important; font-size: 64px !important;
  line-height: 1.05 !important; letter-spacing: -1.5px; color: {INK} !important; }}
[class*="st-key-omnix_brand"] [data-testid="stButton"] button:hover p {{ color: {TEAL_DARK} !important; }}
[class*="st-key-omnix_brand_compact"] [data-testid="stButton"] button p {{ font-size: 30px !important; letter-spacing: -0.6px; }}
[data-testid="stMarkdownContainer"] div.omx-tagline {{ font-family: 'IBM Plex Sans', system-ui, sans-serif; font-size: 21px;
  color: #3E4B48; margin: 4px 0 6px 0; }}
[data-testid="stMarkdownContainer"] div.omx-tagline.compact {{ font-size: 13.5px; margin: 1px 0 0 1px; line-height: 1.2; }}
/* logo sits right next to the wordmark */
.stApp [class*="st-key-omnix_brandrow"] {{ column-gap: 12px !important; }}
.stApp [class*="st-key-omnix_brandrow_compact"] {{ column-gap: 10px !important; }}
[class*="st-key-omnix_brandrow_compact"] [data-testid="stMarkdownContainer"] p {{ margin: 0 !important; }}

/* header nav: text-style menu items (override the platforms' button styling) */
[class*="st-key-omnix_nav"] [data-testid="stButton"] button,
[class*="st-key-omnix_nav"] [data-testid="stPopover"] button {{
  background: transparent !important; border: none !important; border-bottom: 3px solid transparent !important;
  border-radius: 0 !important; box-shadow: none !important; transform: none !important;
  padding: 4px 2px 12px 2px !important; min-height: 0 !important; height: auto !important; }}
[class*="st-key-omnix_nav"] [data-testid="stButton"] button p,
[class*="st-key-omnix_nav"] [data-testid="stPopover"] button p {{
  font-family: 'Sora', system-ui, sans-serif !important; font-size: 19px !important; font-weight: 500 !important;
  color: #17211F !important; white-space: nowrap !important; }}
[class*="st-key-omnix_nav"] [data-testid="stButton"] button:hover,
[class*="st-key-omnix_nav"] [data-testid="stPopover"] button:hover {{
  background: transparent !important; box-shadow: none !important; transform: none !important;
  border-bottom-color: #F5B48A !important; }}
[class*="st-key-omnix_nav"] [data-testid="stButton"] button:hover p,
[class*="st-key-omnix_nav"] [data-testid="stPopover"] button:hover p {{ color: {ORANGE} !important; }}
[class*="st-key-omnix_nav_{active}"] [data-testid="stButton"] button,
[class*="st-key-omnix_navpop_{active}"] [data-testid="stPopover"] button {{ border-bottom-color: {ORANGE} !important; }}
[class*="st-key-omnix_nav_{active}"] [data-testid="stButton"] button p,
[class*="st-key-omnix_navpop_{active}"] [data-testid="stPopover"] button p {{ color: {TEAL_DARK} !important; font-weight: 700 !important; }}
[class*="st-key-omnix_nav"] {{ column-gap: 44px !important; }}
[class*="st-key-omnix_navpop_"] [data-testid="stPopover"] button {{ margin-top: 3px !important; }}
[class*="st-key-omnix_header_compact"] [class*="st-key-omnix_nav"] button p {{ font-size: 16px !important; }}
[class*="st-key-omnix_header_compact"] [class*="st-key-omnix_nav"] {{ column-gap: 30px !important; }}

/* header band */
[class*="st-key-omnix_header"] {{ background: #FFFFFF; border-bottom: 1px solid {LINE}; padding: 8px 4px 0 4px; margin-bottom: 14px; }}

/* popover menus */
[data-testid="stPopoverBody"] {{ min-width: 340px; }}
[data-testid="stPopoverBody"] [data-testid="stButton"] button {{
  background: transparent !important; border: none !important; box-shadow: none !important; transform: none !important;
  justify-content: flex-start !important; text-align: left !important; padding: 8px 6px !important; }}
[data-testid="stPopoverBody"] [data-testid="stButton"] button p {{
  color: {TEAL_DARK} !important; font-family: 'Sora', system-ui, sans-serif !important; font-weight: 600 !important; font-size: 16px !important; }}
[data-testid="stPopoverBody"] [data-testid="stButton"] button:hover p {{ color: {ORANGE} !important; }}
.omx-menu-title {{ font-family: 'Sora', system-ui, sans-serif; font-weight: 700; font-size: 13px; letter-spacing: 0.6px;
  text-transform: uppercase; color: {INK2}; margin: 2px 0 6px 0; }}
.omx-kv {{ display: flex; gap: 12px; align-items: flex-start; margin: 10px 0; font-family: 'IBM Plex Sans', system-ui, sans-serif; }}
.omx-kv .k {{ font-size: 13px; color: {INK2}; }}
.omx-kv .v {{ font-size: 16px; color: #17211F; line-height: 1.4; }}

/* side panel: brand and access status */
.omx-side-brand {{ font-family: 'Sora', system-ui, sans-serif; font-weight: 800; font-size: 30px; color: {TEAL}; letter-spacing: -0.5px; }}
.omx-side-tag {{ font-size: 13px; color: #5B6B68; margin: 2px 0 6px 0; }}
.omx-access {{ border: 1px solid #bfe3dc; background: #FFFFFF; border-radius: 10px; padding: 12px 14px; margin-top: 6px;
  font-family: 'IBM Plex Sans', system-ui, sans-serif; }}
.omx-access .t {{ font-size: 12px; letter-spacing: 0.6px; text-transform: uppercase; color: {INK2}; font-weight: 700; }}
.omx-access .s {{ font-family: 'Sora', system-ui, sans-serif; font-size: 17px; font-weight: 700; margin-top: 4px; }}
.omx-access .s.on {{ color: {TEAL_DARK}; }} .omx-access .s.off {{ color: #8A4B08; }}
.omx-access .d {{ font-size: 13px; color: {INK2}; margin-top: 2px; line-height: 1.4; }}

/* demo-mode banner inside a platform */
[class*="st-key-omnix_demo_banner"] {{ background: #FFF6EC; border: 1px solid #F3D3B0; border-radius: 10px;
  padding: 10px 16px !important; margin-bottom: 10px; }}
.omx-banner {{ font-family: 'IBM Plex Sans', system-ui, sans-serif; font-size: 15px; color: #5A3A12; }}
.omx-banner b {{ color: #8A4B08; }}
[class*="st-key-omnix_demo_banner"] [data-testid="stButton"] button {{ background: #FFFFFF !important; border: 1px solid #E8B98A !important;
  box-shadow: none !important; padding: 5px 14px !important; min-height: 0 !important; }}
[class*="st-key-omnix_demo_banner"] [data-testid="stButton"] button p {{ color: #8A4B08 !important; font-weight: 600 !important; font-size: 14px !important; }}

/* home */
.omx-lead {{ font-size: 19px; color: #2C3836; line-height: 1.6; }}
.omx-about {{ width: 100%; margin: 0 0 22px 0; }}
.omx-about .hl {{ font-family: 'Sora', system-ui, sans-serif; font-size: 26px; font-weight: 700; color: {TEAL_DARK};
  line-height: 1.3; margin: 0 0 12px 0; }}
.omx-about p {{ font-size: 18px; color: #2C3836; line-height: 1.65; margin: 0 0 12px 0; }}
.omx-about p.close {{ font-weight: 600; color: {INK}; margin-bottom: 0; }}
.omx-feature-card {{ background: #FFFFFF; border: 1px solid {LINE}; border-radius: 14px; padding: 22px 22px 20px 22px;
  height: 100%; min-height: 300px; box-shadow: 0 2px 10px rgba(17,24,22,0.04); margin-bottom: 18px; }}
.omx-feature-icon {{ width: 48px; height: 48px; border-radius: 12px; background: #E3F2EF; display: flex;
  align-items: center; justify-content: center; margin-bottom: 12px; }}
.omx-feature-card h3 {{ margin: 0 0 6px 0 !important; padding: 0 !important; font-size: 19px !important; color: {INK}; }}
.omx-feature-card p {{ margin: 0 !important; font-size: 15.5px; color: {INK2}; line-height: 1.55; }}
.omx-value {{ background: #FFFFFF; border: 1px solid {LINE}; border-radius: 12px; padding: 18px 20px; height: 100%; }}
.omx-value h4 {{ margin: 10px 0 6px 0 !important; padding: 0 !important; font-size: 18px !important; }}
.omx-value p {{ margin: 0 !important; font-size: 15px; color: {INK2}; line-height: 1.5; }}
[class*="st-key-omnix_card_"] {{ background: #FFFFFF; border: 1px solid {LINE} !important; border-radius: 14px !important;
  padding: 26px 26px 22px 26px !important; }}
.omx-card-icon {{ width: 52px; height: 52px; border-radius: 12px; background: #E3F2EF; display: flex;
  align-items: center; justify-content: center; margin-bottom: 14px; }}
.omx-card h2 {{ margin: 0 !important; padding: 0 !important; font-size: 26px !important; font-weight: 700 !important; }}
.omx-card .sub {{ font-size: 16px; color: {INK2}; margin: 4px 0 14px 0; }}
.omx-chip {{ display: inline-block; font-size: 14px; font-weight: 500; color: {TEAL_DARK}; background: #E3F2EF;
  border-radius: 999px; padding: 4px 12px; margin: 0 6px 6px 0; }}
.omx-card ul {{ margin: 10px 0 4px 0; padding-left: 20px; }}
.omx-card li {{ font-size: 16px; color: #2C3836; margin: 6px 0; line-height: 1.4; }}
.omx-note {{ font-size: 15px; color: {INK2}; margin-top: 18px; }}

/* portal action buttons (cards, pages) */
[class*="st-key-omnix_open_"] [data-testid="stButton"] button,
[class*="st-key-omnix_action"] [data-testid="stButton"] button {{
  background: {TEAL} !important; border: none !important; border-radius: 9px !important;
  padding: 11px 22px !important; box-shadow: none !important; }}
[class*="st-key-omnix_open_"] [data-testid="stButton"] button p,
[class*="st-key-omnix_action"] [data-testid="stButton"] button p {{
  color: #FFFFFF !important; font-family: 'Sora', system-ui, sans-serif !important; font-weight: 600 !important; font-size: 16px !important; }}
[class*="st-key-omnix_open_"] [data-testid="stButton"] button:hover,
[class*="st-key-omnix_action"] [data-testid="stButton"] button:hover {{ background: {TEAL_DARK} !important; }}
[data-testid="stFormSubmitButton"] button {{ background: {TEAL} !important; border: none !important; padding: 8px 26px !important; }}
[data-testid="stFormSubmitButton"] button p {{ color: #FFFFFF !important; font-weight: 600 !important; }}
[class*="st-key-omnix_quiet"] [data-testid="stButton"] button {{ background: #FFFFFF !important; border: 1px solid #CBD5E1 !important;
  box-shadow: none !important; }}
[class*="st-key-omnix_quiet"] [data-testid="stButton"] button p {{ color: #334155 !important; font-weight: 600 !important; }}

/* help pages */
[class*="st-key-omnix_tut_"] [data-testid="stButton"] button {{ background: #F8FAFC !important; border: 1px solid #E2E8F0 !important;
  border-radius: 8px !important; box-shadow: none !important; padding: 8px 18px !important; }}
[class*="st-key-omnix_tut_"] [data-testid="stButton"] button p {{ color: #334155 !important; font-weight: 600 !important; font-size: 17px !important; }}
[class*="st-key-omnix_tut_{page[5:] if page.startswith('help-') else ''}"] [data-testid="stButton"] button {{ background: {ORANGE} !important; border-color: {ORANGE} !important; }}
[class*="st-key-omnix_tut_{page[5:] if page.startswith('help-') else ''}"] [data-testid="stButton"] button p {{ color: #FFFFFF !important; }}
.omx-step {{ display: flex; gap: 16px; align-items: flex-start; background: #FFFFFF; border: 1px solid {LINE};
  border-radius: 12px; padding: 16px 18px; margin-bottom: 10px; }}
.omx-step .n {{ flex-shrink: 0; width: 34px; height: 34px; border-radius: 50%; background: {TEAL}; color: #FFFFFF;
  font-family: 'Sora', system-ui, sans-serif; font-weight: 700; display: flex; align-items: center; justify-content: center; }}
.omx-step .h {{ font-family: 'Sora', system-ui, sans-serif; font-weight: 700; font-size: 17px; color: {INK}; }}
.omx-step .w {{ font-size: 13px; color: {INK2}; margin: 1px 0 4px 0; }}
.omx-step .b {{ font-size: 15px; color: #2C3836; line-height: 1.5; }}
.omx-badge {{ display: inline-block; font-size: 12px; font-weight: 700; letter-spacing: 0.4px; border-radius: 999px;
  padding: 2px 10px; margin-left: 8px; vertical-align: middle; font-family: 'IBM Plex Sans', system-ui, sans-serif; }}
.omx-badge.req {{ background: #E3F2EF; color: {TEAL_DARK}; }} .omx-badge.opt {{ background: #EEF2F7; color: #475569; }}
[data-testid="stTab"] {{ padding: 6px 14px 10px 0 !important; margin-right: 14px; }}
[data-testid="stTab"] p {{ font-family: 'Sora', system-ui, sans-serif !important; font-size: 18px !important; font-weight: 600 !important; }}
[data-testid="stTab"][data-selected="true"] p {{ color: {TEAL_DARK} !important; }}
table.omx-spec {{ border-collapse: collapse; width: 100%; font-size: 15px; margin: 6px 0 4px 0; }}
table.omx-spec th {{ text-align: left; font-weight: 600; color: {INK2}; border-bottom: 2px solid {LINE}; padding: 7px 10px; }}
table.omx-spec td {{ border-bottom: 1px solid {LINE}; padding: 7px 10px; vertical-align: top; color: #2C3836; }}
table.omx-spec td code {{ font-size: 13px; }}

[class*="st-key-omnix_plan_"] {{ background: #FFFFFF; border-radius: 16px !important; padding: 6px 4px 10px 4px; }}
.st-key-omnix_plan_multi {{ border: 2px solid {TEAL} !important; box-shadow: 0 8px 24px rgba(0,105,92,0.12); }}
.omx-plan h3 {{ margin: 0 !important; padding: 0 !important; font-size: 22px !important; }}
.omx-plan .aud {{ font-size: 14px; color: {INK2}; margin-top: 2px; min-height: 22px; line-height: 1.45; }}
.omx-plan .badge {{ display: inline-block; background: {ORANGE}; color: #FFF; font-size: 12px; font-weight: 700;
  letter-spacing: 0.5px; text-transform: uppercase; border-radius: 999px; padding: 3px 10px; margin-bottom: 8px; }}
.omx-plan .badge.none {{ visibility: hidden; }}
.omx-plan .price {{ font-family: 'Sora', system-ui, sans-serif; font-size: 34px; font-weight: 700; color: {TEAL_DARK};
  margin: 14px 0 0 0; line-height: 1.1; }}
.omx-plan .price span {{ font-size: 15px; font-weight: 500; color: {INK2}; }}
.omx-plan .per {{ font-size: 13.5px; color: {INK2}; min-height: 22px; margin: 4px 0 8px 0; }}
.omx-plan ul {{ margin: 8px 0 6px 0; padding-left: 20px; min-height: 290px; }} .omx-plan li {{ margin: 6px 0; font-size: 15px; }}
.omx-compare {{ width: 100%; border-collapse: collapse; font-size: 15px; background: #FFFFFF; }}
.omx-compare th, .omx-compare td {{ border-bottom: 1px solid {LINE}; padding: 10px 12px; text-align: center; }}
.omx-compare th {{ font-family: 'Sora', system-ui, sans-serif; font-size: 15px; color: {INK}; }}
.omx-compare td:first-child, .omx-compare th:first-child {{ text-align: left; }}
.omx-compare .y {{ color: {TEAL}; font-weight: 700; }} .omx-compare .n {{ color: #B0B7BC; }}
.omx-notes {{ font-size: 14px; color: {INK2}; margin-top: 10px; }} .omx-notes li {{ margin: 3px 0; }}
{account_page.css()}
.omx-tier {{ background: #FFFFFF; border: 1px solid {LINE}; border-radius: 14px; padding: 22px 24px; height: 100%; }}
.omx-tier h3 {{ margin: 0 0 4px 0 !important; padding: 0 !important; font-size: 21px !important; }}
.omx-tier .k {{ font-size: 13px; letter-spacing: 0.6px; text-transform: uppercase; font-weight: 700; color: {INK2}; }}
.omx-tier ul {{ margin: 10px 0 0 0; padding-left: 20px; }} .omx-tier li {{ margin: 5px 0; font-size: 15px; }}
[class*="st-key-omnix_footer"] {{ background: linear-gradient(135deg, #00574D 0%, #00695C 100%); border-radius: 16px; padding: 38px 40px 30px 40px !important;
  margin-top: 56px; }}
[class*="st-key-omnix_footer"] [data-testid="stButton"] button {{ background: transparent !important; border: none !important;
  box-shadow: none !important; padding: 3px 0 !important; min-height: 0 !important; justify-content: flex-start !important; }}
[class*="st-key-omnix_footer"] [data-testid="stButton"] button p {{ color: #E6EEEC !important; font-size: 16px !important;
  font-weight: 400 !important; font-family: 'IBM Plex Sans', system-ui, sans-serif !important; }}
[class*="st-key-omnix_footer"] [data-testid="stButton"] button:hover p {{ color: #7FD3C3 !important; text-decoration: underline; }}
[class*="st-key-omnix_footer"] [data-testid="stVerticalBlock"] {{ gap: 6px !important; }}
.omx-fbrand {{ font-family: 'Sora', system-ui, sans-serif; font-weight: 800; font-size: 34px; color: #FFFFFF; letter-spacing: -0.8px; }}
.omx-fbrand sup {{ font-size: 13px; color: #7FD3C3; font-weight: 600; }}
.omx-ftag {{ color: #7FD3C3; font-size: 14px; margin: 2px 0 20px 0; }}
.omx-fhead {{ font-family: 'Sora', system-ui, sans-serif; font-weight: 700; font-size: 17px; color: #FFFFFF; margin-bottom: 8px; }}
[class*="st-key-omnix_footer"] .ln {{ color: #E6EEEC; font-size: 15.5px; margin: 6px 0; }}
[class*="st-key-omnix_footer"] .ln.addr {{ color: #B9C7C3; line-height: 1.5; }}
.omx-fcopy {{ color: #8FA19C; font-size: 14px; margin-top: 22px; }}
.omx-dp-h {{ font-family: 'Sora', system-ui, sans-serif; font-weight: 700; font-size: 20px; color: {INK}; margin: 22px 0 2px 0; }}
.omx-dp-s {{ font-size: 15px; color: {INK2}; margin-bottom: 8px; }}
.omx-dp-t {{ font-family: 'Sora', system-ui, sans-serif; font-weight: 600; font-size: 16.5px; color: {TEAL_DARK}; }}
.omx-dp-a {{ font-size: 14.5px; color: #2C3836; margin-top: 2px; }}
.omx-dp-f {{ display: block; font-size: 13px; color: {INK2}; margin-top: 3px; }}
[class*="st-key-omnix_dp_"] [data-testid="stDownloadButton"] button {{ background: {TEAL} !important; border: none !important; }}
[class*="st-key-omnix_dp_"] [data-testid="stDownloadButton"] button p,
[class*="st-key-omnix_dp_"] [data-testid="stDownloadButton"] button span {{ color: #FFFFFF !important; font-weight: 600 !important; }}
[class*="st-key-omnix_action_dpzip_"] [data-testid="stDownloadButton"] button {{ background: {TEAL} !important; border: none !important; }}
[class*="st-key-omnix_action_dpzip_"] [data-testid="stDownloadButton"] button p {{ color: #FFFFFF !important; font-weight: 600 !important; }}
.omx-footer {{ margin-top: 48px; padding: 22px 4px; border-top: 1px solid {LINE}; font-size: 15px; color: {INK2};
  font-family: 'IBM Plex Sans', system-ui, sans-serif; }}
.omx-footer b {{ font-family: 'Sora', system-ui, sans-serif; color: {INK}; }}

/* The platform apps style every focused button's text white (they have teal buttons). A header, menu or
   footer item keeps focus after it is clicked, so on a platform page the clicked tab turned white on white
   and "disappeared". These rules (more specific than the platforms') keep the portal's own colors. */
.stApp [class*="st-key-omnix_brand"] [data-testid="stButton"] button:focus:not(:active) p {{ color: {INK} !important; }}
.stApp [class*="st-key-omnix_brand"] [data-testid="stButton"] button:focus:not(:active):hover p {{ color: {TEAL_DARK} !important; }}
.stApp [class*="st-key-omnix_nav"] [data-testid="stButton"] button:focus:not(:active) p {{ color: #17211F !important; }}
.stApp [class*="st-key-omnix_nav"] [data-testid="stButton"] button:focus:not(:active):hover p {{ color: {ORANGE} !important; }}
.stApp [class*="st-key-omnix_nav_{active}"] [data-testid="stButton"] button:focus:not(:active) p {{ color: {TEAL_DARK} !important; }}
.stApp [class*="st-key-omnix_nav_account"] [data-testid="stButton"] button:focus:not(:active) p {{ color: #FFFFFF !important; }}
.stApp [data-testid="stPopoverBody"] [data-testid="stButton"] button:focus:not(:active) p {{ color: {TEAL_DARK} !important; }}
.stApp [data-testid="stPopoverBody"] [data-testid="stButton"] button:focus:not(:active):hover p {{ color: {ORANGE} !important; }}
.stApp [class*="st-key-omnix_footer"] [data-testid="stButton"] button:focus:not(:active) p {{ color: #E6EEEC !important; }}
.stApp [class*="st-key-omnix_footer"] [data-testid="stButton"] button:focus:not(:active):hover p {{ color: #7FD3C3 !important; }}
.stApp [class*="st-key-omnix_demo_banner"] [data-testid="stButton"] button:focus:not(:active) p {{ color: #8A4B08 !important; }}
</style>
""",
        unsafe_allow_html=True,
    )


# ---------------------------------------------------------------------------
# header, menus, side panel
# ---------------------------------------------------------------------------
def _contact_lines():
    """Contact rows from config.CONTACT; empty fields are not shown."""
    c = config.CONTACT
    street, _, city = c.get("address", "").partition(", ")
    address = "<br>".join(html.escape(x) for x in (c.get("organization", ""), street, city) if x)
    rows = [("Email", html.escape(c.get("email", ""))), ("Phone", html.escape(c.get("phone", ""))),
            ("Address", address)]
    return "".join(f'<div class="omx-kv"><div><div class="k">{k}</div><div class="v">{v}</div></div></div>'
                   for k, v in rows if v)


def _nav():
    with st.container(key="omnix_nav", horizontal=True, gap=None, vertical_alignment="bottom"):
        for k in APPS:
            with st.container(key=f"omnix_nav_{k}", width="content"):
                st.button(config.PLATFORMS[k]["name"], key=f"omnix_navbtn_{k}", on_click=go, args=(k,))
        with st.container(key="omnix_nav_subscription", width="content"):
            st.button("Subscription", key="omnix_navbtn_subscription", on_click=go, args=("subscription",))
        with st.container(key="omnix_navpop_contact", width="content"):
            with st.popover("Contact"):
                st.markdown('<div class="omx-menu-title">Contact information</div>' + _contact_lines(),
                            unsafe_allow_html=True)
                st.button("Open the Contact page", key="omnix_menu_contact", on_click=go, args=("contact",))
        with st.container(key="omnix_navpop_help", width="content"):
            with st.popover("Help"):
                st.markdown('<div class="omx-menu-title">Tutorials</div>', unsafe_allow_html=True)
                for k in APPS:
                    st.button(f"About {config.PLATFORMS[k]['name']}", key=f"omnix_menu_help_{k}", on_click=go,
                              args=(f"help-{k}",))
                st.markdown('<div class="omx-menu-title" style="margin-top:10px">Example data</div>',
                            unsafe_allow_html=True)
                st.button("Data preparation (Excel demo files)", key="omnix_menu_dataprep", on_click=go,
                          args=("data-preparation",))
        with st.container(key="omnix_nav_account", width="content"):
            st.button("Account" if account.current_user() else "Login", key="omnix_navbtn_account", on_click=go,
                      args=("account",), icon=":material/account_circle:")


def header(compact=False):
    """Omnix wordmark (click: back to the home page), tagline and the main menu."""
    with st.container(key="omnix_header_compact" if compact else "omnix_header"):
        if compact:
            # logo, then the wordmark with the tagline directly underneath it
            with st.container(horizontal=True, vertical_alignment="center", gap="small", key="omnix_brandrow_compact"):
                if _logo("compact"):
                    st.markdown(_logo("compact"), unsafe_allow_html=True, width="content")
                with st.container(gap=None, width="content"):
                    with st.container(key="omnix_brand_compact", width="content"):
                        st.button(f"{config.BRAND}™", key="omnix_brandbtn_compact", on_click=go, args=("home",),
                                  help="Omnix home")
                    st.markdown(f'<div class="omx-tagline compact">{config.TAGLINE}</div>', unsafe_allow_html=True)
        else:
            with st.container(horizontal=True, vertical_alignment="center", gap="small", key="omnix_brandrow"):
                if _logo("full"):
                    st.markdown(_logo("full"), unsafe_allow_html=True, width="content")
                with st.container(gap=None):
                    with st.container(key="omnix_brand", width="content"):
                        st.button(f"{config.BRAND}™", key="omnix_brandbtn", on_click=go, args=("home",),
                                  help="Omnix home")
                    st.markdown(f'<div class="omx-tagline">{config.TAGLINE}</div>', unsafe_allow_html=True)
        _nav()


def access_panel():
    """Access status in the side panel (every page)."""
    s = access.status()
    if s["subscribed"]:
        detail = " · ".join(x for x in (s.get("plan"), f"until {s['expires']}" if s.get("expires") else "") if x)
        if len(s["platforms"]) < len(APPS):
            detail += f"<br>Your data: {_platform_names(s['platforms'])}"
        who = f"<div class='d'>{html.escape(s['name'])}</div>" if s.get("name") else ""
        st.sidebar.markdown(f"<div class='omx-access'><div class='t'>Your access</div><div class='s on'>{'Owner' if access.is_admin() else 'Subscriber'}</div>"
                            f"<div class='d'>{detail}</div>{who}</div>", unsafe_allow_html=True)
    else:
        st.sidebar.markdown("<div class='omx-access'><div class='t'>Your access</div><div class='s off'>Demo</div>"
                            "<div class='d'>All analyses on the built-in demo datasets. Subscribe to analyze your "
                            "own data.</div></div>", unsafe_allow_html=True)
        with st.sidebar.container(key="omnix_action_side_sub"):
            st.button("Subscribe or enter access key", key="omnix_side_subscribe", on_click=go, args=("subscription",),
                      use_container_width=True)


def portal_sidebar(page):
    """Side panel for the portal's own pages (the platforms bring their own)."""
    st.sidebar.markdown(f"<div class='omx-side-head'>{_logo('side')}<div><div class='omx-side-brand'>{config.BRAND}™</div>"
                        f"<div class='omx-side-tag'>{config.TAGLINE}</div></div></div>",
                        unsafe_allow_html=True)
    items = [("Home", "home"), *[(config.PLATFORMS[k]["name"], k) for k in APPS], ("Subscription", "subscription"),
             ("Contact", "contact"), ("Help", "help-" + (page[5:] if page.startswith("help-") else APPS[0])),
             ("Data preparation", "data-preparation"), ("Account" if account.current_user() else "Login", "account")]
    if access.is_admin():
        items.append(("Owner console", "owner"))
    with st.sidebar.container(key="omnix_side_nav"):
        for label, target in items:
            with st.container(key=f"omnix_sidebtn_{_slug(target.split('-')[0] if target.startswith('help-') else target)}"):
                st.button(label, key=f"omnix_side_{_slug(label)}", on_click=go, args=(target,), use_container_width=True)
    access_panel()


def _platform_names(keys):
    names = [config.PLATFORMS[k]["name"] for k in keys if k in config.PLATFORMS]
    return names[0] if len(names) == 1 else ", ".join(names[:-1]) + " and " + names[-1] if names else "no platform"


def demo_banner(page):
    if access.is_subscriber(page):
        return
    name = config.PLATFORMS[page]["name"]
    s = access.status()
    with st.container(key="omnix_demo_banner"):
        c1, c2 = st.columns([5, 2], vertical_alignment="center")
        if s["subscribed"]:            # subscribed, but the plan does not include this platform
            c1.markdown(f"<div class='omx-banner'><b>Demo mode for {html.escape(name)}.</b> Your "
                        f"{html.escape(s.get('plan') or 'current')} plan covers {html.escape(_platform_names(s['platforms']))}. "
                        f"Upgrade to Multi-Omics to analyze your own {html.escape(name)} data.</div>",
                        unsafe_allow_html=True)
            with c2:
                st.button("Upgrade your plan", key="omnix_banner_upgrade", on_click=go, args=("subscription",))
        else:
            c1.markdown(f"<div class='omx-banner'><b>Demo mode.</b> You can run every {html.escape(name)} analysis on "
                        "the built-in demo datasets. Uploading your own data requires an Omnix subscription.</div>",
                        unsafe_allow_html=True)
            with c2:
                st.button("Subscribe or enter access key", key="omnix_banner_sub", on_click=go, args=("subscription",))


def footer():
    """Site footer: brand and contact on the left, three columns of links (buttons, so the visitor's
    session, e.g. an activated key, is kept)."""
    c = config.CONTACT
    year = dt.date.today().year
    columns = [
        ("Platforms", [(config.PLATFORMS[k]["name"], k) for k in APPS]),
        ("Product", [("Login" if not account.current_user() else "Account", "account"),
                     ("Pricing", "subscription"), ("Get a quote", "contact"), ("FAQs", "subscription")]),
        ("Resources", [("About Omnix", "home"), ("Help center", f"help-{APPS[0]}"),
                       ("Data preparation", "data-preparation"), ("Contact", "contact")]),
    ]
    with st.container(key="omnix_footer"):
        left, *cols = st.columns([1.5, 1, 1, 1], gap="large")
        street, _, city = c.get("address", "").partition(", ")
        lines = "".join(f"<div class='ln'>{html.escape(x).replace('@', '@<wbr>')}</div>"
                        for x in (c.get("email", ""), c.get("phone", "")) if x)
        addr = "<br>".join(html.escape(x) for x in (street, city) if x)
        badge = f"<span class='omx-fbadge'>{_logo('foot')}</span>" if _logo("foot") else ""
        left.markdown(f"<div class='omx-fhead-row'>{badge}<div class='omx-fbrand'>{config.BRAND}<sup>™</sup></div></div>"
                      f"<div class='omx-ftag'>{html.escape(config.TAGLINE)}</div>{lines}"
                      + (f"<div class='ln addr'>{addr}</div>" if addr else "")
                      + f"<div class='omx-fcopy'>© {year} {config.BRAND}. All rights reserved.</div>",
                      unsafe_allow_html=True)
        for col, (title, links) in zip(cols, columns):
            with col:
                st.markdown(f"<div class='omx-fhead'>{title}</div>", unsafe_allow_html=True)
                for i, (label, target) in enumerate(links):
                    if label == "Get a quote":
                        st.button(label, key=f"omnix_foot_{_slug(title)}_{i}", on_click=contact.set_reason,
                                  args=("Request a quote",))
                    else:
                        st.button(label, key=f"omnix_foot_{_slug(title)}_{i}", on_click=go, args=(target,))


# ---------------------------------------------------------------------------
# home
# ---------------------------------------------------------------------------
def page_home():
    st.markdown("## About Omnix")
    about = [
        "Omnix unifies metabolomics, proteomics, and transcriptomics in a single, streamlined analysis platform. "
        "Transform complex datasets into meaningful biological insights through an intuitive, guided workflow—from "
        "data upload and quality control to normalization, statistical analysis, visualization, and pathway "
        "interpretation.",
        "Built for modern life-science research, Omnix makes every analysis transparent, reproducible, and "
        "publication-ready. Understand what happens at every step, explore your results with confidence, and export "
        "high-quality figures and tables directly for your research and manuscripts.",
    ]
    st.markdown("<div class='omx-about'><div class='hl'>One platform. Three omics. From data to discovery.</div>"
                + "".join(f"<p>{html.escape(t)}</p>" for t in about)
                + "<p class='close'>From raw data to biological discovery—Omnix brings the entire analysis workflow "
                  "together.</p></div>", unsafe_allow_html=True)
    features = [
        ("workflow", "One guided workflow",
         "A consistent analysis experience across metabolomics, proteomics, and transcriptomics, from data upload "
         "and quality control to statistics, visualization, and biological interpretation."),
        ("stats", "Statistics you can report",
         "QC-based filtering, established normalization methods, statistical testing, multiple-testing correction, "
         "and clearly documented analysis steps designed for research reporting."),
        ("figure", "Publication-ready output",
         "Export high-resolution figures and analysis tables in commonly used formats, making it easy to move from "
         "results to reports, presentations, and manuscripts."),
    ]
    col1, col2, col3 = st.columns(3, gap="medium")
    for col, (icon, title, text) in zip((col1, col2, col3), features):
        # one line of HTML: indented multi-line HTML would be shown as a code block by Markdown
        col.markdown(f'<div class="omx-feature-card"><div class="omx-feature-icon">{_svg(ICONS[icon], 26)}</div>'
                     f'<h3>{html.escape(title)}</h3><p>{html.escape(text)}</p></div>', unsafe_allow_html=True)

    st.markdown("## Choose an analysis platform")
    cols = st.columns(3, gap="medium")
    for col, k in zip(cols, APPS):
        p = config.PLATFORMS[k]
        with col, st.container(key=f"omnix_card_{k}", border=True, height="stretch"):
            chips = "".join(f'<span class="omx-chip">{html.escape(c)}</span>' for c in p["chips"])
            items = "".join(f"<li>{html.escape(h)}</li>" for h in p["highlights"])
            st.markdown(f'<div class="omx-card"><div class="omx-card-icon">{_svg(ICONS[k])}</div>'
                        f'<h2>{p["name"]}</h2><div class="sub">{html.escape(p["subtitle"])}</div>'
                        f'<div>{chips}</div><ul>{items}</ul></div>', unsafe_allow_html=True)
            with st.container(key=f"omnix_open_{k}"):
                st.button(f"Open {p['name']}  →", key=f"omnix_openbtn_{k}", on_click=go, args=(k,))
    st.markdown("<div class='omx-note'>Try any platform now on its built-in demo datasets. "
                "A subscription lets you upload and analyze your own data.</div>", unsafe_allow_html=True)


# ---------------------------------------------------------------------------
# subscription
# ---------------------------------------------------------------------------
def page_subscription():
    st.markdown("## Subscription")
    st.markdown(f"<p class='omx-lead'>{html.escape(config.SUBSCRIPTION_INTRO)}</p>", unsafe_allow_html=True)

    c1, c2 = st.columns(2, gap="medium")
    c1.markdown("<div class='omx-tier'><div class='k'>Without a subscription</div><h3>Demo access</h3><ul>"
                "<li>Every analysis in all three platforms</li><li>Built-in demo datasets</li>"
                "<li>All figures and tables can be downloaded</li><li>Full tutorials and example files</li></ul></div>",
                unsafe_allow_html=True)
    c2.markdown("<div class='omx-tier'><div class='k'>With a subscription</div><h3>Your own data</h3><ul>"
                "<li>Upload your own data, metadata and annotation files</li>"
                "<li>Metabolomics, Proteomics and Transcriptomics</li>"
                "<li>Everything included in demo access</li></ul></div>", unsafe_allow_html=True)

    st.markdown("### Plans")
    period = st.segmented_control("Billing", ["Annual", "Monthly"], default="Annual", key="omnix_billing",
                                  label_visibility="collapsed") or "Annual"
    st.caption(f"Annual billing: {config.ANNUAL_SAVING}. Academic and non-profit pricing."
               if period == "Annual" else "Switch to annual billing to get " + config.ANNUAL_SAVING + ".")
    plans = config.SUBSCRIPTION_PLANS
    cols = st.columns(len(plans), gap="small")
    for col, plan in zip(cols, plans):
        price, per = _plan_price(plan, period)
        badge = plan.get("badge")
        feats = "".join(f"<li>{html.escape(f)}</li>" for f in plan["features"])
        with col, st.container(key=f"omnix_plan_{plan['id']}", border=True, height="stretch"):
            hide = "" if badge else " style='visibility:hidden'"
            st.markdown(f'<div class="omx-plan"><div class="badge"{hide}>{html.escape(badge or "-")}</div>'
                        f'<h3>{html.escape(plan["name"])}</h3><div class="aud">{html.escape(plan["audience"])}</div>'
                        f'<div class="price">{price}</div><div class="per">{per}</div><ul>{feats}</ul></div>',
                        unsafe_allow_html=True)
            label = "Contact sales" if plan.get("annual") is None else "Choose plan"
            with st.container(key=f"omnix_action_plan_{plan['id']}"):
                st.button(label, key=f"omnix_choose_{plan['id']}", help=f"Choose {plan['name']}", on_click=_choose_plan, args=(plan["id"],),
                          use_container_width=True)

    with st.expander("Compare plans"):
        st.markdown(_compare_table(plans), unsafe_allow_html=True)
    notes = "".join(f"<li>{html.escape(n)}</li>" for n in config.PRICING_NOTES)
    st.markdown(f"<ul class='omx-notes'>{notes}</ul>", unsafe_allow_html=True)

    st.markdown("### Activate your subscription")
    s = access.status()
    if s["subscribed"]:
        detail = " · ".join(x for x in (s.get("plan"), f"active until {s['expires']}" if s.get("expires") else "") if x)
        where = ("every platform" if len(s["platforms"]) == len(APPS)
                 else _platform_names(s["platforms"]))
        st.success(f"Your subscription is active{': ' + detail if detail else ''}. You can upload your own data in "
                   f"{where}.")
        if s.get("via") == "key":
            with st.container(key="omnix_quiet_forget"):
                if st.button("Remove the access key from this session", key="omnix_forget_key"):
                    access.forget_key()
                    st.rerun()
    else:
        st.markdown("Subscribers receive a personal access key. Enter it here to unlock your own data for this "
                    "session.")
        with st.form("omnix_key_form", border=True):
            key = st.text_input("Access key", type="password", placeholder="OMX-XXXX-XXXX-XXXX-XXXX",
                                key="omnix_key_input")
            submitted = st.form_submit_button("Activate")
        if submitted:
            ok, msg = access.redeem_key(key)
            if ok:
                st.success(msg)
                st.rerun()
            else:
                st.error(msg)
    if billing.backend() is not None:
        st.markdown("#### Or subscribe online")
        st.markdown("Sign in to your Omnix account to subscribe, download invoices, and renew or cancel at any time.")
        with st.container(key="omnix_action_to_account"):
            st.button("Go to your account", key="omnix_sub_to_account", on_click=go, args=("account",))

    st.markdown("### Questions")
    for q, a in config.PRICING_FAQ:
        with st.expander(q):
            st.markdown(a)


def _money(v):
    return f"{config.CURRENCY}{v:,.0f}" if float(v).is_integer() else f"{config.CURRENCY}{v:,.2f}"


def _plan_price(plan, billing):
    """(price html, sub-line html) for a plan card under the chosen billing period."""
    m, a = plan.get("monthly"), plan.get("annual")
    if a is None and m is None:
        return "Custom", "Tailored to your organization"
    users = f"{plan['keys']} users" if plan.get("keys", 1) and plan.get("keys", 1) > 1 else "per user"
    if billing == "Monthly" and m is not None:
        return f"{_money(m)}<span> / month</span>", f"{users}, billed monthly"
    if a is None:
        return f"{_money(m)}<span> / month</span>", f"{users}, billed monthly"
    sub = f"{users}, billed annually" if m is None else f"{_money(a / 12)} per month, billed annually"
    return f"{_money(a)}<span> / year</span>", sub


def _compare_table(plans):
    Y, N = "<span class='y'>✓</span>", "<span class='n'>–</span>"
    def users(p):
        return "Unlimited" if p.get("keys") is None else str(p["keys"])
    def plats(p):
        return "1 of 3" if p.get("platforms") == "one" else "All 3"
    rows = [("Demo analyses in all three platforms", lambda p: Y),
            ("Platforms with your own data", plats),
            ("Users (personal access keys)", users),
            ("All analyses, figures and tables", lambda p: Y),
            ("Monthly billing", lambda p: Y if p.get("monthly") is not None else N),
            ("Onboarding session for your team", lambda p: Y if p["id"] in ("lab", "enterprise") else N),
            ("Priority support", lambda p: N if p["id"] == "single" else Y),
            ("Purchase order / invoice", lambda p: Y if p["id"] in ("lab", "enterprise") else N)]
    head = "".join(f"<th>{html.escape(p['name'])}</th>" for p in plans)
    body = "".join(f"<tr><td>{html.escape(label)}</td>" + "".join(f"<td>{fn(p)}</td>" for p in plans) + "</tr>"
                   for label, fn in rows)
    return f"<table class='omx-compare'><tr><th></th>{head}</tr>{body}</table>"


def _choose_plan(plan_id):
    """Self-service plans go to the Account page when online billing is on; otherwise (and Enterprise) to Contact."""
    if plan_id in billing.SELF_SERVICE and billing.backend() is not None:
        st.session_state.omnix_acct_plan = plan_id
        st.session_state.omnix_page = "account"
        return
    name = next(p["name"] for p in config.SUBSCRIPTION_PLANS if p["id"] == plan_id)
    st.session_state.omnix_plan_interest = name
    contact.set_reason("Request a quote", name)


# ---------------------------------------------------------------------------
# owner console (role = "admin" only)
# ---------------------------------------------------------------------------
def page_owner():
    if not access.is_admin():
        st.markdown("## Owner console")
        st.info("This page is only available with an owner access key.")
        with st.container(key="omnix_action_owner_sub"):
            st.button("Enter access key", key="omnix_owner_to_sub", on_click=go, args=("subscription",))
        return
    st.markdown("## Owner console")
    st.markdown("<p class='omx-lead'>Create access keys, see who has access, and watch for key-guessing attempts. "
                "The app cannot change its own secrets: new keys, renewals and revocations take effect when you "
                "edit <b>Settings → Secrets</b> on Streamlit Cloud.</p>", unsafe_allow_html=True)
    be = billing.backend()
    names = (["Subscribers", "Create keys", "Users & trials", "Renew or revoke", "Security"]
             + (["Online subscriptions"] if be else []))
    tabs = st.tabs(names)
    t_subs, t_new, t_users, t_manage, t_sec = tabs[:5]
    with t_users:
        _owner_users_tab()
    if be:
        with tabs[5]:
            try:
                subs = be.all_subscriptions()
            except Exception as e:
                subs = []
                st.error(f"Could not read subscriptions from Stripe: {e}")
            if subs:
                live = [x for x in subs if x["status"] in billing.LIVE]
                c1, c2, c3 = st.columns(3)
                c1.metric("Active subscriptions", len(live))
                c2.metric("Cancelling at period end", sum(x["cancel_at_period_end"] for x in live))
                yearly = sum(x["amount"] * (12 if x["interval"] == "month" else 1) for x in live)
                c3.metric("Annual recurring revenue", f"${yearly:,.0f}")
                st.dataframe(pd.DataFrame([{
                    "Customer": x["owner"], "Plan": x["plan_name"],
                    "Your data in": "All three" if len(x["platforms"]) == len(APPS) else _platform_names(x["platforms"]),
                    "Billing": f"${x['amount']:,.0f} / {x['interval']}", "Status": x["status"].replace("_", " ").title()
                    + (" (cancelling)" if x["cancel_at_period_end"] and x["status"] in billing.LIVE else ""),
                    "Renews / ends": str(x["renews_on"] or ""), "Team": ", ".join(x["members"])} for x in subs]),
                    hide_index=True, use_container_width=True)
            else:
                st.info("No online subscriptions yet.")
            if be.live:
                st.link_button("Open the Stripe dashboard", "https://dashboard.stripe.com/subscriptions")
                st.caption("Refunds, disputes, coupons and tax settings are managed in the Stripe dashboard.")

    with t_subs:
        rows = access.subscribers()
        active = sum(r["active"] for r in rows)
        soon = sum(1 for r in rows if r["active"] and r["days_left"] is not None and r["days_left"] <= 30)
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Keys", len(rows))
        c2.metric("Active", active)
        c3.metric("Expiring in 30 days", soon)
        c4.metric("Expired", len(rows) - active)
        if rows:
            def state(r):
                if not r["active"]:
                    return "Expired"
                if r["days_left"] is not None and r["days_left"] <= 30:
                    return f"Expires in {r['days_left']} days"
                return "Active"
            table = pd.DataFrame([{
                "Name": r["name"], "Plan": r["plan"] + (" (owner)" if r["role"] == "admin" and r["plan"] != "Owner" else ""),
                "Your data in": "All three" if len(r["platforms"]) == len(APPS) else _platform_names(r["platforms"]),
                "Expires": r["expires"] or "No end date", "Status": state(r), "Key ID": r["key_id"],
            } for r in rows]).sort_values(["Status", "Name"])
            st.dataframe(table, hide_index=True, use_container_width=True)
            st.download_button("Download subscriber list (CSV)", table.to_csv(index=False).encode(),
                               "omnix_subscribers.csv", "text/csv", key="omnix_owner_dl")
        else:
            st.info("No access keys yet. Create them in the Create keys tab.")

    with t_new:
        plan_names = [p["name"] for p in config.SUBSCRIPTION_PLANS]
        with st.form("omnix_owner_new", border=True):
            name = st.text_input("Subscriber or lab name", key="omnix_owner_name")
            plan = st.selectbox("Plan", plan_names + ["Owner"], index=plan_names.index("Multi-Omics")
                                if "Multi-Omics" in plan_names else 0, key="omnix_owner_plan")
            plats = st.multiselect("Platforms (Single-Omics: choose one)", APPS, default=APPS,
                                   format_func=lambda k: config.PLATFORMS[k]["name"], key="omnix_owner_plats")
            c1, c2 = st.columns(2)
            count = c1.number_input("Number of keys (one per user)", 1, 200, 1, key="omnix_owner_count")
            expires = c2.date_input("Last day of access", value=dt.date.today() + dt.timedelta(days=365),
                                    min_value=dt.date.today(), key="omnix_owner_expires")
            no_end = st.checkbox("No end date", key="omnix_owner_noend")
            made = st.form_submit_button("Create keys")
        if made:
            problems = []
            if not name.strip():
                problems.append("Enter a name.")
            if not plats:
                problems.append("Choose at least one platform.")
            if plan == "Single-Omics" and len(plats) != 1:
                problems.append("Single-Omics covers exactly one platform.")
            if problems:
                st.error(" ".join(problems))
            else:
                owner = plan == "Owner"
                issued = keys.issue(name.strip(), plan, APPS if owner else plats,
                                    "" if no_end else expires.isoformat(), 1 if owner else int(count),
                                    role="admin" if owner else "")
                st.warning("These keys are shown only now and are not stored anywhere. Copy them before leaving "
                           "this page.")
                st.markdown("**1. Send each key to its user privately**")
                st.code("\n".join(k for k, _ in issued), language=None)
                st.markdown("**2. Paste these lines into Settings → Secrets, under `[omnix.access_keys]`, and save**")
                st.code("\n".join(line for _, line in issued), language="toml", wrap_lines=True)
                if owner:
                    st.caption("Owner keys open this console. Keep them to yourself.")

    with t_manage:
        st.markdown("Each key is one line under `[omnix.access_keys]` in **Settings → Secrets**. Its line starts "
                    "with the **Key ID** shown in the Subscribers tab.")
        st.markdown("- **Renew:** change the line's `expires` date. The user keeps the same key.\n"
                    "- **Upgrade or change platforms:** edit `plan` and `platforms` on the line "
                    "(remove `platforms` for all three).\n"
                    "- **Revoke:** delete the line. The key stops working immediately.\n"
                    "- **Lost key:** revoke the old line and create a new key; keys cannot be recovered.")
        st.markdown("#### Look up a key")
        with st.form("omnix_owner_lookup", border=True):
            probe = st.text_input("Access key", type="password", key="omnix_owner_probe")
            look = st.form_submit_button("Find")
        if look:
            kid = access.find_key(probe)
            if kid:
                st.success(f"This key is the line starting with **{kid}**.")
            else:
                st.error("This key is not in the secrets.")

    with t_sec:
        st.markdown(f"After **{access.MAX_FAILS}** wrong keys within **{access.LOCK_MINUTES} minutes**, activation is "
                    "paused for that visitor. Attempts are counted since the app last restarted.")
        rows, sessions = access.failure_report()
        if rows:
            st.dataframe(pd.DataFrame(rows, columns=["IP address", f"Wrong keys (last {access.LOCK_MINUTES} min)",
                                                     "Paused"]), hide_index=True, use_container_width=True)
        else:
            st.success("No wrong-key attempts from any IP address in the last "
                       f"{access.LOCK_MINUTES} minutes.")
        if sessions:
            st.caption(f"Browser sessions with wrong keys in that window: {sessions}.")


def _toml_q(v):
    return '"' + str(v).replace("\\", "\\\\").replace('"', '\\"') + '"'


def _owner_users_tab():
    """Owner console: give signed-in people access by email, and create trial / demo accounts."""
    st.markdown("People can get access in two ways: an **access key** (Create keys tab), or by **signing in** with "
                "Google/Microsoft using an email you list here. Paying online through Stripe needs neither: "
                "access follows their subscription automatically.")
    rows = access.email_subscribers()
    st.markdown("#### People with access by email")
    if rows:
        st.dataframe(pd.DataFrame([{
            "Email": r["email"], "Plan": r["plan"],
            "Your data in": "All three" if len(r["platforms"]) == len(APPS) else _platform_names(r["platforms"]),
            "Expires": r["expires"] or "No end date",
            "Status": "Expired" if not r["active"] else (f"Expires in {r['days_left']} days"
                                                         if r["days_left"] is not None and r["days_left"] <= 30
                                                         else "Active")} for r in rows]),
            hide_index=True, use_container_width=True)
    else:
        st.caption("Nobody yet.")

    st.markdown("#### Give someone access, or create a trial / demo account")
    plan_names = [p["name"] for p in config.SUBSCRIPTION_PLANS]
    with st.form("omnix_owner_email_form", border=True):
        c1, c2 = st.columns(2)
        who = c1.text_input("Their email (the account they sign in with)", key="omnix_owner_u_email")
        plan = c2.selectbox("Plan", ["Trial"] + plan_names, key="omnix_owner_u_plan",
                            help="Trial: full access for a limited time, for demos and evaluations.")
        plats = st.multiselect("Platforms", APPS, default=APPS, format_func=lambda k: config.PLATFORMS[k]["name"],
                               key="omnix_owner_u_plats")
        c3, c4 = st.columns(2)
        days = c3.selectbox("Length of access", [7, 14, 30, 90, 365], index=1, key="omnix_owner_u_days",
                            format_func=lambda d: f"{d} days")
        also_key = c4.checkbox("Also create an access key (for people without a Google/Microsoft account)",
                               key="omnix_owner_u_key")
        made = st.form_submit_button("Create")
    if made:
        who = who.strip().lower()
        if not account_page.EMAIL_RE.match(who):
            st.error("Enter a valid email address.")
        elif not plats:
            st.error("Choose at least one platform.")
        else:
            expires = (dt.date.today() + dt.timedelta(days=int(days))).isoformat()
            fields = [f"plan = {_toml_q(plan)}"]
            if len(plats) < len(APPS):
                fields.append("platforms = [" + ", ".join(_toml_q(p_) for p_ in plats) + "]")
            fields.append(f"expires = {_toml_q(expires)}")
            st.markdown("**Paste this line into Settings → Secrets, under `[omnix.subscriber_emails]`, and save.** "
                        f"{who} then signs in with that email and has access until {expires}.")
            st.code(f"{_toml_q(who)} = {{ {', '.join(fields)} }}", language="toml", wrap_lines=True)
            if also_key:
                issued = keys.issue(who, plan, plats, expires, 1)
                st.markdown("**Access key** (send it privately; shown only now):")
                st.code(issued[0][0], language=None)
                st.markdown("**And this line under `[omnix.access_keys]`:**")
                st.code(issued[0][1], language="toml", wrap_lines=True)
    with st.expander("How to manage people, payments and plan changes"):
        st.markdown(
            "- **Online subscribers (Stripe):** everything is automatic. Payments, invoices, receipts, renewals and "
            "failed-payment retries run in Stripe; see them in the *Online subscriptions* tab and the Stripe "
            "dashboard. Refunds and plan changes for a customer are made in the Stripe dashboard.\n"
            "- **Invoiced / purchase-order customers:** create a Stripe invoice (dashboard → Invoices → Create) or "
            "your own invoice, and when it is paid give access here (by email or key) with an end date matching "
            "the paid period.\n"
            "- **Trials and demos:** create a Trial for 7–30 days here. It ends by itself.\n"
            "- **Change a plan or renew:** edit the person's line in Secrets (plan, platforms, expires).\n"
            "- **Remove access:** delete the line in Secrets.")


# ---------------------------------------------------------------------------
# contact
# ---------------------------------------------------------------------------
def page_contact():
    st.markdown("## Contact")
    plan = st.session_state.pop("omnix_plan_interest", None)
    if plan:
        st.info(f"You selected the **{plan}** plan. Tell us the number of users and, for Single-Omics, the platform "
                "you need; we will send your quote or invoice and your access key.")
    c1, c2 = st.columns([1.6, 1], gap="large")
    with c1:
        contact.form()
    with c2:
        with st.container(border=True):
            st.markdown('<div class="omx-menu-title">Contact information</div>' + _contact_lines(),
                        unsafe_allow_html=True)
        st.markdown("#### Help and tutorials")
        st.markdown("Step-by-step guides for each platform, with example data files:")
        for k in APPS:
            with st.container(key=f"omnix_action_help_{k}"):
                st.button(f"About {config.PLATFORMS[k]['name']}", key=f"omnix_contact_help_{k}", on_click=go,
                          args=(f"help-{k}",))
        with st.container(key="omnix_action_contact_dp"):
            st.button("Data preparation (Excel example files)", key="omnix_contact_dp", on_click=go,
                      args=("data-preparation",))


# ---------------------------------------------------------------------------
# data preparation: Excel example files for every platform
# ---------------------------------------------------------------------------
_DP_TABS = {"metabolomics": "Metabolomics / Lipidomics", "proteomics": "Proteomics", "transcriptomics": "Transcriptomics"}
_DP_WHAT = {"Matrix": "The measurements: one row per feature, one column per sample.",
            "Metadata": "One row per sample: its group, whether it is a QC, its batch, and any other variables.",
            "Row annotation": "Optional information per feature, used for heatmap row tracks and pathway analysis."}


def page_data_preparation():
    st.markdown("## Data preparation")
    st.markdown("<p class='omx-lead'>Download ready-made example files in Excel for every platform. Each file is "
                "the demo dataset the platform itself uses, laid out exactly as Omnix expects. Open one, compare it "
                "with your own data, or use it as a template: replace the values and upload it as it is.</p>",
                unsafe_allow_html=True)
    tabs = st.tabs(list(_DP_TABS.values()))
    for tab, key in zip(tabs, _DP_TABS):
        with tab:
            app_path = find_app(config.PLATFORMS[key]["folder"])
            if not app_path:
                st.warning("The example files for this platform are not available on this server.")
                continue
            app_dir = os.path.dirname(app_path)
            c1, c2 = st.columns([3, 1.3], vertical_alignment="center")
            c1.markdown(f"Every workbook has two sheets: **Data**, the table to upload (Omnix reads the first "
                        f"sheet), and **Read me**, which explains each column.")
            with c2, st.container(key=f"omnix_action_dpzip_{key}"):
                st.download_button("Download all (ZIP)", demo_files.zip_all(key, app_dir),
                                   file_name=f"Omnix_{config.PLATFORMS[key]['name']}_example_files.zip",
                                   mime="application/zip", key=f"omnix_dp_zip_{key}", use_container_width=True)
            for n, (section, items) in enumerate(demo_files.DEMO[key], 1):
                req = section != "Row annotation"
                badge = ("<span class='omx-badge req'>REQUIRED</span>" if req
                         else "<span class='omx-badge opt'>OPTIONAL</span>")
                st.markdown(f"<div class='omx-dp-h'>{n}. {html.escape(_DP_TABS[key])}: {section} {badge}</div>"
                            f"<div class='omx-dp-s'>{html.escape(_DP_WHAT[section])}</div>", unsafe_allow_html=True)
                for it in items:
                    data, df = demo_files.workbook(key, it["id"], app_dir)
                    fname = demo_files.filename(key, it)
                    with st.container(border=True, key=f"omnix_dp_{it['id']}"):
                        a, b = st.columns([4, 1.3], vertical_alignment="center")
                        a.markdown(f"<div class='omx-dp-t'>{html.escape(it['title'])}</div>"
                                   f"<div class='omx-dp-a'>{html.escape(it['about'])} "
                                   f"<span class='omx-dp-f'>{html.escape(fname)} · {df.shape[0]:,} rows × "
                                   f"{df.shape[1]} columns</span></div>", unsafe_allow_html=True)
                        b.download_button("Download Excel", data, file_name=fname, key=f"omnix_dpx_{it['id']}",
                                          mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                                          icon=":material/download:", use_container_width=True)
                        with st.expander("Preview"):
                            st.dataframe(df.iloc[:6, :8], hide_index=True, width="stretch")
            st.markdown("#### Before you upload: checklist")
            st.markdown("\n".join(f"- {c}" for c in help_content.CHECKLIST))
            with st.container(key=f"omnix_action_dp_help_{key}"):
                st.button(f"Read the {config.PLATFORMS[key]['name']} guide  →", key=f"omnix_dp_help_{key}",
                          on_click=go, args=(f"help-{key}",))


# ---------------------------------------------------------------------------
# help
# ---------------------------------------------------------------------------
@functools.lru_cache(maxsize=8)
def _tutorial_sections(key, mtime):
    """The tutorial split at its '## ' headings: [(title, markdown), ...] in order."""
    with open(os.path.join(TUTORIALS, f"{key}.md"), encoding="utf-8") as fh:
        text = fh.read()
    sections, title, buf = [], None, []
    for line in text.splitlines():
        if line.startswith("## "):
            if title is not None:
                sections.append((title, "\n".join(buf).strip()))
            title, buf = line[3:].strip(), []
        else:
            buf.append(line)
    if title is not None:
        sections.append((title, "\n".join(buf).strip()))
    return sections


@functools.lru_cache(maxsize=32)
def _example(path, mtime):
    raw = open(path, "rb").read()
    df = pd.read_csv(path, nrows=8)
    for col in df.columns:                  # show True/False as written in the file, not as checkboxes
        if df[col].dtype == bool:
            df[col] = df[col].map({True: "True", False: "False"})
    return raw, df


def _spec_table(rows):
    body = "".join(f"<tr><td><b>{html.escape(c)}</b></td><td>{html.escape(r)}</td><td>{html.escape(d)}</td>"
                   f"<td><code>{html.escape(e)}</code></td></tr>" for c, r, d, e in rows)
    return ("<table class='omx-spec'><thead><tr><th>Column</th><th>Required</th><th>What goes here</th>"
            f"<th>Example</th></tr></thead><tbody>{body}</tbody></table>")


def page_help(key):
    p = config.PLATFORMS[key]
    name = p["name"]
    with st.container(horizontal=True, gap="small"):
        for k in APPS:
            with st.container(key=f"omnix_tut_{k}", width="content"):
                st.button(f"About {config.PLATFORMS[k]['name']}", key=f"omnix_tutbtn_{k}", on_click=go,
                          args=(f"help-{k}",))
    st.markdown(f"# About {name}")
    sections = _tutorial_sections(key, os.path.getmtime(os.path.join(TUTORIALS, f"{key}.md")))
    by_title = dict(sections)
    overview = by_title.get("Overview", "")
    lead = overview.split("\n\n", 1)[0]
    st.markdown(lead)
    with st.container(key=f"omnix_action_open_{key}"):
        st.button(f"Open {name}  →", key=f"omnix_help_open_{key}", on_click=go, args=(key,))

    t_quick, t_data, t_guide, t_tips = st.tabs(["Quick start", "Prepare your data", "Step-by-step guide",
                                                "Tips & glossary"])
    with t_quick:
        st.markdown(f"The {name} workflow runs top to bottom through the steps in the side panel. Each step's pages "
                    "appear as buttons at the top of the page.")
        for i, (step, pages, text) in enumerate(help_content.QUICK_START[key], 1):
            st.markdown(f"<div class='omx-step'><div class='n'>{i}</div><div><div class='h'>{html.escape(step)}</div>"
                        f"<div class='w'>{html.escape(pages)}</div></div></div>", unsafe_allow_html=True)
            st.markdown(text)
        st.info("New to the platform? Open it and use the built-in demo dataset: every step works on the demo data, "
                "so you can follow this guide before preparing your own files.")

    with t_data:
        app_path = find_app(p["folder"])
        app_dir = os.path.dirname(app_path) if app_path else None
        st.markdown("Prepare these files before uploading. Each example below is the demo dataset the platform "
                    "itself uses, so you can open it, compare it with your own file, and use it as a template.")
        st.info("All example files are also available as Excel workbooks, each with a Read me sheet, on the "
                "**Data preparation** page.")
        with st.container(key=f"omnix_action_help_dp_{key}"):
            st.button("Open Data preparation (Excel files)  →", key=f"omnix_help_dp_{key}", on_click=go,
                      args=("data-preparation",))
        for spec in help_content.FILES[key]:
            badge = "<span class='omx-badge req'>REQUIRED</span>" if spec["required"] else "<span class='omx-badge opt'>OPTIONAL</span>"
            with st.container(border=True):
                st.markdown(f"#### {html.escape(spec['title'])} {badge}", unsafe_allow_html=True)
                st.markdown(spec["purpose"])
                st.markdown(_spec_table(spec["columns"]), unsafe_allow_html=True)
                for label, fname in spec["examples"]:
                    path = os.path.join(app_dir, fname) if app_dir else None
                    if not path or not os.path.isfile(path):
                        continue
                    raw, df = _example(path, os.path.getmtime(path))
                    shown = df.iloc[:, :8]
                    st.markdown(f"**{label}** · `{fname}` · first {len(shown)} rows"
                                + (f" and {shown.shape[1]} of {df.shape[1]} columns" if df.shape[1] > 8 else ""))
                    st.dataframe(shown, hide_index=True, width="stretch")
                    st.download_button(f"Download {fname}", raw, file_name=fname, mime="text/csv",
                                       key=f"omnix_ex_{key}_{_slug(fname)}")
        st.markdown("#### Before you upload: checklist")
        st.markdown("\n".join(f"- {c}" for c in help_content.CHECKLIST))

    with t_guide:
        skip = {"Before you start: input files", "Tips & troubleshooting", "Glossary"}
        steps = [(t, body) for t, body in sections if t not in skip]
        st.markdown("Open a step to see exactly what each page does, the options it offers and what you can "
                    "download.")
        for i, (title, body) in enumerate(steps):
            with st.expander(title, expanded=(i == 0)):
                st.markdown(body)

    with t_tips:
        if "Tips & troubleshooting" in by_title:
            st.markdown("### Tips & troubleshooting")
            st.markdown(by_title["Tips & troubleshooting"])
        if "Glossary" in by_title:
            st.markdown("### Glossary")
            st.markdown(by_title["Glossary"])


# ---------------------------------------------------------------------------
# the analysis platforms
# ---------------------------------------------------------------------------
@functools.lru_cache(maxsize=8)
def _compiled_app(path, mtime):
    with open(path, encoding="utf-8") as fh:
        return compile(fh.read(), path, "exec")


@functools.lru_cache(maxsize=8)
def find_app(folder):
    """Path of <folder>/app.py. Normally apps/<folder>/app.py next to omnix_app.py, but the platforms
    may sit elsewhere in the repository (e.g. uploaded to GitHub in parts, each part in its own
    top-level folder), so the repository is searched too. A match must hold the platform's own code
    package (<folder>_modules), so an unrelated app.py is never picked up. None if not found."""
    direct = os.path.join(ROOT, "apps", folder, "app.py")
    if os.path.isfile(direct):
        return direct
    marker = f"{folder}_modules"
    bases = [ROOT, os.path.dirname(ROOT), os.path.dirname(os.path.dirname(ROOT))]
    for base in dict.fromkeys(bases):
        base_depth = base.rstrip(os.sep).count(os.sep)
        for dirpath, dirnames, filenames in os.walk(base):
            dirnames[:] = [d for d in dirnames if not d.startswith((".", "__")) and d not in ("node_modules", "venv", ".venv", "site-packages")]
            if dirpath.count(os.sep) - base_depth >= 5:
                dirnames[:] = []
            if os.path.basename(dirpath) == folder and "app.py" in filenames and marker in dirnames:
                return os.path.join(dirpath, "app.py")
    return None


def run_app(key):
    """Run one platform (apps/<folder>/app.py) exactly as `streamlit run` would, as a fresh script
    namespace on every rerun, with OMNIX set so the portal owns the page title."""
    folder = config.PLATFORMS[key]["folder"]
    path = find_app(folder)
    if path is None:
        name = config.PLATFORMS[key]["name"]
        st.error(f"The {name} platform could not be found. Omnix expects the folder `apps/{folder}` "
                 f"(containing `app.py` and `{folder}_modules/`) in the same repository as `omnix_app.py`.")
        st.markdown("Expected layout:\n```\nomnix_app.py\nomnix_portal/\napps/\n    metabolomics/\n"
                    "    proteomics/\n    transcriptomics/\n```")
        return
    code = _compiled_app(path, os.path.getmtime(path))
    namespace = {"__name__": "__main__", "__file__": path, "__builtins__": __builtins__, "OMNIX": True}
    try:
        exec(code, namespace)
    finally:
        access_panel()          # below the platform's own side panel


# ---------------------------------------------------------------------------
# entry
# ---------------------------------------------------------------------------
def main():
    access.install_upload_gate()
    page = _current_page()
    in_app = page in APPS
    isolation.activate(page if in_app else None)
    access.set_demo_only(not access.is_subscriber(page if in_app else None))
    _styles(page, portal_page=not in_app)
    header(compact=in_app)
    if in_app:
        demo_banner(page)
        run_app(page)
        return
    portal_sidebar(page)
    if page == "home":
        page_home()
    elif page == "subscription":
        page_subscription()
    elif page == "contact":
        page_contact()
    elif page == "account":
        account_page.page()
    elif page == "owner":
        page_owner()
    elif page == "data-preparation":
        page_data_preparation()
    elif page.startswith("help-"):
        page_help(page[5:])
    footer()
