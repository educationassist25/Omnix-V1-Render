"""
Keeps each app's data separate inside one Omnix session.

The three apps were written as independent programs and use the same session-state names
(log2_data, meta, stats_result, nav_group, ...). When the user moves from one app to another,
the leaving app's data is set aside (stashed) and the arriving app's own data is put back, so
Proteomics never sees Metabolomics' matrix and each app is exactly as the user left it.

Widget values are not stashed (Streamlit rebuilds them); only the app's data is. Portal state
uses keys that start with "omnix_" and is never touched.
"""

import streamlit as st

PORTAL_PREFIX = "omnix_"
_STASH = "omnix_app_stash"
_ACTIVE = "omnix_active_app"


def _widget_keys():
    """Keys that belong to widgets in this session (None if Streamlit's layout is unknown)."""
    try:
        from streamlit.runtime.state.session_state_proxy import get_session_state
        return set(get_session_state()._state._key_id_mapper._key_id_mapping)
    except Exception:
        return None


def activate(app_key):
    """Make `app_key`'s data live (None = a portal page, which has no app data)."""
    ss = st.session_state
    current = ss.get(_ACTIVE)
    if current == app_key:
        return
    stash = ss.get(_STASH)
    if stash is None:
        stash = {}
        ss[_STASH] = stash
    widgets = _widget_keys()
    saved = {}
    for key in list(ss.keys()):
        if str(key).startswith(PORTAL_PREFIX):
            continue
        if widgets is not None and key not in widgets and current is not None:
            saved[key] = ss[key]
        del ss[key]
    if current is not None and widgets is not None:
        stash[current] = saved
    for key, value in stash.pop(app_key, {}).items() if app_key is not None else ():
        ss[key] = value
    ss[_ACTIVE] = app_key
