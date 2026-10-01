"""
Access-key creation, shared by the Owner console and tools/make_access_key.py (no Streamlit needed).

A key looks like OMX-7KQ4-M2XP-9HTD-WZ3R: 16 characters from a 32-character alphabet (80 bits),
drawn with the `secrets` module (a cryptographically secure generator). Only its SHA-256 hash is
stored in the app's secrets, so a key cannot be recovered from the secrets.
"""

import datetime as _dt
import hashlib
import secrets

ALPHABET = "ABCDEFGHJKLMNPQRSTUVWXYZ23456789"      # no 0/O/1/I, easy to read and type
PLATFORMS = ("metabolomics", "proteomics", "transcriptomics")


def new_key() -> str:
    groups = ["".join(secrets.choice(ALPHABET) for _ in range(4)) for _ in range(4)]
    return "OMX-" + "-".join(groups)


def hash_key(key: str) -> str:
    return hashlib.sha256(key.strip().upper().encode("utf-8")).hexdigest()


def _q(text) -> str:
    return '"' + str(text).replace("\\", "\\\\").replace('"', '\\"') + '"'


def secrets_line(key: str, name: str, plan: str, platforms=None, expires: str = "", role: str = "") -> str:
    """The TOML line to paste under [omnix.access_keys] for this key."""
    fields = [f"name = {_q(name)}", f"plan = {_q(plan)}"]
    plats = [p for p in (platforms or []) if p in PLATFORMS]
    if plats and len(plats) < len(PLATFORMS):
        fields.append("platforms = [" + ", ".join(_q(p) for p in plats) + "]")
    if expires:
        _dt.date.fromisoformat(str(expires))           # validates the date
        fields.append(f"expires = {_q(expires)}")
    if role:
        fields.append(f"role = {_q(role)}")
    return f'"{hash_key(key)}" = {{ {", ".join(fields)} }}'


def issue(name: str, plan: str, platforms=None, expires: str = "", count: int = 1, role: str = ""):
    """Create `count` keys. Returns [(key, secrets line), ...]."""
    out = []
    for i in range(max(1, int(count))):
        key = new_key()
        who = name if count == 1 else f"{name} · user {i + 1}"
        out.append((key, secrets_line(key, who, plan, platforms, expires, role)))
    return out
