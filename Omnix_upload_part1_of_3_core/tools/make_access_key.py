"""
Create Omnix access keys. Run this on your own computer, never on a shared machine.

    python tools/make_access_key.py --owner --name "Your Name"
    python tools/make_access_key.py --name "Jane Doe" --plan "Multi-Omics" --expires 2027-09-30
    python tools/make_access_key.py --name "Sam Lee" --plan "Single-Omics" --platforms proteomics --expires 2027-09-30
    python tools/make_access_key.py --name "Smith Lab" --plan "Lab" --count 10 --expires 2027-09-30

--owner makes YOUR key: full access to all platforms plus the Owner console (subscriber list, key
creation, lock-out monitor). Make it once with this script; after that you can create subscriber
keys from the Owner console in the app.

Plans (omnix_portal/config.py): Single-Omics = one platform (--platforms), Multi-Omics = all three,
Lab = all three with up to 10 keys (--count 10), Enterprise = all three, as many keys as agreed.

Prints the key(s) to send, and the line(s) to paste into the app's secrets under
[omnix.access_keys] (Streamlit Cloud: App settings -> Secrets; locally: .streamlit/secrets.toml).
Only the hash goes into the secrets, so the key itself is never kept on the server.
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from omnix_portal import keys  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--name", required=True, help="subscriber or lab name (shown in the app)")
    ap.add_argument("--plan", default="Multi-Omics", help="plan name (shown in the app)")
    ap.add_argument("--platforms", default="all",
                    help="comma-separated: metabolomics, proteomics, transcriptomics, or all (default)")
    ap.add_argument("--count", type=int, default=1, help="number of keys, one per user (e.g. 10 for Lab)")
    ap.add_argument("--expires", default="", help="last day of access, YYYY-MM-DD (optional)")
    ap.add_argument("--owner", action="store_true", help="owner key: all platforms + Owner console")
    args = ap.parse_args()

    plats = [p.strip().lower() for p in args.platforms.split(",") if p.strip()]
    if plats == ["all"]:
        plats = list(keys.PLATFORMS)
    bad = [p for p in plats if p not in keys.PLATFORMS]
    if bad or not plats:
        ap.error(f"unknown platform(s) {bad}; choose from {', '.join(keys.PLATFORMS)} or all")
    if args.owner:
        args.plan, plats, args.count = "Owner", list(keys.PLATFORMS), 1
    issued = keys.issue(args.name, args.plan, plats, args.expires, args.count, role="admin" if args.owner else "")

    print("Access key(s) — send each to its user privately (this is the only time they are shown):\n")
    for key, _ in issued:
        print(f"    {key}")
    print("\nLine(s) for the app's secrets, under [omnix.access_keys]:\n")
    for _, line in issued:
        print(f"    {line}")


if __name__ == "__main__":
    main()
