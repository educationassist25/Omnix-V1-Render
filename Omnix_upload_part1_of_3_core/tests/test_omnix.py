"""
Omnix portal tests: every page renders, deep links work, each platform runs inside the portal,
and the platforms' data stays separate when moving between them.

    python -m pytest tests/test_omnix.py -q
"""

import os
import warnings

import pytest
from streamlit.testing.v1 import AppTest

warnings.filterwarnings("ignore")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PORTAL = os.path.join(ROOT, "omnix_app.py")
APPS = {"metabolomics": "Metabolomics", "lipidomics": "Lipidomics", "proteomics": "Proteomics",
        "transcriptomics": "Transcriptomics", "fluxomics": "Fluxomics"}


def _open(at, page):
    at.session_state["omnix_page"] = page
    return at.run()


def _step(at, group, page):
    at.session_state["nav_group"] = group
    at.session_state["nav_page"] = page
    return at.run()


def _click(at, label):
    b = next(b for b in at.button if b.label.startswith(label) or label in b.label)
    return b.click().run()


def test_home_and_portal_pages_render():
    at = AppTest.from_file(PORTAL, default_timeout=120).run()
    assert not at.exception
    labels = [b.label for b in at.button]
    for name in APPS.values():
        assert name in labels and f"Open {name}  →" in labels and f"About {name}" in labels
    assert "Subscription" in labels
    assert not at.sidebar.title                      # portal pages have no app sidebar
    for page in ["subscription", "contact", "help-metabolomics", "help-proteomics", "help-transcriptomics",
                 "help-fluxomics"]:
        _open(at, page)
        assert not at.exception, page
    text = " ".join(m.value for m in at.markdown)
    assert [t.label for t in at.tabs] == ["Quick start", "Prepare your data", "Step-by-step guide", "Tips & glossary"]


def test_menu_buttons_navigate():
    at = AppTest.from_file(PORTAL, default_timeout=120).run()
    _click(at, "About Proteomics")
    assert at.session_state["omnix_page"] == "help-proteomics" and not at.exception
    _click(at, "Open Proteomics")
    assert at.session_state["omnix_page"] == "proteomics"
    assert [t.value for t in at.sidebar.title] == ["Proteomics"]


@pytest.mark.parametrize("page", ["metabolomics", "proteomics", "transcriptomics", "fluxomics", "help-transcriptomics",
                                  "help-fluxomics", "contact"])
def test_deep_link(page):
    at = AppTest.from_file(PORTAL, default_timeout=120)
    at.query_params["page"] = page
    at.run()
    assert not at.exception and at.session_state["omnix_page"] == page


@pytest.mark.parametrize("key,name", list(APPS.items()))
def test_each_platform_runs_inside_omnix_under_its_new_name(key, name):
    at = AppTest.from_file(PORTAL, default_timeout=300).run()
    _open(at, key)
    assert not at.exception
    assert [t.value for t in at.sidebar.title] == [name]
    assert "Setup" in [b.label for b in at.sidebar.button]


def test_platform_data_stays_separate_and_is_restored():
    at = AppTest.from_file(PORTAL, default_timeout=300).run()
    # Metabolomics: load its demo
    _open(at, "metabolomics")
    _click(at, "Load & Validate Data")
    assert not at.exception
    metab = at.session_state["raw_peak_df"]
    assert metab is not None
    metab_type = at.session_state["data_type"]
    _step(at, "Normalization", "Normalization")
    # Proteomics starts clean: nothing from Metabolomics leaks in
    _open(at, "proteomics")
    assert not at.exception
    assert at.session_state["raw_peak_df"] is None
    assert at.session_state["nav_group"] == "Setup"
    _click(at, "Load & Validate Data")
    prot = at.session_state["raw_peak_df"]
    assert prot is not None and prot.shape != metab.shape
    # a portal page in between, then back to Metabolomics: its own data and page are back
    _open(at, "help-metabolomics")
    _open(at, "metabolomics")
    assert not at.exception
    assert at.session_state["raw_peak_df"].equals(metab)
    assert at.session_state["data_type"] == metab_type
    assert at.session_state["nav_group"] == "Normalization"
    # and Proteomics still has its own
    _open(at, "proteomics")
    assert at.session_state["raw_peak_df"].equals(prot)


def test_proteomics_workflow_inside_omnix():
    at = AppTest.from_file(PORTAL, default_timeout=600).run()
    _open(at, "proteomics")
    _click(at, "Load & Validate Data")
    _step(at, "Preprocessing", "Data Cleaning & Imputation")
    _click(at, "Apply Missingness Filter")
    _click(at, "Apply Imputation")
    _step(at, "Preprocessing", "QC Validation")
    _click(at, "Confirm QC & Proceed")
    _step(at, "Normalization", "Normalization")
    _click(at, "Apply Log2 Transformation")
    _click(at, "Apply Median Centering Normalization")
    _step(at, "Statistical Analysis", "Statistics")
    _click(at, "Run Two-Group Test")
    assert not at.exception
    _step(at, "Exploratory Analysis", "Volcano Plot")
    assert not at.exception
    labels = [d.label for d in at.get("download_button")]
    assert any(l.endswith("_VolcanoPlot.png") for l in labels), labels


def test_transcriptomics_workflow_inside_omnix():
    at = AppTest.from_file(PORTAL, default_timeout=900).run()
    _open(at, "transcriptomics")
    next(c for c in at.checkbox if c.label.startswith("Use built-in Rich Clinical")).check().run()
    _click(at, "Load Dataset")
    _step(at, "Preprocessing", "Data Cleaning")
    _click(at, "Run Gene Filtering")
    _step(at, "Normalization", "Normalization")
    _click(at, "Run DESeq2 Normalization")
    _step(at, "Statistical Analysis", "Statistics")
    _click(at, "Run Two-Group Test")
    assert not at.exception
    _step(at, "Exploratory Analysis", "Volcano Plot")
    assert not at.exception
    labels = [d.label for d in at.get("download_button")]
    assert any(l.endswith("_VolcanoPlot.png") for l in labels), labels


def test_metabolomics_pipeline_to_network_inside_omnix():
    at = AppTest.from_file(PORTAL, default_timeout=600).run()
    _open(at, "metabolomics")
    at.radio[0].set_value("Multiple datasets (combine methods/modes)").run()
    at.radio(key="multi_data_source").set_value("Demo Data").run()
    _click(at, "Load Demo Datasets")
    _click(at, "Auto-Process All Datasets")
    _step(at, "Normalization", "Combined Normalization Data (Multi-Dataset)")
    _click(at, "Generate Combined Normalized Data")
    assert not at.exception and at.session_state["log2_data"] is not None
    _step(at, "Network Analysis", "Correlation Rewiring Map")
    at.number_input(key="rewire_main_nperm").set_value(100)
    at.number_input(key="rewire_main_nboot").set_value(30)
    _click(at, "Run Rewiring Analysis")
    assert not at.exception and at.session_state["rewire_main_result"] is not None


# ---------------------------------------------------------------------------
# subscription access
# ---------------------------------------------------------------------------
import hashlib  # noqa: E402

TEST_KEY = "OMX-TEST-KEYS-1234-ABCD"
TEST_HASH = hashlib.sha256(TEST_KEY.encode()).hexdigest()


@pytest.fixture(autouse=True)
def _reset_guess_counter():
    # AppTest gives every test the same session id, so start each test with a clean counter
    from omnix_portal import access, billing
    access._fails.clear()
    billing._store.clear()
    yield


def _portal(secrets=None):
    at = AppTest.from_file(PORTAL, default_timeout=300)
    if secrets is not None:
        at.secrets["omnix"] = secrets
    return at.run()


def _activate(at, key):
    _open(at, "subscription")
    at.text_input(key="omnix_key_input").input(key)
    next(b for b in at.button if b.label == "Activate").click().run()
    return at


@pytest.mark.parametrize("key", list(APPS))
def test_demo_mode_locks_every_upload(key):
    at = _portal({"access_keys": {TEST_HASH: {"name": "Lab", "plan": "Standard"}}})
    _open(at, key)
    ups = at.get("file_uploader")
    assert ups and all(u.disabled for u in ups)
    assert all("subscription required" in u.label for u in ups)
    assert "Subscribe or enter access key" in [b.label for b in at.button]      # the demo banner


def test_access_key_unlocks_uploads_in_every_platform():
    at = _portal({"access_keys": {TEST_HASH: {"name": "Lab", "plan": "Standard", "expires": "2099-12-31"}}})
    _activate(at, "wrong-key")
    assert any("not recognized" in e.value for e in at.error)
    _activate(at, TEST_KEY.lower())                                            # case-insensitive
    assert any("subscription is active" in s.value for s in at.success)
    for key in APPS:
        _open(at, key)
        ups = at.get("file_uploader")
        assert ups and not any(u.disabled for u in ups), key
        assert "Subscribe or enter access key" not in [b.label for b in at.button]


def test_expired_key_is_refused():
    at = _portal({"access_keys": {TEST_HASH: {"name": "Lab", "plan": "Standard", "expires": "2020-01-01"}}})
    _activate(at, TEST_KEY)
    assert any("expired" in e.value for e in at.error)
    _open(at, "proteomics")
    assert all(u.disabled for u in at.get("file_uploader"))


def test_open_access_setting_and_no_secrets():
    at = _portal({"open_access": True})
    _open(at, "transcriptomics")
    assert not any(u.disabled for u in at.get("file_uploader"))
    at = _portal()                                                              # no secrets at all: demo only
    _open(at, "transcriptomics")
    assert all(u.disabled for u in at.get("file_uploader"))


def test_logo_returns_home_and_help_tabs():
    at = _portal()
    _open(at, "proteomics")
    next(b for b in at.button if b.label == "Omnix™").click().run()
    assert at.session_state["omnix_page"] == "home" and not at.exception
    _open(at, "help-transcriptomics")
    labels = [d.label for d in at.get("download_button")]
    assert "Download sample_counts_matrix_richdemo.csv" in labels and "Download sample_metadata_richdemo.csv" in labels
    assert len(at.dataframe) >= 3 and len(at.expander) >= 6


def test_subscription_page_shows_plans_and_billing_toggle():
    at = _portal()
    _open(at, "subscription")
    assert not at.exception
    text = " ".join(m.value for m in at.markdown)
    for name in ["Single-Omics", "Multi-Omics", "Lab", "Enterprise"]:
        assert name in text
    assert "$490" in text and "$990" in text and "$3,900" in text and "Custom" in text
    at.segmented_control(key="omnix_billing").set_value("Monthly").run()
    text = " ".join(m.value for m in at.markdown)
    assert "$49<span> / month" in text and "$99<span> / month" in text
    at.button(key="omnix_choose_multi").click().run()
    assert at.session_state["omnix_page"] == "contact"
    assert any("Multi-Omics" in i.value for i in at.info)


def test_single_platform_key_unlocks_only_that_platform():
    at = _portal({"access_keys": {TEST_HASH: {"name": "Sam", "plan": "Single-Omics", "platforms": ["proteomics"]}}})
    _activate(at, TEST_KEY)
    assert any("Proteomics" in s.value for s in at.success)
    _open(at, "proteomics")
    assert not any(u.disabled for u in at.get("file_uploader"))
    assert "Upgrade your plan" not in [b.label for b in at.button]
    for other in ["metabolomics", "transcriptomics", "fluxomics"]:
        _open(at, other)
        assert all(u.disabled for u in at.get("file_uploader")), other
        assert "Upgrade your plan" in [b.label for b in at.button]


OWNER_KEY = "OMX-OWNR-KEYS-5678-WXYZ"
OWNER_HASH = hashlib.sha256(OWNER_KEY.encode()).hexdigest()


def test_repeated_wrong_keys_pause_activation():
    at = _portal({"access_keys": {TEST_HASH: {"name": "Lab", "plan": "Multi-Omics"}}})
    for i in range(5):
        _activate(at, f"OMX-WRNG-{i}")
    _activate(at, TEST_KEY)                                                     # even the right key waits
    assert any("Too many incorrect" in e.value for e in at.error)
    _open(at, "proteomics")
    assert all(u.disabled for u in at.get("file_uploader"))


def test_owner_console_only_for_owner_and_creates_valid_keys():
    secrets = {"access_keys": {TEST_HASH: {"name": "Lab", "plan": "Multi-Omics", "expires": "2099-01-01"},
                               OWNER_HASH: {"name": "Owner", "plan": "Owner", "role": "admin"}}}
    at = _portal(secrets)
    _activate(at, TEST_KEY)                                                     # subscriber: no console
    assert "Owner console" not in [b.label for b in at.sidebar.button]
    _open(at, "owner")
    assert any("only available with an owner" in i.value for i in at.info)
    at = _portal(secrets)
    _activate(at, OWNER_KEY)
    assert "Owner console" in [b.label for b in at.sidebar.button]
    _open(at, "owner")
    assert not at.exception and len(at.dataframe) >= 1
    at.text_input(key="omnix_owner_name").input("Sam Lee")
    at.selectbox(key="omnix_owner_plan").set_value("Single-Omics")
    at.multiselect(key="omnix_owner_plats").set_value(["proteomics"])
    next(b for b in at.button if b.label == "Create keys").click().run()
    assert not at.exception, at.exception
    new_key, line = at.code[0].value.strip(), at.code[1].value.strip()
    assert new_key.startswith("OMX-") and hashlib.sha256(new_key.encode()).hexdigest() in line
    assert 'platforms = ["proteomics"]' in line and 'plan = "Single-Omics"' in line
    # the new line works once pasted into the secrets
    import tomllib
    entry = tomllib.loads("[k]\n" + line)["k"]
    secrets["access_keys"].update(entry)
    at2 = _portal(secrets)
    _activate(at2, new_key)
    _open(at2, "proteomics")
    assert not any(u.disabled for u in at2.get("file_uploader"))
    _open(at2, "metabolomics")
    assert all(u.disabled for u in at2.get("file_uploader"))


def test_open_access_is_not_owner():
    at = _portal({"open_access": True})
    _open(at, "owner")
    assert any("only available with an owner" in i.value for i in at.info)


# ---------------------------------------------------------------------------
# accounts and online billing (preview backend; Stripe backend with a stand-in client)
# ---------------------------------------------------------------------------
PREVIEW = {"billing_preview": True}
DEMO_EMAIL = "demo.researcher@example.edu"


def _sign_in_preview(at):
    _open(at, "account")
    next(b for b in at.button if b.label == "Preview with a demo account").click().run()
    assert not at.exception, at.exception
    return at


def _subscribe(at, plan_id, period="Annual", platform=None):
    at.radio(key="omnix_acct_pick_new").set_value(plan_id).run()
    at.segmented_control(key="omnix_acct_per_new").set_value(period).run()
    if platform:
        at.selectbox(key="omnix_acct_plat_new").set_value(platform).run()
    at.button(key="omnix_acct_checkout_new").click().run()
    assert not at.exception, at.exception


def test_login_tab_and_signed_out_account_page():
    at = _portal()
    labels = [b.label for b in at.button]
    assert labels.index("Login") > labels.index("Subscription")           # last item of the main menu
    _click(at, "Login")
    assert at.session_state["omnix_page"] == "account" and not at.exception
    text = " ".join(m.value for m in at.markdown)
    assert "Omnix Account" in text and "Sign in or register" in text
    assert any("being set up" in i.value for i in at.info)                  # no providers, no preview
    at = _portal(PREVIEW)
    _open(at, "account")
    assert "Preview with a demo account" in [b.label for b in at.button]


def test_preview_account_subscribe_invoice_cancel_resume():
    at = _sign_in_preview(_portal(PREVIEW))
    assert "Account" in [b.label for b in at.button]
    _open(at, "proteomics")
    assert all(u.disabled for u in at.get("file_uploader"))                # signed in, not subscribed
    _open(at, "account")
    _subscribe(at, "multi", "Monthly")
    assert any("Subscribed to Multi-Omics" in x.value for x in at.success)
    for key in APPS:
        _open(at, key)
        assert not any(u.disabled for u in at.get("file_uploader")), key
    _open(at, "account")
    labels = [d.label for d in at.get("download_button")]
    assert any(l.startswith("Download invoice OMX-") for l in labels)
    at.button(key="omnix_cancel_sub_preview_1").click().run()
    at.button(key="omnix_cancel_ok_sub_preview_1").click().run()
    assert any("Cancellation confirmed" in x.value for x in at.success)
    _open(at, "metabolomics")                                                # still has access until period end
    assert not any(u.disabled for u in at.get("file_uploader"))
    _open(at, "account")
    at.button(key="omnix_resume_sub_preview_1").click().run()
    assert any("renew automatically" in x.value for x in at.success)
    next(b for b in at.button if b.label == "Sign out").click().run()
    _open(at, "proteomics")
    assert all(u.disabled for u in at.get("file_uploader"))


def test_preview_single_omics_and_lab_team():
    at = _sign_in_preview(_portal(PREVIEW))
    _subscribe(at, "single", platform="transcriptomics")
    _open(at, "transcriptomics")
    assert not any(u.disabled for u in at.get("file_uploader"))
    _open(at, "proteomics")
    assert all(u.disabled for u in at.get("file_uploader"))
    assert "Upgrade your plan" in [b.label for b in at.button]
    # someone else's Lab plan lists the demo user as a member -> access through the team
    from omnix_portal import billing
    billing.PreviewBilling().checkout("pi@lab.edu", "PI", "lab", "year", "", "")
    sid = billing._store["pi@lab.edu"]["subs"][0]["id"]
    billing.PreviewBilling().set_members("pi@lab.edu", sid, [DEMO_EMAIL])
    billing.forget_cache  # noqa: B018
    at.session_state["omnix_page"] = "account"
    for k in [k for k in at.session_state if str(k).startswith("omnix_bill_")]:
        del at.session_state[k]
    at.run()
    _open(at, "proteomics")
    assert not any(u.disabled for u in at.get("file_uploader"))


def test_lab_owner_adds_members_and_plan_buttons_route_to_account():
    at = _portal(PREVIEW)
    _open(at, "subscription")
    at.button(key="omnix_choose_lab").click().run()
    assert at.session_state["omnix_page"] == "account"
    _sign_in_preview(at)
    _subscribe(at, "lab")
    at.text_input(key="omnix_team_email").input("postdoc@university.edu")
    next(b for b in at.button if b.label == "Add member").click().run()
    assert any("postdoc@university.edu now has access" in x.value for x in at.success)
    from omnix_portal import billing
    assert billing.PreviewBilling().memberships("postdoc@university.edu")


def test_owner_email_gets_owner_console_with_online_subscriptions():
    at = _sign_in_preview(_portal({**PREVIEW, "owner_emails": [DEMO_EMAIL]}))
    assert "Owner console" in [b.label for b in at.sidebar.button]
    _open(at, "owner")
    assert "Online subscriptions" in [t.label for t in at.tabs]


class _Obj(dict):
    def __getattr__(self, k):
        return self[k]


class _FakeStripe:
    """Stand-in for the stripe module: records calls, returns Stripe-shaped data."""
    def __init__(self):
        self.calls = []
        sub = _Obj(id="sub_1", customer="cus_A", status="active", cancel_at_period_end=False, canceled_at=None,
                   metadata={"omnix_plan": "single", "omnix_platform": "proteomics"},
                   items={"data": [{"current_period_end": 1893456000,
                                    "price": {"id": "price_single_y", "unit_amount": 49000, "currency": "usd",
                                              "recurring": {"interval": "year"}}}]})
        fake = self

        class Customer:
            @staticmethod
            def list(email, limit):
                return _Obj(data=[_Obj(id="cus_A")] if email == "a@x.edu" else [])

            @staticmethod
            def create(**kw):
                fake.calls.append(("customer.create", kw))
                return _Obj(id="cus_NEW")

        class Subscription:
            @staticmethod
            def list(**kw):
                return _Obj(data=[sub])

            @staticmethod
            def retrieve(sid):
                return sub

            @staticmethod
            def modify(sid, **kw):
                fake.calls.append(("subscription.modify", sid, kw))

            @staticmethod
            def search(**kw):
                return _Obj(data=[])

        class _Session:
            @staticmethod
            def create(**kw):
                fake.calls.append(("checkout", kw))
                return _Obj(url="https://checkout.stripe.test/s")

        class checkout:
            Session = _Session

        self.Customer, self.Subscription, self.checkout = Customer, Subscription, checkout


def test_stripe_backend_reads_plans_and_protects_other_accounts():
    from omnix_portal import billing
    be = billing.StripeBilling.__new__(billing.StripeBilling)
    be.stripe = _FakeStripe()
    be.prices = {"single_annual": "price_single_y", "multi_annual": "price_multi_y"}
    be.price_to_plan = {"price_single_y": ("single", "year"), "price_multi_y": ("multi", "year")}
    be.automatic_tax = False
    subs = be.subscriptions("a@x.edu")
    assert subs[0]["plan"] == "single" and subs[0]["platforms"] == ["proteomics"] and subs[0]["amount"] == 490
    assert str(subs[0]["renews_on"]) == "2030-01-01"
    url = be.checkout("new@x.edu", "New", "multi", "year", "", "https://omnix.app/?page=account")
    call = next(c for c in be.stripe.calls if c[0] == "checkout")[1]
    assert url.startswith("https://checkout.stripe") and call["line_items"][0]["price"] == "price_multi_y"
    assert call["customer"] == "cus_NEW" and call["success_url"].endswith("checkout=success")
    be.set_cancel("a@x.edu", "sub_1", True)
    assert ("subscription.modify", "sub_1", {"cancel_at_period_end": True}) in be.stripe.calls
    with pytest.raises(PermissionError):
        be.set_cancel("someone.else@x.edu", "sub_1", True)                   # not their subscription
    be.set_members("a@x.edu", "sub_1", ["m@x.edu"])
    md = be.stripe.calls[-1][2]["metadata"]
    assert md[billing.member_tag("m@x.edu")] == "1" and md["omnix_members"] == "m@x.edu"
    with pytest.raises(ValueError):
        be.checkout("a@x.edu", "A", "lab", "month", "", "u")                 # no such price


def test_signin_misconfiguration_is_explained_not_crashing():
    at = AppTest.from_file(PORTAL, default_timeout=120)
    at.secrets["auth"] = {"redirect_uri": "https://omnix.streamlit.app/oauth2callback",
                          "google": {"client_id": "123-abc.apps.googleusercontent.com",
                                     "server_metadata_url": "https://accounts.google.com/.well-known/openid-configuration"}}
    at.run()
    _open(at, "account")
    next(b for b in at.button if b.label == "Continue with Google").click().run()
    assert not at.exception
    msg = " ".join(e.value for e in at.error)
    assert "missing cookie_secret" in msg and "missing client_secret" in msg


def test_data_preparation_page_and_teal_footer():
    at = _portal()
    labels = [b.label for b in at.button]
    for link in ["Pricing", "Get a quote", "FAQs", "About Omnix", "Help center", "Data preparation"]:
        assert link in labels, link
    next(b for b in at.button if b.label == "Data preparation" and "foot" in str(b.key)).click().run()
    assert at.session_state["omnix_page"] == "data-preparation" and not at.exception
    assert [t.label for t in at.tabs][:5] == ["Metabolomics", "Lipidomics", "Proteomics", "Transcriptomics",
                                              "Fluxomics (isotope tracing)"]
    names = [d.label for d in at.get("download_button")]
    assert names.count("Download Excel") == 5 + 5 + 5 + 3 + 7
    assert "Download all (ZIP)" in names


def test_excel_examples_load_like_the_csv_demo_files():
    import importlib.util
    import io
    import sys
    sys.path.insert(0, ROOT)
    from omnix_portal import config, demo_files
    for plat in demo_files.DEMO:
        if plat == "fluxomics":                  # own reader: see test_fluxomics_excel_examples_read_like_csv
            continue
        folder = config.PLATFORMS[plat]["folder"]
        app_dir = os.path.join(ROOT, "apps", folder)
        spec = importlib.util.spec_from_file_location(f"u_{plat}", os.path.join(app_dir, f"{folder}_modules", "utils.py"))
        utils = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(utils)
        for _, it in demo_files.items(plat):
            data, _ = demo_files.workbook(plat, it["id"], app_dir)
            back = utils.load_table(io.BytesIO(data), demo_files.filename(plat, it))
            ref = demo_files._table(app_dir, it)
            assert back.shape == ref.shape and list(back.columns) == list(ref.columns), it["id"]


def test_contact_form_prefills_and_builds_email():
    at = _portal()
    next(b for b in at.button if b.label == "Get a quote").click().run()
    assert at.session_state["omnix_page"] == "contact" and not at.exception
    assert at.selectbox(key="omnix_contact_reason").value == "Request a quote"
    assert at.text_input(key="omnix_contact_subject").value.startswith("Quote request")
    at.selectbox(key="omnix_contact_reason").set_value("Technical support").run()
    assert at.text_input(key="omnix_contact_subject").value == "Technical support request"
    at.text_input(key="omnix_contact_name").input("Ada Lovelace")
    at.text_input(key="omnix_contact_email").input("ada@uni.edu")
    at.text_area(key="omnix_contact_message").input("The volcano plot does not load.")
    next(b for b in at.button if b.label == "Send message").click().run()
    assert not at.exception
    links = at.get("link_button")
    assert links and links[0].proto.url.startswith("mailto:education.assist25@gmail.com?subject=")


def test_contact_form_validates():
    at = _portal()
    _open(at, "contact")
    next(b for b in at.button if b.label == "Send message").click().run()
    assert any("Please add" in e.value for e in at.error)


def test_owner_creates_trial_by_email_and_it_unlocks_uploads():
    owner = {**PREVIEW, "owner_emails": [DEMO_EMAIL]}
    at = _sign_in_preview(_portal(owner))
    _open(at, "owner")
    assert "Users & trials" in [t.label for t in at.tabs]
    at.text_input(key="omnix_owner_u_email").input("trial.user@lab.org")
    next(b for b in at.button if b.label == "Create").click().run()
    assert not at.exception
    line = next(c.value for c in at.code if "trial.user@lab.org" in c.value)
    assert 'plan = "Trial"' in line and "expires" in line
    import tomllib
    entry = tomllib.loads("[e]\n" + line)["e"]
    at2 = _sign_in_preview(_portal({**PREVIEW, "subscriber_emails": {DEMO_EMAIL: entry["trial.user@lab.org"]}}))
    _open(at2, "metabolomics")
    assert not any(u.disabled for u in at2.get("file_uploader"))


# ---------------------------------------------------------------------------
# Fluxomics
# ---------------------------------------------------------------------------
def _flux_step(at, step):
    at.session_state["step"] = step
    return at.run()


def test_fluxomics_layout_matches_the_other_platforms():
    at = _portal()
    _open(at, "fluxomics")
    assert not at.exception
    assert [t.value for t in at.sidebar.title] == ["Fluxomics"]
    side = [b.label for b in at.sidebar.button]
    assert side[:7] == ["Setup", "Normalization", "Isotope Labeling", "Statistical Analysis", "Exploratory Analysis",
                        "Pathway Analysis", "Reset all steps"]                 # no Export & Report group
    main = [b.label for b in at.main.button]
    assert "Upload & inspect" in main and "Experimental design" in main     # step names, no numbers
    # a sidebar group opens its first step; the step row follows
    at.button(key="navgrp__exploratory-analysis").click().run()
    assert at.session_state["step"] == "boxplots" and not at.exception
    assert "Heatmaps" in [b.label for b in at.main.button]


def test_fluxomics_workflow_inside_omnix():
    at = AppTest.from_file(PORTAL, default_timeout=600).run()
    _open(at, "fluxomics")
    at.radio(key="in_route").set_value("Generic isotopologue table + metadata").run()
    _click(at, "Load synthetic example data")
    assert not at.exception and at.session_state["ds"] is not None
    _flux_step(at, "design")
    _click(at, "Confirm design")
    assert at.session_state["design_confirmed"]
    _flux_step(at, "normalize")
    assert not at.exception and at.session_state["norm"] is not None
    for step in ["incorporation", "boxplots", "kinetics", "pathways", "heatmap"]:
        _flux_step(at, step)
        assert not at.exception, step
    _flux_step(at, "stats")
    _click(at, "Run statistics")
    assert not at.exception and at.session_state["stats"]
    _flux_step(at, "upload")
    main = [b.label for b in at.main.button]
    assert "Upload & inspect" in main and not any(l.startswith("✓") for l in main)       # no check marks
    assert not any("Export" in l for l in main + [b.label for b in at.sidebar.button])
    _flux_step(at, "export")                                  # an old link to the removed step opens Upload
    assert not at.exception and at.session_state["step"] == "upload"


def test_fluxomics_compound_discoverer_demo_loads():
    at = AppTest.from_file(PORTAL, default_timeout=600).run()
    _open(at, "fluxomics")
    _click(at, "Load Compound Discoverer example")
    assert not at.exception and at.session_state["ds"] is not None and at.session_state["ds"].source == "cd"


def test_fluxomics_data_stays_separate_and_reset_keeps_access():
    at = _portal({"access_keys": {TEST_HASH: {"name": "Lab", "plan": "Multi-Omics", "expires": "2099-12-31"}}})
    _activate(at, TEST_KEY)
    _open(at, "fluxomics")
    assert not any(u.disabled for u in at.get("file_uploader"))           # Multi-Omics covers Fluxomics
    at.radio(key="in_route").set_value("Generic isotopologue table + metadata").run()
    _click(at, "Load synthetic example data")
    ds = at.session_state["ds"]
    _open(at, "metabolomics")
    assert "ds" not in at.session_state                                      # nothing leaks into Metabolomics
    _open(at, "fluxomics")
    assert at.session_state["ds"] is ds
    at.button(key="reset").click().run()
    assert not at.exception and at.session_state["ds"] is None
    assert at.session_state["omnix_page"] == "fluxomics"                    # still in Fluxomics ...
    assert not any(u.disabled for u in at.get("file_uploader"))           # ... with the subscription kept


def test_single_omics_fluxomics_key():
    at = _portal({"access_keys": {TEST_HASH: {"name": "Sam", "plan": "Single-Omics", "platforms": ["fluxomics"]}}})
    _activate(at, TEST_KEY)
    _open(at, "fluxomics")
    assert not any(u.disabled for u in at.get("file_uploader"))
    _open(at, "proteomics")
    assert all(u.disabled for u in at.get("file_uploader"))
    assert "Upgrade your plan" in [b.label for b in at.button]


def test_fluxomics_help_page_and_key_tool():
    import subprocess
    import sys
    at = _portal()
    _open(at, "help-fluxomics")
    assert not at.exception
    labels = [d.label for d in at.get("download_button")]
    assert "Download synthetic_isotopologue_table.csv" in labels and "Download synthetic_sample_metadata.csv" in labels
    assert len(at.expander) >= 7
    out = subprocess.run([sys.executable, os.path.join(ROOT, "tools", "make_access_key.py"), "--name", "X",
                          "--plan", "Single-Omics", "--platforms", "fluxomics"], capture_output=True, text=True)
    assert out.returncode == 0 and 'platforms = ["fluxomics"]' in out.stdout


def test_fluxomics_excel_examples_read_like_csv():
    import io
    import sys
    sys.path.insert(0, ROOT)
    from omnix_portal import demo_files, portal
    app_dir = os.path.dirname(portal.find_app("fluxomics"))
    sys.path.insert(0, app_dir)
    from core import dataset as dsm
    tables = [it for s, its in demo_files.DEMO["fluxomics"] if s == "Isotopologue table" for it in its]
    meta = next(it for s, its in demo_files.DEMO["fluxomics"] if s == "Sample metadata" for it in its)
    meta_x, _ = demo_files.workbook("fluxomics", meta["id"], app_dir)
    meta_c = open(os.path.join(app_dir, meta["csv"]), "rb").read()
    for it in tables:
        x, _ = demo_files.workbook("fluxomics", it["id"], app_dir)
        c = open(os.path.join(app_dir, it["csv"]), "rb").read()
        a = dsm.read_generic(io.BytesIO(x), "t.xlsx", io.BytesIO(meta_x), "m.xlsx", "13C")
        b = dsm.read_generic(io.BytesIO(c), "t.csv", io.BytesIO(meta_c), "m.csv", "13C")
        assert a.CT.shape == b.CT.shape and len(a.S) == len(b.S), it["id"]
    for _, it in demo_files.items("fluxomics"):
        if it.get("file"):
            data, df = demo_files.workbook("fluxomics", it["id"], app_dir)
            assert data[:2] == b"PK" and "Relative Area 0" in df.columns


# ---------------------------------------------------------------------------
# Metabolomics and Lipidomics: one app, two tabs
# ---------------------------------------------------------------------------
def test_metabolomics_and_lipidomics_tabs_offer_only_their_own_types():
    at = _portal()
    _open(at, "metabolomics")
    assert not at.exception and [t.value for t in at.sidebar.title] == ["Metabolomics"]
    types = next(s for s in at.selectbox if s.label == "Analysis type").options
    assert types == ["Untargeted Metabolomics", "Targeted Metabolomics"]
    _open(at, "lipidomics")
    assert not at.exception and [t.value for t in at.sidebar.title] == ["Lipidomics"]
    assert next(s for s in at.selectbox if s.label == "Analysis type").options == ["Untargeted Lipidomics", "Targeted Lipidomics"]


def test_lipidomics_data_stays_separate_from_metabolomics():
    at = _portal()
    _open(at, "metabolomics")
    _click(at, "Load & Validate Data")
    metab = at.session_state["raw_peak_df"]
    _open(at, "lipidomics")
    assert at.session_state["raw_peak_df"] is None and at.session_state["data_type"] == "Untargeted Lipidomics"
    _click(at, "Load & Validate Data")
    lipid = at.session_state["raw_peak_df"]
    assert lipid is not None and lipid.shape != metab.shape
    _open(at, "metabolomics")
    assert at.session_state["raw_peak_df"].equals(metab)


def test_single_omics_lipidomics_key():
    at = _portal({"access_keys": {TEST_HASH: {"name": "Sam", "plan": "Single-Omics", "platforms": ["lipidomics"]}}})
    _activate(at, TEST_KEY)
    _open(at, "lipidomics")
    assert not any(u.disabled for u in at.get("file_uploader"))
    _open(at, "metabolomics")
    assert all(u.disabled for u in at.get("file_uploader"))


@pytest.mark.parametrize("page,dtype", [("lipidomics", "Targeted Lipidomics"), ("metabolomics", "Targeted Metabolomics")])
def test_targeted_type_survives_revisiting_data_upload(page, dtype):
    """Regression: revisiting Data Upload used to reset the Analysis type to Untargeted, so loaded
    targeted data were then normalized with the untargeted Log2 -> Median-IQR pipeline."""
    at = _portal()
    _open(at, page)
    next(s for s in at.selectbox if s.label == "Analysis type").set_value(dtype).run()
    _click(at, "Load & Validate Data")
    _step(at, "Normalization", "Normalization")
    _click(at, "Apply ISTD Normalization")
    assert at.session_state["istd_normalized"] is not None
    _step(at, "Setup", "Data Upload")                       # look at the upload page again
    assert at.session_state["data_type"] == dtype
    assert next(s for s in at.selectbox if s.label == "Analysis type").value == dtype
    _step(at, "Normalization", "Normalization")
    assert any("ISTD Normalization → Log2" in h.value for h in at.subheader)
    _click(at, "Apply Strict Log2 Transformation")
    L, R = at.session_state["log2_data"], at.session_state["istd_normalized"]
    import numpy as np
    assert L is not None
    expect = np.log2(R.where(R > 0).astype(float))            # strict log2: values <= 0 become missing
    assert np.allclose(L.values, expect.values, equal_nan=True)
