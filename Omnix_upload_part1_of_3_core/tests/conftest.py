"""
Compatibility for Streamlit's test tool (AppTest) in Streamlit 1.52, the version Omnix runs on Render.

Two limits of that AppTest version (fixed in later Streamlit releases) would make tests fail although the app
works in the browser:

* a single-choice ``st.segmented_control`` (Annual / Monthly) is read as a list, so "Annual" was taken letter
  by letter;
* a radio button that is not on the current page any more raised KeyError when the next run read its state.

These patches only change how the test tool reads those widgets; the app code is not touched.
"""
from streamlit import __version__ as _st_version
from streamlit.testing.v1 import element_tree as _et

OLD_APPTEST = tuple(int(x) for x in _st_version.split(".")[:2]) < (1, 58)

_bg_value = _et.ButtonGroup.value.fget


def _button_group_value(self):
    v = _bg_value(self)
    if v is not None and not isinstance(v, (list, tuple)):
        return [v]                       # single-choice segmented control: one selected value
    return v


if OLD_APPTEST:
    _et.ButtonGroup.value = property(_button_group_value)

_radio_index = _et.Radio.index.fget
_radio_value = _et.Radio.value.fget


def _radio_index_safe(self):
    try:
        return _radio_index(self)
    except KeyError:                     # widget no longer in the session: report its default choice
        return int(self.proto.default) if self.options else None


def _radio_value_safe(self):
    try:
        return _radio_value(self)
    except KeyError:
        return None


_et.Radio.index = property(_radio_index_safe)
_et.Radio.value = property(_radio_value_safe)

for _cls in (_et.Selectbox,):
    _idx, _val = _cls.index.fget, _cls.value.fget

    def _index_safe(self, _idx=_idx):
        try:
            return _idx(self)
        except KeyError:
            return int(self.proto.default) if self.options and self.proto.HasField("default") else None

    def _value_safe(self, _val=_val):
        try:
            return _val(self)
        except KeyError:
            return None

    _cls.index = property(_index_safe)
    _cls.value = property(_value_safe)

# st.segmented_control: newer AppTest has `at.segmented_control`; 1.52 lists it under the button groups
from streamlit.testing.v1 import app_test as _at  # noqa: E402

if not hasattr(_at.AppTest, "segmented_control"):
    _at.AppTest.segmented_control = property(lambda self: _et.WidgetList(self._tree.get("button_group")))
if not hasattr(_et.Block, "segmented_control"):
    _et.Block.segmented_control = property(lambda self: _et.WidgetList(self.get("button_group")))
