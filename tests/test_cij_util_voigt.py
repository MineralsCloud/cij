import pytest

from cij.util import c_


@pytest.mark.parametrize( "a, b, result", [
    (c_(11),        c_(11),     True ),
    (c_(11),        c_(1111),   True ),
    (c_("11"),      c_(11),     True ),
    (c_("1111") ,   c_(11),     True ),
    (c_(11),        c_(22),     False),
    (c_(34),        c_(43),     True ),
    (c_(1112),      c_(16),     True ),
    (c_(1121),      c_(16),     True ),
    (c_(1211),      c_(16),     True ),
    (c_(56),        c_(65),     True ),
])
def test_expression(a, b, result):
    assert (a == b) == result