import numpy
import pytest


from cij.util.units import (
    _from_ang3, _to_ang3,
    _from_gpa, _to_gpa,
    _from_ev, _to_ev
)


@pytest.mark.parametrize("x", numpy.random.rand(3, 3))
def test_units(x):

    # angstrom3 <-> bohr3
    assert numpy.allclose(_from_ang3(x), x / 0.529177249 ** 3)
    assert numpy.allclose(  _to_ang3(x), x * 0.529177249 ** 3)

    # gpa <-> ryd / bohr3
    assert numpy.allclose(_from_gpa(x), x / 14710.5076)
    assert numpy.allclose(  _to_gpa(x), x * 14710.5076)

    # eV <-> ryd
    assert numpy.allclose(_from_ev(x), x / 13.6056980659)
    assert numpy.allclose(  _to_ev(x), x * 13.6056980659)

