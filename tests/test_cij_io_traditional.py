import pytest
from glob import glob

from cij.io.traditional import read_energy, read_elast_data

input01s = [*glob("examples/*/input01")]
input02s = [*glob("examples/*/elast.dat")]

@pytest.fixture
def input01(request):
    return read_energy(request.param)

@pytest.fixture
def input02(request):
    return read_elast_data(request.param)

@pytest.mark.parametrize("input01", input01s, indirect=True)
def test_validate_input01(input01):
    assert input01
    assert len(input01.volumes) == input01.nv
    for v in input01.volumes:
        assert len(v.q_points) == input01.nq
        for q in v.q_points:
            assert len(q.modes) == input01.np
    assert len(input01.weights) == input01.nq
    # read_energy(input01)

@pytest.mark.parametrize("input02", input02s, indirect=True)
def test_validate_input02(input02):
    assert input02
    assert len(input02.volumes) == input02.nv
    # read_elast_data(input02)