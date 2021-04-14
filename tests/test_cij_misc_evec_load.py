import pytest
from pathlib import Path
import numpy

from cij.misc.evec_load import evec_load
from cij.misc.evec_sort import evec_sort


test_files = {
    "eig": "pwscf.eig",
    "vec": "pwscf.vec"
}


@pytest.mark.parametrize("fname", [Path(__file__).parent / "data" / f for f in test_files.values()])
@pytest.mark.parametrize("nq, nbnd", [(2, 60)])
def test_evec_load_basic(fname, nq, nbnd):
    evecs = evec_load(fname, nq, nbnd)
    assert isinstance(evecs, list)
    assert len(evecs) == nq
    for q_coord, modes in evecs:
        assert len(modes) == nbnd
        for (mode_id, thz, cm_1), vec in modes:
            assert len(vec) == nbnd
