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


@pytest.mark.parametrize("key, fname", [(k, Path(__file__).parent / "data" / f) for k, f in test_files.items()])
@pytest.mark.parametrize("nq, nbnd", [(2, 60)])
def test_evec_load_normality(key, fname, nq, nbnd):
    evecs = evec_load(fname, nq, nbnd)
    for q_coord, modes in evecs:
        eigs = numpy.array([modes[ibnd][1] for ibnd in range(nbnd)])
        A = numpy.abs(numpy.conj(eigs) @ eigs.T)
        assert numpy.allclose(numpy.diag(A), 1)


@pytest.mark.parametrize("key, fname", [(k, Path(__file__).parent / "data" / f) for k, f in test_files.items()])
@pytest.mark.parametrize("nq, nbnd", [(2, 60)])
def test_evec_load_orthogonal(key, fname, nq, nbnd):
    evecs = evec_load(fname, nq, nbnd)
    for q_coord, modes in evecs:
        eigs = numpy.array([modes[ibnd][1] for ibnd in range(nbnd)])
        A = numpy.abs(numpy.conj(eigs) @ eigs.T)
        if key == "eig":
            assert numpy.allclose(numpy.max(A - numpy.eye(nbnd)), 0, atol=1e-4)
        elif key == "vec":
            assert not numpy.allclose(numpy.max(A - numpy.eye(nbnd)), 0, atol=1e-4)
