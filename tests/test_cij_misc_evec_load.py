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

@pytest.fixture
def mass():
    mmap = {1: 24.305, 2: 40.078, 3: 28.0855, 4: 15.999}
    atypes = [2, 2, 1, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3]
    mass = [mmap[i] for i in atypes]
    return mass

@pytest.mark.parametrize("nq, nbnd", [(2, 60)])
def test_evec_dispvec_disp2eig_3Nx3N(nq, nbnd, mass):

    vecs = evec_load(Path(__file__).parent / "data" / test_files["vec"], nq, nbnd)
    eigs = evec_load(Path(__file__).parent / "data" / test_files["eig"], nq, nbnd)

    import numpy
    from cij.misc.evec_disp2eig import evec_disp2eig

    for iq in range(nq):

        a1 = numpy.array([vecs[iq][1][ibnd][1] for ibnd in range(nbnd)])
        a2 = numpy.array([eigs[iq][1][ibnd][1] for ibnd in range(nbnd)])
        a1 = evec_disp2eig(a1, mass)

        assert numpy.allclose(a1, a2, atol=1e-4)


@pytest.mark.parametrize("nq, nbnd", [(2, 60)])
def test_evec_dispvec_disp2eig_1x3N(nq, nbnd, mass):

    vecs = evec_load(Path(__file__).parent / "data" / test_files["vec"], nq, nbnd)
    eigs = evec_load(Path(__file__).parent / "data" / test_files["eig"], nq, nbnd)

    import numpy
    from cij.misc.evec_disp2eig import evec_disp2eig

    for iq in range(nq):

        a1 = numpy.array([vecs[iq][1][ibnd][1] for ibnd in range(nbnd)])
        a2 = numpy.array([eigs[iq][1][ibnd][1] for ibnd in range(nbnd)])

        for ibnd in range(nbnd):

            a1 = numpy.array([vecs[iq][1][ibnd][1]])
            a2 = numpy.array([eigs[iq][1][ibnd][1]])
            a1 = evec_disp2eig(a1, mass)

            assert numpy.allclose(a1[0], a2, atol=1e-4)
