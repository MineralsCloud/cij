import numpy
import random
import pytest

from cij.misc.evec_sort import evec_sort

@pytest.mark.parametrize("A", numpy.random.rand(4, 3, 3))
def test_evec_sort(A):
    A = A.T @ A
    evals1, evecs1 = numpy.linalg.eig(A)
    evecs1 = evecs1
    for _ in range(3):
        pairs = list(zip(evals1, evecs1)) 
        random.shuffle(pairs)
        evals2, evecs2 = zip(*pairs)
        evals2, evecs2 = zip(*evec_sort(pairs, evecs2, evecs1))
        assert numpy.allclose(evals1, evals2)
        assert numpy.allclose(evecs1, evecs2)
