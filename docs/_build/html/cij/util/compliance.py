import itertools
import numpy
from typing import List, Dict
from cij.util import c_

def reverse_moduli(
    stiffness: Dict[c_, numpy.array]
) -> numpy.ndarray:
    '''
    Calculate compliance from stiffness or stiffness from compliance

    .. math::
        S_{ij}(T, V) = C{ij}^{-1}(T, V)

    :param `modulus_keys`: array of modulus
    :param `stiffness`: 

    '''

    # Flatten C[ij](T, V) to C[T, V, i, j] form

    dims = stiffness.values()[0].shape

    _stiffness = numpy.zeros((*dims, 6, 6))

    for key in stiffness.keys():
        for i, j in set(itertools.permutations(key.voigt, 2)):
            _stiffness[:, :, i-1, j-1] = stiffness[key][:, :]

    # Matrix inverse

    _compliance = numpy.linalg.inv(_stiffness)

    # Reconstruct with S[ij](T, V) from S[T, V, i, j] form

    compliance = {}

    for i, j in itertools.product(range(6), range(6)):
        if i > j: continue
        if numpy.allclose(_compliance[:, :, i, j], 0): continue
        compliance[c_(i+1, j+1)][:, :] = _compliance[:, :, i, j]
    
    return compliance
    
