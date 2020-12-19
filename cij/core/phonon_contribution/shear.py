import numpy
import itertools
from typing import Union, Tuple, List
from logging import getLogger

from lazy_property import LazyProperty

from cij.util import C_, c_


logger = getLogger(__name__)


def calculate_fictitious_strain_energy(
    fictitious_strain: numpy.ndarray,
    resolve_elastic_modulus: callable,
    target: Union[C_, None] = None
) -> numpy.ndarray:
    '''Calculate the strain energy for under given coordinate system, ignore
    the unknown.

    :param fictitious_strain: the fictitious strain where strain energy
        is calculated
    :param resolve_elastic_modulus: the function used to get elastic modulus
    :param target: the key for the elastic modulus
    '''

    _energy = 0

    nz = numpy.argwhere(
        numpy.logical_not(numpy.isclose(fictitious_strain, 0))
    )

    for (i, j), (k, l) in itertools.product(nz, nz):

        key = c_(i+1, j+1, k+1, l+1)

        if target and key == target: continue

        _moduli = resolve_elastic_modulus(key)
        _energy += _moduli * fictitious_strain[i,j] * fictitious_strain[k,l] / 2

    return _energy

def get_fictitious_strain_energy_keys(
    fictitious_strain: numpy.ndarray,
    target: Union[C_, None] = None
) -> numpy.ndarray:

    nz = numpy.argwhere(
        numpy.logical_not(numpy.isclose(fictitious_strain, 0))
    )

    _keys = []

    for (i, j), (k, l) in itertools.product(nz, nz):

        key = c_(i+1, j+1, k+1, l+1)

        if target and key == target: continue

        _keys.append(key)

    return _keys

class ShearElasticModulusPhononContribution:
    '''
    '''

    def __init__(self, strain: Tuple[float, float, float], key: C_, calculator=None):
        self.key = key
        self.strain = strain

        self.modulus_isothermal = dict()
        self.modulus_isothermal_rotated = dict()
        self.calculator = calculator

    @LazyProperty
    def fictitious_strain(self):
        '''The fictitious strain
        ''' 
        e = numpy.zeros((3, 3))

        e[self.key.i[0] - 1, self.key.i[1] - 1] = 1
        e[self.key.i[1] - 1, self.key.i[0] - 1] = 1

        e[self.key.j[0] - 1, self.key.j[1] - 1] = 1
        e[self.key.j[1] - 1, self.key.j[0] - 1] = 1

        return e
    
    @LazyProperty
    def fictitious_strain_rotated(self):
        '''The fictitious strain in the rotated coordinate system
        '''

        logger.debug(
            "Fictitious rotated -> " + \
            str(self.key) + \
            str(numpy.diag(numpy.linalg.eig(self.fictitious_strain)[0])))

        return numpy.diag(numpy.linalg.eig(self.fictitious_strain)[0])

    @LazyProperty
    def transformation_matrix(self):
        '''The transformation matrix
        '''

        logger.debug(
            "Transformation matrix -> " + \
            str(self.key) + \
            str(numpy.linalg.eig(self.fictitious_strain)[1]))

        return numpy.linalg.eig(self.fictitious_strain)[1]

    @property
    def fictitious_strain_energy(self) -> numpy.ndarray:
        '''The strain energy for the fictitious strain under the original
        coordinate system except for the unknown
        '''
        return calculate_fictitious_strain_energy(
            self.fictitious_strain,
            lambda _key: self.get_elastic_modulus(_key),
            self.key
        )

    @property
    def fictitious_strain_energy_rotated(self) -> numpy.ndarray:
        '''The strain energy for the fictitious strain under the rotated
        coordinate system
        '''
        return calculate_fictitious_strain_energy(
            self.fictitious_strain_rotated,
            lambda _key: self.get_elastic_modulus_rotated(_key)
        )
    
    @LazyProperty
    def strain_rotated(self):

        # Dimension of self.strain is (ntv, 3), now we build a (ntv, 3, 3)
        # matrix for the calculation purpose
        strain = numpy.zeros((*self.strain.shape, 3))

        # Blow up the matrix to form a diagonal (in last two indices) matrix
        # * left hand side: get a reference of the diagonal (1+1D) of the above
        #   matrix (1+2D),
        # * right hand side: strain (1+1D)
        # the diagonal is writable with [...] here
        numpy.einsum('...ii -> ...i', strain)[...] = self.strain

        # Calculate the rotated strain
        strain_rotated = self.transformation_matrix.T @ strain @ self.transformation_matrix

        # Get the diagonal in the last step
        return numpy.diagonal(strain_rotated, axis1=-2, axis2=-1)

    def get_target_elastic_modulus(self) -> numpy.ndarray:
        '''The elastic modulus need to be calculated from the difference in
        fictitious strain energy of the rotated coordinate system and the known
        terms in the original coordinate.
        '''
        i, j, k, l = self.key.standard

        strain_energy_difference = \
            self.fictitious_strain_energy_rotated - self.fictitious_strain_energy

        return 2 * strain_energy_difference / (
            self.fictitious_strain[i-1,j-1] * self.fictitious_strain[k-1,l-1]
        ) / self.key.multiplicity
    
    def get_modulus_keys(self):
        return get_fictitious_strain_energy_keys(self.fictitious_strain, self.key)

    def get_modulus_keys_rotated(self):
        return get_fictitious_strain_energy_keys(self.fictitious_strain_rotated)
    
    def get_elastic_modulus(self, key: C_) -> numpy.ndarray:
        return self.modulus[key]

    def get_elastic_modulus_rotated(self, key: C_) -> numpy.ndarray:
        return self.modulus_rotated[key]

    @LazyProperty
    def value_isothermal(self):
        return self.get_target_elastic_modulus()

    @property
    def value_adiabatic(self):
        return self.value_isothermal
