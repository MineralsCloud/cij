from qha.fitting import polynomial_least_square_fitting
from qha.grid_interpolation import calculate_eulerian_strain
import numpy
from numpy import newaxis as nax
from cij.util import units, C_, E_, c_, e_, _to_gpa, _from_gpa
import itertools
import logging
import scipy.interpolate
from lazy_property import LazyProperty
from typing import Dict, Union

from .phonon_contribution import (
    LongitudinalElasticModulusPhononContribution,
    OffDiagonalElasticModulusPhononContribution,
    ShearElasticModulusPhononContribution
)

logger = logging.getLogger(__name__)

class ElasticModulusWorker:

    def __init__(self, calculator):
        self.calculator = calculator
        self.elast_data = self.calculator.elast_data
        self.pool = {
            LongitudinalElasticModulusPhononContribution: [],
            OffDiagonalElasticModulusPhononContribution: [],
        }

    @property
    def volumes(self) -> numpy.ndarray:
        '''The array of original volume points :math:`V` from input.
        '''
        return numpy.array([volume.volume for volume in self.elast_data.volumes])

    @property
    def v_array(self) -> numpy.ndarray:
        '''The array of interpolated volume points :math:`V`.
        '''
        return self.calculator.v_array

    def fit_modulus(self, moduli: numpy.ndarray, order: int = 2) -> numpy.ndarray:
        '''Interpolate static elastic constants :math:`c^\\text{st}_{ij}(V)` as
        a function of volume with polynomial least square fitting

        :param moduli: The static elastic moduli :math:`c^\\text{st}_{ij}(V)`
            to be fitted
        :param order: The order for least square fitting

        :returns: The interpolated modulus :math:`c^\\text{st}_{ij}(V)`
        '''
        # c_of_v = scipy.interpolate.UnivariateSpline(
        #     numpy.flip(self.volumes, axis=0),
        #     numpy.flip(moduli, axis=0)
        # )
        # return c_of_v(self.v_array)

        strains = calculate_eulerian_strain(self.volumes[0], self.volumes)
        strain_array = calculate_eulerian_strain(self.volumes[0], self.v_array)
        _, modulus_array = polynomial_least_square_fitting(
            strains, moduli, strain_array,
            order=order
        )
        return modulus_array
    
    def get_static_modulus(self, key: C_):
        '''Get interpolated static elastic modulus :math:`c^\\text{st}_{ij}(V)`
        as a function of volume.

        :param key: the subscript symbol for static modulus
        '''
        static_moduli = numpy.array([
            volume.static_elastic_modulus[key] for volume in self.elast_data.volumes
        ])
        static_moduli = _from_gpa(static_moduli)
        static_modulus_array = self.fit_modulus(static_moduli)
        return static_modulus_array

    
    def find_or_create_elastic_modulus_phonon_contribution(self, cls, e):
        try:
            m = next(
                m for m in  self.pool[cls]
                if numpy.allclose(list(sorted(e)), m.e)
            )
        except StopIteration:
            m = cls(self.calculator, e)
            self.pool[cls].append(m)
        return m

    def get_phonon_modulus_isothermal(self, key: C_):
        if key.is_longitudinal:
            modulus = self.find_or_create_elastic_modulus_phonon_contribution(
                LongitudinalElasticModulusPhononContribution, (1/3, 1/3))
            return modulus.value_isothermal
        elif key.is_off_diagonal:
            modulus = self.find_or_create_elastic_modulus_phonon_contribution(
                OffDiagonalElasticModulusPhononContribution, (1/3, 1/3))
            return modulus.value_isothermal
        else:
            value = self.get_shear_phonon_modulus(key, "isothermal")
            return value

    def get_phonon_modulus_adiabatic(self, key: C_):
        if key.is_longitudinal:
            modulus = self.find_or_create_elastic_modulus_phonon_contribution(
                LongitudinalElasticModulusPhononContribution, (1/3, 1/3))
            return modulus.value_adiabatic
        elif key.is_off_diagonal:
            modulus = self.find_or_create_elastic_modulus_phonon_contribution(
                OffDiagonalElasticModulusPhononContribution, (1/3, 1/3))
            return modulus.value_adiabatic
        else:
            value = self.get_shear_phonon_modulus(key, "adiabatic")
            return value

    def get_modulus_adiabatic(self, key: C_):
        if self.calculator.config["settings"]["disable_phonon_contribution"]:
            return numpy.broadcast_to(self.get_static_modulus(key)[nax,:], (self.calculator.qha_calculator.volume_base.pressures.shape))
        return self.get_static_modulus(key)[nax,:] + self.get_phonon_modulus_adiabatic(key)

    def get_modulus_isothermal(self, key: C_):
        if self.calculator.config["settings"]["disable_phonon_contribution"]:
            return numpy.broadcast_to(self.get_static_modulus(key)[nax,:], (self.calculator.qha_calculator.volume_base.pressures.shape))
        return self.get_static_modulus(key)[nax,:] + self.get_phonon_modulus_isothermal(key)

    def get_shear_phonon_modulus(self, key: C_, flag: str):

        if self.calculator.config["settings"]["disable_phonon_contribution"]:
            return numpy.zeros(self.calculator.qha_calculator.volume_base.pressures.shape)

        logger.debug(f"Start calculating {key}.")

        dims = (
            self.calculator.qha_calculator.t_array.shape[0],
            self.calculator.qha_calculator.v_array.shape[0]
        )

        if key.i != key.j: raise NotImplementedError()

        e = numpy.zeros((3, 3))

        e[key.i[0] - 1, key.i[1] - 1] = 1
        e[key.i[1] - 1, key.i[0] - 1] = 1

        eig, t = numpy.linalg.eig(e)

        eig_new = numpy.diag(t @ numpy.eye(3) @ t.T)

        nz, = numpy.nonzero(eig)

        rhs = numpy.zeros(dims)
        
        for i, j in itertools.product(nz, nz):
            
            if i > j: continue

            _key: C_ = c_(i+1, j+1)
            #if not in allowed_keys:
            if _key.is_longitudinal:
                modulus = self.find_or_create_elastic_modulus_phonon_contribution(
                    LongitudinalElasticModulusPhononContribution, 
                    (eig_new[i] / 3, eig_new[j] / 3))
            elif _key.is_off_diagonal:
                modulus = self.find_or_create_elastic_modulus_phonon_contribution(
                    OffDiagonalElasticModulusPhononContribution, 
                    (eig_new[i] / 3, eig_new[j] / 3))
            #else:
                #c_shear()

            #value = getattr(modulus, f"value_{flag}")
            value = modulus.value_isothermal

            logger.debug(f"{_key.multiplicity * eig[i] * eig[j]} * c_{_key}")

            rhs += _key.multiplicity * eig[i] * eig[j] * value

        logger.debug(f"Dividing by {key.multiplicity}")

        return (rhs - 0) / key.multiplicity

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


class ShearPhononCalculator:

    def __init__(self, key: C_, e: tuple):
        self.key = key

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
        return numpy.diag(numpy.linalg.eig(self.fictitious_strain)[0])

    @LazyProperty
    def transformation_matrix(self):
        '''The transformation matrix
        '''
        return numpy.linalg.eig(self.fictitious_strain)[1]

    @property
    def fictitious_strain_energy(self) -> numpy.ndarray:
        '''The strain energy for the fictitious strain under the original
        coordinate system except for the unknown
        '''
        resolve_elastic_modulus = ...
        calculate_fictitious_strain_energy(
            self.fictitious_strain,
            resolve_elastic_modulus,
            self.key,
        )

    @property
    def fictitious_strain_energy_rotated(self) -> numpy.ndarray:
        '''The strain energy for the fictitious strain under the rotated
        coordinate system
        '''
        resolve_elastic_modulus = ...
        calculate_fictitious_strain_energy(
            numpy.diag(numpy.diag(self.fictitious_strain_rotated)),
            resolve_elastic_modulus
        )

    def get_target_elastic_modulus(self) -> numpy.ndarray:
        '''The elastic modulus need to be calculated from the difference in
        fictitious strain energy of the rotated coordinate system and the known
        terms in the original coordinate.
        '''
        i, j, k, l = self.key.standard
        strain_energy_difference = self.fictitious_strain_energy_rotated - self.fictitious_strain_energy
        return 2 * strain_energy_difference / (
            self.fictitious_strain[i,j] * self.fictitious_strain[k,l]
        ) / self.key.multiplicity
    
    def get_elastic_modulus(self, key: C_) -> numpy.ndarray:
        pass

    def get_elastic_modulus_rotated(self, key: C_) -> numpy.ndarray:
        pass