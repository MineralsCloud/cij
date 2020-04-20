from qha.fitting import polynomial_least_square_fitting
from qha.grid_interpolation import calculate_eulerian_strain
import numpy
from numpy import newaxis as nax
from cij.util import units, C_, E_, c_, e_, _to_gpa, _from_gpa
import itertools
import logging
import scipy.interpolate

from .phonon_contribution import (
    LogitudinalElasticModulusPhononContribution,
    OffDiagnonalElasticModulusPhononContribution,
    ShearElasticModulusPhononContribution
)

logger = logging.getLogger(__name__)

class ElasticModulusWorker:
    def __init__(self, calculator):
        self.calculator = calculator
        self.elast_data = self.calculator.elast_data
        self.pool = {
            LogitudinalElasticModulusPhononContribution: [],
            OffDiagnonalElasticModulusPhononContribution: [],
        }

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
    
    def get_static_modulus(self, key: tuple):
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
                LogitudinalElasticModulusPhononContribution, (1/3, 1/3))
            return modulus.value_isothermal
        elif key.is_off_diagonal:
            modulus = self.find_or_create_elastic_modulus_phonon_contribution(
                OffDiagnonalElasticModulusPhononContribution, (1/3, 1/3))
            return modulus.value_isothermal
        else:
            value = self.get_shear_phonon_modulus(key, "isothermal")
            return value

    def get_phonon_modulus_adiabatic(self, key: C_):
        if key.is_longitudinal:
            modulus = self.find_or_create_elastic_modulus_phonon_contribution(
                LogitudinalElasticModulusPhononContribution, (1/3, 1/3))
            return modulus.value_adiabatic
        elif key.is_off_diagonal:
            modulus = self.find_or_create_elastic_modulus_phonon_contribution(
                OffDiagnonalElasticModulusPhononContribution, (1/3, 1/3))
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
                    LogitudinalElasticModulusPhononContribution, 
                    (eig_new[i] / 3, eig_new[j] / 3))
            elif _key.is_off_diagonal:
                modulus = self.find_or_create_elastic_modulus_phonon_contribution(
                    OffDiagnonalElasticModulusPhononContribution, 
                    (eig_new[i] / 3, eig_new[j] / 3))
            #else:
                #c_shear()

            #value = getattr(modulus, f"value_{flag}")
            value = modulus.value_isothermal

            logger.debug(f"{_key.multiplicity * eig[i] * eig[j]} * c_{_key}")

            rhs += _key.multiplicity * eig[i] * eig[j] * value

        logger.debug(f"Dividing by {key.multiplicity}")

        return (rhs - 0) / key.multiplicity
