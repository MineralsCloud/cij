import itertools
import logging
import itertools
from typing import Dict, Union, List

import numpy
from numpy import newaxis as nax
import scipy.interpolate
from lazy_property import LazyProperty

from qha.fitting import polynomial_least_square_fitting
from qha.grid_interpolation import calculate_eulerian_strain

from cij.util import units, C_, E_, c_, e_, _to_gpa, _from_gpa

from .tasks import PhononContributionTaskList


logger = logging.getLogger(__name__)

class FullThermalElasticModulus:

    def __init__(self, calculator: 'cij.core.Calculator'):
        self.calculator = calculator
        self.elast_data = self.calculator.elast_data
        self.calculate_phonon_contribution()
    
    @property
    def modulus_keys(self) -> List[C_]:
        return self.calculator.modulus_keys

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
    
    def calculate_phonon_contribution(self):
        self._phonon_contribution_task_list = PhononContributionTaskList(self.calculator)
        self._phonon_contribution_task_list.resolve((1/3, 1/3, 1/3), self.modulus_keys)
        self._phonon_contribution_task_list.calculate()
        self._adiabatic_phonon_contribution = self._phonon_contribution_task_list.get_adiabatic_results()
        self._isothermal_phonon_contribution = self._phonon_contribution_task_list.get_isothermal_results()

    @LazyProperty
    def modulus_adiabatic(self) -> None:
        results = dict()
        for key in self.modulus_keys:
            results[key] = self.get_static_modulus(key)[nax,:] + self._adiabatic_phonon_contribution[key]
        return results

    @LazyProperty
    def modulus_isothermal(self) -> None:
        results = dict()
        for key in self.modulus_keys:
            results[key] = self.get_static_modulus(key)[nax,:] + self._isothermal_phonon_contribution[key]
        return results

