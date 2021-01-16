import itertools
import logging
import itertools
from typing import Dict, Union, List

import numpy
from numpy import newaxis as nax
import scipy.interpolate
from lazy_property import LazyProperty

# from qha.fitting import polynomial_least_square_fitting
from qha.grid_interpolation import calculate_eulerian_strain

from cij.util import units, C_, E_, c_, e_, _to_gpa, _from_gpa

from .tasks import PhononContributionTaskList


logger = logging.getLogger(__name__)

class FullThermalElasticModulus:
    '''
    Calculates the static and thermal part of the elastic constants.
    '''

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

        strains = calculate_eulerian_strain(self.volumes[0], self.volumes)
        strain_array = calculate_eulerian_strain(self.volumes[0], self.v_array)
        p = numpy.polyfit(strains, self.volumes * moduli, deg = order + 1)
        modulus_array = numpy.polyval(p, strain_array) / self.v_array
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
    
    def _get_init_strain(self) -> tuple:
        if "init_strain" in self.calculator.config["elast"]["settings"]:
            strain = self.calculator.config["elast"]["settings"]["init_strain"]
            _sum = sum(strain)
            return tuple(x / _sum for x in strain)
        return (1/3, 1/3, 1/3)
    
    def get_axial_strains(self) -> numpy.ndarray:

        # TODO: need to write a clearer implementation

        ntv = self.v_array.shape[0]

        if len(self.elast_data.lattice_parmeters) == 0:
            return numpy.ones((ntv, 3))

        lattice_params = numpy.array(self.elast_data.lattice_parmeters)
        strains = numpy.zeros((ntv, 3))

        for i in range(3):
            params = self.fit_modulus(lattice_params[:, i])
            tmp = params[[0, *range(len(params)), -1]]
            strains[:,i] = (tmp[2:] - tmp[:-2]) / (tmp[2:] + tmp[:-2])
        
        # strains = strains / strains[0,:] # TODO: not really need this line
        strains = strains / numpy.sum(strains, axis=1, keepdims=True)

        return strains
    
    def calculate_phonon_contribution(self) -> None:

        logging.debug("init_strain -> " + repr(self._get_init_strain()))
        axial_strains = self.get_axial_strains()
        # print(axial_strains)

        self._phonon_contribution_task_list = PhononContributionTaskList(self.calculator)
        # self._phonon_contribution_task_list.resolve(self._get_init_strain(), self.modulus_keys)
        self._phonon_contribution_task_list.resolve(axial_strains, self.modulus_keys)
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

