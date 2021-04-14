import re
import itertools
from typing import List, Tuple, Union, Iterable
from pathlib import Path
import numpy
import scipy.constants
from lazy_property import LazyProperty
from collections import UserDict

from qha.v2p import v2p
from qha.fitting import polynomial_least_square_fitting
from qha.grid_interpolation import calculate_eulerian_strain

import cij.io
from cij.util import c_, C_, units, _to_gpa, _to_ang3
from cij.io.traditional.qha_output import save_x_tv, save_x_tp
from cij.io.traditional.elast_dat import apply_symetry_on_elast_data
from cij.io.output import ResultsWriter

from .mode_gamma import interpolate_modes
from .qha_adapter import QHACalculatorAdapter
# from .modulus_worker import ElasticModulusWorker
from .full_modulus import FullThermalElasticModulus

import logging

logger = logging.getLogger(__name__)

REGEX_CIJ = r'^(c|s)_?([1-6]{2,2}|[1-3]{4,4})(s|t)?$'

class Calculator:
    '''The main entrance for QHA calculator

    :param config_fname: the location of the configuration file
    '''

    def __init__(self, config_fname: str):
        self._load(config_fname)
        self._apply_elastic_constants_symmetry()
        self._interpolate_modes()

        self.nv = self.qha_input.nv
        self.np = self.qha_input.np
        self.nq = self.qha_input.nq
        self.na = self.qha_input.na

        self.volume_based_result = CijVolumeBaseInterface(self)
        self.pressure_based_result = CijPressureBaseInterface(self)

        self._calculate_pressure_static()
        self._process_cij()
        self._calculate_compliances()
        # self._calc_velocities()

    def _load(self, config_fname: str):

        config_fname = Path(config_fname)
        work_dir = config_fname.parent

        self.config = cij.io.read_config(config_fname)
        self.config = cij.io.apply_default_config(self.config)
        self.qha_input = cij.io.traditional.read_energy(work_dir / self.config["qha"]["input"])
        self.elast_data = cij.io.traditional.read_elast_data(work_dir / self.config["elast"]["input"])
        self.qha_calculator = QHACalculatorAdapter(
            self.config["qha"]["settings"],
            self.qha_input
        )
    
    def _apply_elastic_constants_symmetry(self):

        symmetry = self.config["elast"]["settings"]["symmetry"]
        system = symmetry.get("system", None)

        if system == None or system == "triclinic":
            logger.warning(f"Symmetry constraints check not performed! Make sure to fill in all non-zero terms for correct VRH averages!")
        else:
            apply_symetry_on_elast_data(self.elast_data, symmetry)

    def _interpolate_modes(self):
        interp_freq, gamma_i, vdr_dv = interpolate_modes(
            self.qha_input, self.qha_calculator.v_array,
            method=self.config["elast"]["settings"]["mode_gamma"]["interpolator"],
            order=self.config["elast"]["settings"]["mode_gamma"]["order"]
        )
        self.freq_array = interp_freq
        self.mode_gamma = [vdr_dv, gamma_i, gamma_i**2]
    
    @LazyProperty
    def modulus_keys(self) -> List[C_]:
        '''Elastic coefficient keys
        '''
        return list(self.elast_data.volumes[0].static_elastic_modulus.keys())

    # def _process_cij(self):
    #     self.modulus_adiabatic = {}
    #     self.modulus_isothermal = {}
    #     self.modulus_worker = ElasticModulusWorker(self)
    #     for key in self.modulus_keys:
    #         self.modulus_adiabatic[key] = self.modulus_worker.get_modulus_adiabatic(key)
    #         self.modulus_isothermal[key] = self.modulus_worker.get_modulus_isothermal(key)

    def _process_cij(self):
        self._full_modulus = FullThermalElasticModulus(self)
        self.modulus_adiabatic = self._full_modulus.modulus_adiabatic
        self.modulus_isothermal = self._full_modulus.modulus_isothermal

    def _calculate_pressure_static(self, order: int = 3):

        volumes = numpy.array([volume.volume for volume in self.qha_input.volumes])
        static_energies = numpy.array([volume.energy for volume in self.qha_input.volumes])

        strains = calculate_eulerian_strain(volumes[0], volumes)
        strain_array = calculate_eulerian_strain(volumes[0], self.v_array)
        _, static_energy_array = polynomial_least_square_fitting(
            strains, static_energies, strain_array,
            order=order
        )
        
        self.static_p_array = - numpy.gradient(static_energy_array) / numpy.gradient(self.v_array)
    
    @property
    def dims(self) -> Tuple[int, int]:
        '''The dimension of the temperature-volume :math:`(T, V)` grid
        '''
        nt = self.qha_calculator.t_array.shape[0]
        ntv = self.qha_calculator.v_array.shape[0]

        return (nt, ntv)
    
    def _calculate_compliances(self):

        elastic_moduli = numpy.zeros((*self.dims, 6, 6))
        self._compliances = {}
        for key in self.modulus_keys:
            for i, j in set(itertools.permutations(key.voigt, 2)):
                elastic_moduli[:, :, i-1, j-1] = self.modulus_adiabatic[key]
        
        compliances = numpy.linalg.inv(elastic_moduli)

        for i, j in itertools.product(range(6), range(6)):
            if i > j: continue
            if numpy.allclose(compliances[:, :, i, j], 0): continue
            self._compliances[c_(i+1, j+1)] = compliances[:, :, i, j]
        
    @property
    def volume_base(self) -> 'CijVolumeBaseInterface':
        return self.volume_based_result

    @property
    def pressure_base(self) -> 'CijPressureBaseInterface':
        return self.pressure_based_result
 
    def __getattr__(self, prop):
        return getattr(self.qha_calculator, prop)
    
    def write_output(self):

        output_config = self.config["output"]

        if "pressure_base" in output_config.keys():
            self.pressure_base.write_variables(output_config["pressure_base"])

        if "volume_base" in output_config.keys():
            self.volume_base.write_variables(output_config["volume_base"])


class CijVolumeBaseInterface:
    '''Elastic and accoustic properties calculated at the volume-temperature
    :math:`(T, V)` grid.
    '''

    _base_name = "tv"

    def __init__(self, calculator: Calculator):
        self.calculator = calculator

    @property
    def v_array(self) -> numpy.ndarray:
        '''The array of volume points :math:`V` of the temperature-volume
        :math:`(T, V)` grid
        '''
        return self.calculator.qha_calculator.volume_base.v_array

    @property
    def t_array(self) -> numpy.ndarray:
        '''The array of temperature points :math:`T` of the temperature-volume
        :math:`(T, V)` grid
        '''
        return self.calculator.qha_calculator.volume_base.t_array

    def __getattr__(self, name):
        res = re.search(REGEX_CIJ, name)
        if res:
            if res.group(1) == "c":
                key = c_(res.group(2))
                if key not in self.calculator.modulus_keys:
                    raise AttributeError()
                if res.group(3) == 't':
                    return self.calculator.modulus_isothermal[key]
                else:
                    return self.calculator.modulus_adiabatic[key]
            if res.group(1) == "s":
                key = c_(res.group(2))
                if key not in self.calculator._compliances.keys():
                    raise AttributeError()
                if res.group(1) == 't':
                    raise AttributeError()
                return self.calculator._compliances[key]
        raise AttributeError(name)

    @property
    def modulus_adiabatic(self) -> dict:
        return self.calculator.modulus_adiabatic

    @property
    def modulus_isothermal(self) -> dict:
        return self.calculator.modulus_isothermal

    @property
    def bulk_modulus_voigt(self) -> numpy.ndarray:
        '''Voigt average of bulk modulus :math:`K_\\text{V}(T, V)` as a function
        of temperature and volume.

        .. math::
            K_\\text{V} = [(c_{11}+c_{22}+c_{33}) + 2(c_{12}+c_{23}+c_{31})]/9
        '''
        return (self.c11 + self.c22 + self.c33 \
                + 2 * (self.c12 + self.c23 + self.c13)) / 9

    @property
    def bulk_modulus_reuss(self) -> numpy.ndarray:
        '''Reuss average of bulk modulus :math:`K_\\text{R}(T, V)` as a function
        of temperature and volume.

        .. math::
            K_\\text{R} = [(s_{11}+s_{22}+s_{33})+2(s_{12}+s_{23}+s_{31})]^{-1}
        '''
        return 1 / (self.s11 + self.s22 + self.s33 \
                    + 2 * (self.s12 + self.s23 + self.s13))
    
    @property
    def bulk_modulus_voigt_reuss_hill(self) -> numpy.ndarray:
        '''Voigt-Reuss-Hill average of bulk modulus :math:`K_\\text{VRH}(T, V)`
        as a function of temperature and volume.

        .. math::
            K_\\text{VRH} = (K_\\text{V} + K_\\text{R}) / 2
        '''
        return (self.bulk_modulus_reuss + self.bulk_modulus_voigt) / 2
    
    @property
    def shear_modulus_voigt(self) -> numpy.ndarray:
        '''Voigt average of shear modulus :math:`G_\\text{V}(T, V)` as a function
        of temperature and volume.

        .. math::
            G_\\text{V} = [
                    (c_{11} + c_{22} + c_{33})
                -   (c_{12} + c_{23} + c_{31})
                + 3 (c_{44} + c_{55} + c_{66})
            ] / 15
        '''
        return ( \
            + (self.c11 + self.c22 + self.c33) \
            - (self.c12 + self.c23 + self.c13)
            + 3 * (self.c44 + self.c55 + self.c66)) / 15

    @property
    def shear_modulus_reuss(self) -> numpy.ndarray:
        '''Reuss average of shear modulus :math:`K_\\text{R}(T, V)` as a function
        of temperature and volume.

        .. math::
            G_\\text{R} = 15 / [
                  4 (s_{11} + s_{22} + s_{33})
                + 2 (s_{12} + s_{23} + s_{31})
                + 3 (s_{44} + s_{55} + s_{66})
            ]
        ''' 
        return 15 / ( \
            + 4 * (self.s11 + self.s22 + self.s33) \
            - 4 * (self.s12 + self.s23 + self.s13) \
            + 3 * (self.s44 + self.s55 + self.s66))

    @property
    def shear_modulus_voigt_reuss_hill(self) -> numpy.ndarray:
        '''Voigt-Reuss-Hill average of shear modulus :math:`G_\\text{VRH}(T, V)`
        as a function of temperature and volume.

        .. math::
            G_\\text{VRH} = (G_\\text{V} + G_\\text{R}) / 2
        '''
        return (self.shear_modulus_reuss + self.shear_modulus_voigt) / 2
        
    @property
    def mass(self) -> float:
        '''The mass per cell :math:`m` in kilogram.'''
        m = self.calculator.elast_data.cellmass 
        N = scipy.constants.physical_constants["Avogadro constant"][0]
        return m * 1e-3 / N
    
    @property
    def primary_velocities(self) -> numpy.ndarray:
        '''Primary accoustic wave velocity :math:`v_\\text{p}(T, V)` as a
        function of temperature and volume.

        .. math::
            v_\\text{p} = \\sqrt{\\frac{K_\\text{VRH} + 3/4 \, G_\\text{VRH} }{\\rho}}
        '''
        e = units.Quantity((self.bulk_modulus_voigt_reuss_hill + 4 / 3 * self.shear_modulus_voigt_reuss_hill) * self.v_array, units.rydberg).to(units.kg * units.km ** 2 / units.s ** 2).magnitude
        return numpy.sqrt(e / self.mass)

    @property
    def secondary_velocities(self) -> numpy.ndarray:
        '''Secondary accoustic wave velocity :math:`v_\\text{s}(T, V)` as a
        function of temperature and volume.

        .. math::
            v_\\text{s} = \\sqrt{\\frac{G_\\text{VRH}}{\\rho}}
        '''
        e = units.Quantity(self.shear_modulus_voigt_reuss_hill * self.v_array, units.rydberg).to(units.kg * units.km ** 2 / units.s ** 2).magnitude
        return numpy.sqrt(e / self.mass)
    
    @property
    def pressures(self) -> numpy.ndarray:
        return self.calculator.qha_calculator.volume_base.pressures

    def write_table(self, fname: str, value: numpy.ndarray) -> None:
        '''Write variable as functions of temperature and volume in QHA
        output format, intended to be called by ``ResultsWriter`` only.
        '''
        v_array = _to_ang3(self.v_array)
        save_x_tv(value, self.t_array, v_array, self.t_array, fname)

    def write_variables(self, variables: Iterable[Union[str, dict]]):
        '''Write variables to files

        :param variables: List of varables to be written to file, see output
            file description for a detailed description
        '''
        writer = ResultsWriter(self)
        for c in variables:
            writer.write(c)

class CijPressureBaseModulusInterface:

    def __init__(self, modulus, v2p: callable):
        self.modulus = modulus
        self.v2p = v2p

    def items(self):
        for key in self.modulus.keys():
            yield key, self[key]

    def __getitem__(self, key: str) -> numpy.ndarray:
        return self.v2p(self.modulus[key])

    
class CijPressureBaseInterface:
    '''Elastic and accoustic properties calculated at the temperature-pressure
    :math:`(T, P)` grid.
    '''

    _base_name = "tp"

    def __init__(self, calculator: Calculator):
        self.calculator = calculator

    def v2p(self, func_of_t_v: numpy.ndarray) -> numpy.ndarray:
        '''The conversion function from :math:`(T, V)` to :math:`(T, P)` grid

        .. math::
            f(T, V) \\rightarrow f(T, P)

        :param func_of_t_v: the input function :math:`f(T, V)` under the :math:`(T, V)` grid
        :returns: the output function :math:`f(T, P)` under the :math:`(T, P)` grid
        '''
        return v2p(func_of_t_v, self.calculator.qha_calculator.volume_base.pressures, self.p_array)
    
    @property
    def p_array(self) -> numpy.ndarray:
        '''The array of pressure points :math:`P` of the pressure-temperature
        :math:`(T, P)` grid
        '''
        return self.calculator.qha_calculator.pressure_base.p_array

    @property
    def t_array(self) -> numpy.ndarray:
        '''The array of temperature points :math:`T` of the pressure-temperature
        :math:`(T, P)` grid
        '''
        return self.calculator.qha_calculator.pressure_base.t_array

    @property
    def modulus_adiabatic(self) -> CijPressureBaseModulusInterface:
        '''Adiabatic elastic modulus :math:`c^S_{ij}(T, P)`  as a function of
        temperature and pressure.
        '''
        return CijPressureBaseModulusInterface(
            self.calculator.modulus_adiabatic,
            self.v2p
        )

    @property
    def modulus_isothermal(self) -> CijPressureBaseModulusInterface:
        '''Isothermal elastic modulus :math:`c^T_{ij}(T, P)`  as a function of
        temperature and pressure.
        '''
        return CijPressureBaseModulusInterface(
            self.calculator.modulus_isothermal,
            self.v2p
        )

    @property
    def bulk_modulus_voigt(self) -> numpy.ndarray:
        '''Voigt average of bulk modulus :math:`K_\\text{V}(T, P)` as a function
        of temperature and pressure.

        .. math::
            K_\\text{V} = [(c_{11}+c_{22}+c_{33}) + 2(c_{12}+c_{23}+c_{31})]/9
        '''
        return self.v2p(self.calculator.volume_base.bulk_modulus_voigt)

    @property
    def bulk_modulus_reuss(self) -> numpy.ndarray:
        '''Reuss average of bulk modulus :math:`K_\\text{R}(T, P)` as a function
        of temperature and pressure.

        .. math::
            K_\\text{R} = [(s_{11}+s_{22}+s_{33})+2(s_{12}+s_{23}+s_{31})]^{-1}
        '''
        return self.v2p(self.calculator.volume_base.bulk_modulus_reuss)
   
    @property
    def bulk_modulus_voigt_reuss_hill(self) -> numpy.ndarray:
        '''Voigt-Reuss-Hill average of bulk modulus :math:`K_\\text{VRH}(T, P)`
        as a function of temperature and pressure.

        .. math::
            K_\\text{VRH} = (K_\\text{V} + K_\\text{R}) / 2
        '''
        return self.v2p(self.calculator.volume_base.bulk_modulus_voigt_reuss_hill)
    
    @property
    def shear_modulus_voigt(self) -> numpy.ndarray:
        '''Voigt average of shear modulus :math:`G_\\text{V}(T, P)` as a function
        of temperature and pressure.

        .. math::
            G_\\text{V} = [
                    (c_{11} + c_{22} + c_{33})
                -   (c_{12} + c_{23} + c_{31})
                + 3 (c_{44} + c_{55} + c_{66})
            ] / 15
        '''
        return self.v2p(self.calculator.volume_base.shear_modulus_voigt)

    @property
    def shear_modulus_reuss(self) -> numpy.ndarray:
        '''Reuss average of shear modulus :math:`G_\\text{R}(T, P)` as a function
        of temperature and pressure.

        .. math::
            G_\\text{R} = 15 / [
                  4 (s_{11} + s_{22} + s_{33})
                + 2 (s_{12} + s_{23} + s_{31})
                + 3 (s_{44} + s_{55} + s_{66})
            ] 
        '''
        return self.v2p(self.calculator.volume_base.shear_modulus_reuss)

    @property
    def shear_modulus_voigt_reuss_hill(self) -> numpy.ndarray:
        '''Voigt-Reuss-Hill average of shear modulus :math:`G_\\text{VRH}(T, P)`
        as a function of temperature and pressure.

        .. math::
            G_\\text{VRH} = (G_\\text{V} + G_\\text{R}) / 2
        '''
        return self.v2p(self.calculator.volume_base.shear_modulus_voigt_reuss_hill)
        
    @property
    def mass(self) -> float:
        '''The mass per cell :math:`m` in kilogram.'''
        return self.calculator.volume_base.mass
    
    @property
    def primary_velocities(self) -> numpy.ndarray:
        '''Primary accoustic wave velocity :math:`v_\\text{p}(T, P)` as a
        function of temperature and pressure.

        .. math::
            v_\\text{p} = \\sqrt{\\frac{K_\\text{VRH} + 3/4 \, G_\\text{VRH} }{\\rho}}
        '''
        return self.v2p(self.calculator.volume_base.primary_velocities)

    @property
    def secondary_velocities(self) -> numpy.ndarray:
        '''Secondary accoustic wave velocity :math:`v_\\text{s}(T, P)` as a
        function of temperature and pressure.

        .. math::
            v_\\text{s} = \\sqrt{\\frac{G_\\text{VRH}}{\\rho}}
        '''
 
        return self.v2p(self.calculator.volume_base.secondary_velocities)

    @property
    def volumes(self) -> numpy.ndarray:
        return self.calculator.qha_calculator.pressure_base.volumes

    def __getattr__(self, name):
        func_of_t_v = getattr(self.calculator.volume_base, name)
        func_of_t_p = self.v2p(func_of_t_v)
        return func_of_t_p
    
    def write_table(self, fname: str, value: numpy.ndarray) -> None:
        '''Write variable as functions of temperature and volume in QHA
        output format, intended to be called by ``ResultsWriter`` only.
        '''
        p_array = _to_gpa(self.p_array)
        save_x_tp(value, self.t_array, p_array, p_array, fname)
    
    def write_variables(self, variables: Iterable[Union[str, dict]]):
        '''Write variables to files

        :param variables: List of varables to be written to file, see output
            file description for a detailed description
        '''
        writer = ResultsWriter(self)
        for c in variables:
            writer.write(c)
