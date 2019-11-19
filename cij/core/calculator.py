from cij.core.mode_gamma_interpolate import interpolate_modes
import cij.io
from cij.core.phonon_modulus import LogitudinalElasticModulusPhononContribution
from cij.core.qha_adapter import QHACalculatorAdapter
from qha.v2p import v2p
from cij.util import c_, units
from cij.core.modulus_worker import ElasticModulusWorker
from qha.fitting import polynomial_least_square_fitting
from qha.grid_interpolation import calculate_eulerian_strain

from lazy_property import LazyProperty
import re
import numpy
import itertools
import scipy.constants

import logging

logger = logging.Logger(__name__)

REGEX_CIJ = r'^(c|s)_?([1-6]{2,2}|[1-3]{4,4})(s|t)?$'

class Calculator:

    def __init__(self, config_fname: str):
        self._load(config_fname)
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
        self.config = cij.io.read_config(config_fname)
        self.qha_input = cij.io.traditional.read_energy(self.config["qha"]["input"])
        self.elast_data = cij.io.traditional.read_elast_data(self.config["elast"]["input"])
        self.qha_calculator = QHACalculatorAdapter(
            self.config["qha"]["settings"],
            self.qha_input
        )

    def _interpolate_modes(self):
        interp_freq, gamma_i, vdr_dv = interpolate_modes(
            self.qha_input, self.qha_calculator.v_array,
            method=self.config["settings"]["mode_gamma"]["interpolator"],
            order=self.config["settings"]["mode_gamma"]["order"]
        )
        self.freq_array = interp_freq
        self.mode_gamma = [vdr_dv, gamma_i, gamma_i**2]
    
    @LazyProperty
    def modulus_keys(self):
        return [c_(key) for key in self.config["settings"]["cij_keys"]]

    def _process_cij(self):
        self.modulus_adiabatic = {}
        self.modulus_isothermal = {}
        self.modulus_worker = ElasticModulusWorker(self)
        for key in self.modulus_keys:
            self.modulus_adiabatic[key] = self.modulus_worker.get_modulus_adiabatic(key)
            self.modulus_isothermal[key] = self.modulus_worker.get_modulus_isothermal(key)

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
    def dims(self):
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
    def volume_base(self):
        return self.volume_based_result

    @property
    def pressure_base(self):
        return self.pressure_based_result
 
    def __getattr__(self, prop):
        return getattr(self.qha_calculator, prop)
    
class CijPressureBaseInterface:
    def __init__(self, calculator: Calculator):
        self.calculator = calculator

    def v2p(self, func_of_t_v):
        return v2p(func_of_t_v, self.calculator.qha_calculator.volume_base.pressures, self.p_array)
    
    @property
    def p_array(self):
        return self.calculator.qha_calculator.pressure_base.p_array

    @property
    def t_array(self):
        return self.calculator.qha_calculator.pressure_base.t_array

    def __getattr__(self, name):
        func_of_t_v = getattr(self.calculator.volume_base, name)
        func_of_t_p = self.v2p(func_of_t_v)
        return func_of_t_p

class CijVolumeBaseInterface:
    def __init__(self, calculator: Calculator):
        self.calculator = calculator

    @property
    def v_array(self):
        return self.calculator.qha_calculator.volume_base.v_array

    @property
    def t_array(self):
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
    def bulk_modulus_voigt(self):
        return (self.c11 + self.c22 + self.c33 \
                + 2 * (self.c12 + self.c23 + self.c13)) / 9

    @property
    def bulk_modulus_reuss(self):
        return 1 / (self.s11 + self.s22 + self.s33 \
                    + 2 * (self.s12 + self.s23 + self.s13))
    
    @property
    def bulk_modulus_voigt_reuss_hill(self):
        return (self.bulk_modulus_reuss + self.bulk_modulus_voigt) / 2
    
    @property
    def shear_modulus_voigt(self):
        return ( \
            + (self.c11 + self.c22 + self.c33) \
            - (self.c12 + self.c23 + self.c13)
            + 3 * (self.c44 + self.c55 + self.c66)) / 15

    @property
    def shear_modulus_reuss(self):
        return 15 / ( \
            + 4 * (self.s11 + self.s22 + self.s33) \
            - 4 * (self.s12 + self.s23 + self.s13) \
            + 3 * (self.s44 + self.s55 + self.s66))

    @property
    def shear_modulus_voigt_reuss_hill(self):
        return (self.shear_modulus_reuss + self.shear_modulus_voigt) / 2
        
    @property
    def mass(self):
        m = self.calculator.elast_data.cellmass 
        N = scipy.constants.physical_constants["Avogadro constant"][0]
        return m * 1e-3 / N
    
    @property
    def primary_velocities(self):
        e = units.Quantity((self.bulk_modulus_voigt_reuss_hill + 4 / 3 * self.shear_modulus_voigt_reuss_hill) * self.v_array, units.rydberg).to(units.kg * units.km ** 2 / units.s ** 2).magnitude
        #print(e, self.mass)
        return numpy.sqrt(e / self.mass)

    @property
    def secondary_velocities(self):
        e = units.Quantity(self.shear_modulus_voigt_reuss_hill * self.v_array, units.rydberg).to(units.kg * units.km ** 2 / units.s ** 2).magnitude
        return numpy.sqrt(e / self.mass)
