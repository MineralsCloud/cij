from cij.core.mode_gamma_interpolate import interpolate_modes
import cij.io
from cij.core.phonon_modulus import LogitudinalElasticModulusPhononContribution
from cij.core.qha_adapter import QHACalculatorAdapter
from qha.v2p import v2p
from cij.util import c_, units
from cij.core.modulus_worker import ElasticModulusWorker

from lazy_property import LazyProperty
import re
import numpy

REGEX_CIJ = r'^c_?([1-6]{2,2}|[1-3]{4,4})(s|t)?$'

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

        self._process_cij()
        self._calc_velocities()

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
            self.modulus_isothermal[key] = self.modulus_worker.get_modulus_adiabatic(key)
            #LogitudinalElasticModulusPhononContribution(self, (1/3, 1/3))

    def _calc_velocities(self):

        def au_to_gpa(p):
            return units.Quantity(
                p,
                units.rydberg / units.bohr ** 3).to(units.GPa).magnitude

        t_amount = self.qha_calculator.t_array.shape[0]
        v_amount = self.qha_calculator.v_array.shape[0]

        s11 = numpy.zeros((t_amount, v_amount))
        s22 = numpy.zeros((t_amount, v_amount))
        s33 = numpy.zeros((t_amount, v_amount))
        s44 = numpy.zeros((t_amount, v_amount))
        s55 = numpy.zeros((t_amount, v_amount))
        s66 = numpy.zeros((t_amount, v_amount))
        s12 = numpy.zeros((t_amount, v_amount))
        s13 = numpy.zeros((t_amount, v_amount))
        s23 = numpy.zeros((t_amount, v_amount))

        for i in range(t_amount):
            for j in range(v_amount):
                c11, c22, c33, c44, c55, c66, c12, c13, c23 = (
                    au_to_gpa(self.volume_base.c11s[i,j]),
                    au_to_gpa(self.volume_base.c22s[i,j]),
                    au_to_gpa(self.volume_base.c33s[i,j]),
                    au_to_gpa(self.volume_base.c44s[i,j]),
                    au_to_gpa(self.volume_base.c55s[i,j]),
                    au_to_gpa(self.volume_base.c66s[i,j]),
                    au_to_gpa(self.volume_base.c12s[i,j]),
                    au_to_gpa(self.volume_base.c13s[i,j]),
                    au_to_gpa(self.volume_base.c23s[i,j]),
                )

                cij = numpy.array([
                    [c11, c12, c13,   0,   0,   0],
                    [c12, c22, c23,   0,   0,   0],
                    [c13, c23, c33,   0,   0,   0],
                    [  0,   0,   0, c44,   0,   0],
                    [  0,   0,   0,   0, c55,   0],
                    [  0,   0,   0,   0,   0, c66],
                ])
                sij = numpy.linalg.inv(cij)

                s11[i,j] = sij[0,0]
                s22[i,j] = sij[1,1]
                s33[i,j] = sij[2,2]
                s44[i,j] = sij[3,3]
                s55[i,j] = sij[4,4]
                s66[i,j] = sij[5,5]
                s12[i,j] = sij[1,2]
                s13[i,j] = sij[1,3]
                s23[i,j] = sij[2,3]



        c11 = au_to_gpa(self.volume_base.c11s)
        c22 = au_to_gpa(self.volume_base.c22s)
        c33 = au_to_gpa(self.volume_base.c33s)
        c44 = au_to_gpa(self.volume_base.c44s)
        c55 = au_to_gpa(self.volume_base.c66s)
        c66 = au_to_gpa(self.volume_base.c66s)
        c12 = au_to_gpa(self.volume_base.c12s)
        c23 = au_to_gpa(self.volume_base.c23s)
        c13 = au_to_gpa(self.volume_base.c13s)

        Bm_V = (c11 + c22 + c33 + 2 * (c12 + c23 + c13)) / 9
        Bm_R = 1 / (s11 + s22 + s33 + 2 * (s12 + s23 + s13))
        Bm_VRH = (Bm_V + Bm_R) / 2

        G_V = ((c11 + c22 + c33) - (c12 + c23 + c13) + 3 * (c44 + c55 + c66)) / 15
        G_R = 15 / (4 * (s11 + s22 + s33) - 4 * (s12 + s23 + s13) + 3 * (s44 + s55 + s66))
        G_VRH = (G_V + G_R) / 2

        import scipy.constants

        m = self.elast_data.cellmass 
        v = units.Quantity(self.qha_calculator.v_array[numpy.newaxis, :], units.bohr**3).to(units.angstrom**3).magnitude
        rho = m / v
        N = scipy.constants.physical_constants["Avogadro constant"][0]

        self._primary_velocities = numpy.sqrt((Bm_VRH + 4 / 3 * G_VRH) / rho * 1e-18 * N) * 1e-3
        self._secondary_velocities = numpy.sqrt(G_VRH / rho * 1e-18 * N) * 1e-3

        #self._born = {}
        #self._born[1] = c11 * c22 - c12 ** 2
        #self._born[2] = c11 * c22 * c33 + 2 * c12 * c13 * c23 - c11 * c23 **2 - c22 * c13 ** 2 - c33 * c12 ** 2
    
    @LazyProperty
    def volume_base(self):
        return self.volume_based_result

    @LazyProperty
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
            key = c_(res.group(1))
            if key not in self.calculator.modulus_keys:
                raise AttributeError()
            if res.group(2) == 't':
                return self.calculator.modulus_isothermal[key]
            else:
                return self.calculator.modulus_adiabatic[key]
        raise AttributeError()

    @property
    def primary_velocities(self):
        return self.calculator._primary_velocities

    @property
    def secondary_velocities(self):
        return self.calculator._secondary_velocities