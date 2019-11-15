from typing import Tuple, List
import numpy
from numpy import newaxis as nax
import scipy.constants
from lazy_property import LazyProperty
import logging
from cij.util import units

logger = logging.getLogger(__name__)

_h = scipy.constants.physical_constants["molar Planck constant times c"][0] \
    / scipy.constants.physical_constants["Avogadro constant"][0]
_k = scipy.constants.physical_constants["Boltzmann constant in eV/K"][0]

h_div_k =  units.Quantity(_h / _k,  units.J * units.m / units.eV * units.K).to(units.cm * units.K).magnitude


class ElasticModulus:
    pass

def average_over_modes(amount: numpy.ndarray, q_weights: numpy.ndarray):
    dims = len(amount.shape)
    _amount = amount.copy()
    clear_gamma_point(_amount)

    return numpy.average(
        numpy.average(_amount, axis=dims - 1),
        weights=q_weights,
        axis=dims - 2
    )

def clear_gamma_point(mat: numpy.ndarray):
    dims = len(mat.shape)
    indices = tuple([slice(None)] * (dims - 2) + [0, slice(0, 3)])
    mat[indices] = 0

class LogitudinalElasticModulusPhononContribution(ElasticModulus):

    # type: calculator = cij.core.calculator.Calculator

    def __init__(self, calculator, e: Tuple[float, float]):
        # e ought to be the e / dialation
        self.e = e
        self.calculator = calculator
        self.qha_calculator = self.calculator.qha_calculator

        self.nv = self.calculator.nv
        self.np = self.calculator.np
        self.nq = self.calculator.nq
        self.na = self.calculator.na

        if not numpy.isclose(e[0], e[1]): raise RuntimeError()

    @property
    def v_array(self):
        return self.calculator.v_array

    @property
    def t_array(self):
        return self.calculator.t_array
    
    @property
    def freq_array(self):
        return self.calculator.freq_array

    @property
    def q_weights(self):
        return numpy.array([ weight for coord, weight in self.calculator.qha_input.weights ])

    @LazyProperty
    def prefactors(self):
        return (
            1 / 5 / numpy.prod(self.e),
            (1 / 3 / self.e[0], 1 / 3 / self.e[1]),
            1 / 5 / numpy.prod(self.e)
        )

    @LazyProperty
    def mode_gamma(self):
        return (
            self.prefactors[0]    * self.calculator.mode_gamma[0],
            self.prefactors[1][0] * self.calculator.mode_gamma[1],
            self.prefactors[2]    * self.calculator.mode_gamma[2]
        )

    @LazyProperty
    def Q(self):
        return  h_div_k * (self.freq_array[nax,:,:,:] / self.t_array[:,nax,nax,nax])

    @LazyProperty
    def Q1(self):
        return self.Q / (numpy.exp(self.Q) - 1)

    @LazyProperty
    def Q2(self):
        return self.Q ** 2 * numpy.exp(self.Q) / (numpy.exp(self.Q) - 1) ** 2


    @LazyProperty
    def zero_point_contribution(self):
        h = units.Quantity(_h, units.J * units.m).to(units.rydberg * units.cm).magnitude
        return h / 2 / self.v_array \
            * self.average_over_modes(
                + self.mode_gamma[2] * self.freq_array
                - self.mode_gamma[0] * self.freq_array
                + self.mode_gamma[1] * self.freq_array
            ) * 3 * self.na

    @LazyProperty
    def thermal_contribution(self):

        k = units.Quantity(_k, units.eV / units.K).to(units.rydberg / units.K).magnitude

        ret = k * self.t_array[:, nax] / self.v_array[nax, :] \
            * self.average_over_modes(
                - self.Q2 * \
                    self.mode_gamma[2][nax,:,:,:] \
                + self.Q1 * (
                    + self.mode_gamma[2][nax,:,:,:]
                    - self.mode_gamma[0][nax,:,:,:]
                    + self.mode_gamma[1][nax,:,:,:]
                )
            ) * 3 * self.na

        ret[numpy.where(self.t_array == 0),:] = 0
        
        return ret

    @LazyProperty
    def value_isothermal(self):
        return self.zero_point_contribution + self.thermal_contribution

    @LazyProperty
    def isothermal_to_adiabatic(self):
        
        k = units.Quantity(_k, units.eV / units.K).to(units.rydberg / units.K).magnitude

        ret = self.t_array[:, nax] / self.v_array[nax, :] \
            / self.qha_calculator.volume_base.heat_capacity \
            * self.average_over_modes(self.Q2 * self.mode_gamma[1][0]) \
            * self.average_over_modes(self.Q2 * self.mode_gamma[1][1]) \
            * (3 * k * self.na) ** 2

        ret[numpy.where(self.t_array == 0), :] = 0

        return ret

    @property
    def value_adiabatic(self):
        return self.value_isothermal + self.isothermal_to_adiabatic

    def average_over_modes(self, amount):
        return average_over_modes(amount, self.q_weights)


class OffDiagnonalElasticModulusPhononContribution(LogitudinalElasticModulusPhononContribution):

    @LazyProperty
    def prefactors(self):
        return (
            1 / 15 / numpy.prod(self.e),
            (1 / 3 / self.e[0], 1 / 3 / self.e[1]),
            1 / 15 / numpy.prod(self.e)
        )

    @LazyProperty
    def mode_gamma(self):
        return (
            self.prefactors[0]    * self.calculator.mode_gamma[0],
            numpy.zeros(self.calculator.mode_gamma[1].shape),
            self.prefactors[2]    * self.calculator.mode_gamma[2]
        )

    @LazyProperty
    def zero_point_contribution(self):
        h = units.Quantity(_h, units.J * units.m).to(units.rydberg * units.cm).magnitude
        return h / 2 / self.v_array \
            * self.average_over_modes(
                + self.mode_gamma[2] * self.freq_array
                - self.mode_gamma[0] * self.freq_array
            ) * 3 * self.na

    @LazyProperty
    def thermal_contribution(self):

        k = units.Quantity(_k, units.eV / units.K).to(units.rydberg / units.K).magnitude

        ret = k * self.t_array[:, nax] / self.v_array[nax, :] \
            * self.average_over_modes(
                - self.Q2 * \
                    self.mode_gamma[2][nax,:,:,:] \
                + self.Q1 * (
                    + self.mode_gamma[2][nax,:,:,:]
                    - self.mode_gamma[0][nax,:,:,:]
                )
            ) * 3 * self.na

        ret[numpy.where(self.t_array == 0),:] = 0
        
        return ret

    @LazyProperty
    def value_isothermal(self):
        return self.zero_point_contribution + self.thermal_contribution + \
            self.qha_calculator.volume_base.pressures

class ShearElasticModulusPhononContribution(ElasticModulus):
    def __init__(self):
        self.value_adiabatic = None
        self.value_isothermal = None

