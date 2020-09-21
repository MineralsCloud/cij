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

def average_over_modes(amount: numpy.ndarray, q_weights: numpy.ndarray) -> numpy.ndarray:
    '''Calculate sum of physical quantity over :math:`3N` phonon modes and
    :math:`N_q` :math:`q`-points (:math:`\\sum_{qm}`) except for accoustic modes
    at \\Gamma point (first :math:`q`-point and first three modes).

    :param amount: :math:`X_{qm}`
    :param q_weights: :math:`q`-point multiplicities :math:`w_q`

    :returns: Sum :math:`\\bar X = \sum_{qm} X_{qm} w_q`
    '''
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

class LongitudinalElasticModulusPhononContribution(ElasticModulus):
    '''Represents the phonon part of the longitudinal thermal
    elastic modulus :math:`c^\\text{ph}_{ii}(T, V)`

    .. math::
        c^\\text{ph}_{ii}(T, V)
            = c^\\text{zpm}_{ii}(V) + c^\\text{th}_{ii}(T, V)

    :param calculator:
    :param e: the "strain" corresponding to the subscript 
        (:math:`e_i/\\Delta` and `e_j/\\Delta`) for calculating strain-Gruneisen
        parameter :math:`\\gamma^{ii}_{qm}` from mode-Gruneisen parameter
        :math:`\\gamma_{qm}`.
    '''

    # type: calculator = cij.core.calculator.Calculator

    def __init__(self, calculator: 'cij.core.calculator.Calculator', e: Tuple[numpy.ndarray, numpy.ndarray]):
        # e[i] is of dimension [nax]
        # e ought to be the e / dialation
        self.e = e
        self.calculator = calculator
        self.qha_calculator = self.calculator.qha_calculator

        self.nv = self.calculator.nv
        self.np = self.calculator.np
        self.nq = self.calculator.nq
        self.na = self.calculator.na


    @property
    def v_array(self) -> numpy.ndarray:
        return self.calculator.v_array

    @property
    def t_array(self) -> numpy.ndarray:
        return self.calculator.t_array
    
    @property
    def freq_array(self) -> numpy.ndarray:
        return self.calculator.freq_array

    @property
    def q_weights(self) -> numpy.ndarray:
        '''The :math:`q`-points multiplicities or weights :math:`w_{q}`
        '''
        return numpy.array([ weight for coord, weight in self.calculator.qha_input.weights ])

    @LazyProperty
    def prefactors(self) -> Tuple[numpy.ndarray, Tuple[numpy.ndarray, numpy.ndarray], numpy.ndarray]:
        return (
            1 / 5 / numpy.prod(self.e, axis=0),
            (1 / 3 / self.e[0], 1 / 3 / self.e[1]),
            1 / 5 / numpy.prod(self.e, axis=0)
        )

    @LazyProperty
    def mode_gamma(self) -> tuple:
        '''Values related to the strain-Grüneisen parameter
        '''
        return (
            self.prefactors[0][:,nax,nax] * self.calculator.mode_gamma[0],
            (
                self.prefactors[1][0][:,nax,nax] * self.calculator.mode_gamma[1],
                self.prefactors[1][1][:,nax,nax] * self.calculator.mode_gamma[1],
            ),
            self.prefactors[2][:,nax,nax] * self.calculator.mode_gamma[2]
        )

    @LazyProperty
    def Q(self) -> numpy.ndarray:
        '''Value of expression :math:`Q_{qm}(T, V)`

        .. math::
            Q_{qm}(T, V) = \\frac{\hbar\omega_{qm}(V)}{k_\\text{B}T}
        '''
        return  h_div_k * (self.freq_array[nax,:,:,:] / self.t_array[:,nax,nax,nax])

    @LazyProperty
    def Q1(self) -> numpy.ndarray:
        '''Value of expression

        .. math::
            \\frac{Q_{qm}(T, V)}{\exp Q_{qm}(T, V) - 1}

        where :math:`Q_{qm}(T, V) = \\frac{\hbar\omega_{qm}(V)}{k_\\text{B}T}`
        '''
        return self.Q / (numpy.exp(self.Q) - 1)

    @LazyProperty
    def Q2(self) -> numpy.ndarray:
        '''Value of expression

        .. math::
            \\frac{Q_{qm}^2(T, V) \exp Q_{qm}(T, V) }{(\exp Q_{qm}(T, V) - 1) ^ 2}

        where :math:`Q_{qm}(T, V) = \\frac{\hbar\omega_{qm}(V)}{k_\\text{B}T}`
        '''
        return self.Q ** 2 * numpy.exp(self.Q) / (numpy.exp(self.Q) - 1) ** 2


    @LazyProperty
    def zero_point_contribution(self) -> numpy.ndarray:
        '''The zero-point motion contribution to the
        vibrational part of the longitudinal elastic modulus
        (:math:`c_{iiii}`, :math:`i = 1-3`).

        .. math::
            c^{\\text{zpm}}_{iiii}
                = \\frac{\hbar}{2V}\\sum_{qm}
                    \\left(\\frac{\\partial^2\\omega_{qm} (V)}{\\partial e_{ii} ^ 2}\\right)
                = \\frac{\hbar}{2V}\\sum_{qm}
                    \\left(\\gamma^{ii}_{qm}\\gamma^{ii}_{qm} - \\frac{\\partial \\gamma^{ii}_{qm}}{\\partial e_{ii}} + \gamma^{ii}_{qm}\\right) \\omega_{qm}
        '''
        h = units.Quantity(_h, units.J * units.m).to(units.rydberg * units.cm).magnitude
        return h / 2 / self.v_array \
            * self.average_over_modes(
                + self.mode_gamma[2] * self.freq_array
                - self.mode_gamma[0] * self.freq_array
                + self.mode_gamma[1][0] * self.freq_array
            ) * 3 * self.na

    @LazyProperty
    def thermal_contribution(self) -> numpy.ndarray:
        '''The thermal contribution :math:`c^\\text{th}` to the
        vibrational part of the longitudinal elastic modulus (:math:`c_{iiii}`).
        
        .. math::
        '''
        k = units.Quantity(_k, units.eV / units.K).to(units.rydberg / units.K).magnitude

        ret = k * self.t_array[:, nax] / self.v_array[nax, :] \
            * self.average_over_modes(
                - self.Q2 * \
                    self.mode_gamma[2][nax,:,:,:] \
                + self.Q1 * (
                    + self.mode_gamma[2][nax,:,:,:]
                    - self.mode_gamma[0][nax,:,:,:]
                    + self.mode_gamma[1][0][nax,:,:,:]
                )
            ) * 3 * self.na

        ret[numpy.where(self.t_array == 0),:] = 0

        return ret

    @LazyProperty
    def value_isothermal(self) -> numpy.ndarray:
        '''The phonon part of the elastic constant :math:`c^\\text{ph}_{ij}`
        as a function of temeprature and volume

        .. math::
            c^\\text{ph}_{ij}(T, V)
                = c^\\text{zpm}_{ij}(V) + c^\\text{th}_{ij}(T, V)

        '''
        return self.zero_point_contribution + self.thermal_contribution

    @LazyProperty
    def isothermal_to_adiabatic(self):
        '''The difference between isothermal :math:`c^\\text{T}_{ij}` and
        adiabatic :math:`c^\\text{S}_{ij}` thermal elastic modulus

        .. math::
            c^\\text{S}_{ij} - c^\\text{T}_{ij} = 
                \\frac{T}{V C_{V}}
                \\frac{\\partial S}{\\partial e_{i i}}
                \\frac{\\partial S}{\\partial e_{j j}}
        '''

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
        '''The adiabatic elastic modulus :math:`c^\\text{S}_{ij}(T, V)` as a
        function of temperature and volume, valid only for longitudinal and
        off-diagonal elastic modulus (:math:`i,j=1,2,3`)

        .. math::
            c^\\text{S}_{ij} =
                \\frac{T}{V C_{V}}
                \\frac{\\partial S}{\\partial e_{i i}}
                \\frac{\\partial S}{\\partial e_{j j}}
                + c^\\text{T}_{ij}
        '''
        return self.value_isothermal + self.isothermal_to_adiabatic

    def average_over_modes(self, amount):
        return average_over_modes(amount, self.q_weights)


class OffDiagonalElasticModulusPhononContribution(LongitudinalElasticModulusPhononContribution):

    @LazyProperty
    def prefactors(self) -> Tuple[numpy.ndarray, Tuple[numpy.ndarray, numpy.ndarray], numpy.ndarray]:
        return (
            1 / 15 / numpy.prod(self.e, axis=0),
            (1 / 3 / self.e[0], 1 / 3 / self.e[1]),
            1 / 15 / numpy.prod(self.e, axis=0)
        )

    @LazyProperty
    def mode_gamma(self):
        # print(self.prefactors[0][2], self.prefactors[1][0][2], self.prefactors[1][1][2], self.prefactors[2][2])
        return (
            self.prefactors[0][:,nax,nax] * self.calculator.mode_gamma[0],
            (
                self.prefactors[1][0][:,nax,nax] * self.calculator.mode_gamma[1],
                self.prefactors[1][1][:,nax,nax] * self.calculator.mode_gamma[1]
            ),
            self.prefactors[2][:,nax,nax] * self.calculator.mode_gamma[2]
        )

    @LazyProperty
    def zero_point_contribution(self):
        '''The zero-point motion contribution to the
        vibrational part of the longitudinal elastic modulus
        (:math:`c^\\text{zpm}_{iijj}`, :math:`i,j = 1-3`, :math:`i \\neq j`).

        .. math::
            c^{\\text{zpm}}_{iijj}
                = \\frac{\hbar}{2V}\\sum_{qm}
                    \\left(\\frac{\\partial^2\\omega_{qm} (V)}{\\partial e_{ii} \\partial e_{jj}}\\right)
                = \\frac{\hbar}{2V}\\sum_{qm}
                    \\left(\\gamma^{ii}_{qm}\\gamma^{jj}_{qm} - \\frac{\\partial \\gamma^{ii}_{qm}}{\\partial e_{jj}}\\right) \\omega_{qm}
        '''
        h = units.Quantity(_h, units.J * units.m).to(units.rydberg * units.cm).magnitude
        return h / 2 / self.v_array \
            * self.average_over_modes(
                + self.mode_gamma[2] * self.freq_array
                - self.mode_gamma[0] * self.freq_array
            ) * 3 * self.na

    @LazyProperty
    def thermal_contribution(self):
        '''The thermal contribution to the
        vibrational part of the longitudinal elastic modulus
        (:math:`c^\\text{th}_{iijj}`, :math:`i,j = 1-3`, :math:`i \\neq j`).

        .. math::
            c^{\\text{th}}_{iijj}
                = \\frac{k_{\\mathrm{B}} T}{V} \\sum_{q m} \\frac{\\partial^{2}\\left[\\ln \\left(1-e^{-Q_{q m}}\\right)\\right]}{\\partial e_{i i} \\partial e_{j j}}
                = 
        '''
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
        return self.zero_point_contribution + self.thermal_contribution + (
            + self.qha_calculator.volume_base.pressures
            - self.calculator.static_p_array[nax, :]
        )


