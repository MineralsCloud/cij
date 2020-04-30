'''The plotter that helps you to plot physical properties inside ``qha-cij``'s
calulator, without the worries of dealing with boundaries and units.
'''

from cij.core.calculator import Calculator
from cij.util import C_, c_, units
import numpy
import matplotlib.axes
from cij.util import _from_gpa, _to_gpa


class Plotter:
    '''The plotter for QHA's calculator
    '''
    def __init__(self, calculator: Calculator):
        self.calculator = calculator

    def index_t(self, t: float):
        t_array = self.calculator.pressure_base.t_array
        t_index = numpy.argmin(numpy.abs(t_array - t))
        return t_index

    def index_p(self, p: float):
        p_array = self.calculator.pressure_base.p_array
        p_index = numpy.argmin(numpy.abs(p_array - p))
        return p_index

    @property
    def p_indices(self):
        return self.calculator.qha_calculator.pressure_base.volumes < numpy.max([
            volume.volume for volume in self.calculator.qha_input.volumes
        ])

    def plot_cij_p(self, ax: matplotlib.axes.Axes, cij_key: int, t: float, *argv, **kwargs):
        
        key = c_(cij_key)
        t_index = self.index_t(t)

        # p_array = self.calculator.pressure_base.p_array[p_indices[t_index]]
        # c_array = self.calculator.modulus_adiabatic[key][t_index, p_indices[t_index]]

        p_array = self.calculator.pressure_base.p_array
        c_array = self.calculator.pressure_base.modulus_adiabatic[key][t_index, :]

        line, = ax.plot(
            _to_gpa(p_array),
            _to_gpa(c_array),
            *argv,
            **kwargs
        )

        return line

    def plot_cij_p_with(self, handler: callable, cij_key: int, t: float):
        
        key = c_(cij_key)
        t_index = self.index_t(t)

        # p_array = self.calculator.pressure_base.p_array[p_indices[t_index]]
        # c_array = self.calculator.modulus_adiabatic[key][t_index, p_indices[t_index]]

        p_array = self.calculator.pressure_base.p_array
        c_array = self.calculator.pressure_base.modulus_adiabatic[key][t_index, :]

        return handler(
            _to_gpa(p_array),
            _to_gpa(c_array),
        )


    def plot_cij_t_with(self, handler: callable, cij_key: int, p: float):
        
        key = c_(cij_key)
        p_index = self.index_p(_from_gpa(p))

        # p_array = self.calculator.pressure_base.p_array[p_indices[t_index]]
        # c_array = self.calculator.modulus_adiabatic[key][t_index, p_indices[t_index]]

        t_array = self.calculator.pressure_base.t_array
        c_array = self.calculator.pressure_base.modulus_adiabatic[key][:, p_index]

        return handler(
            t_array,
            _to_gpa(c_array),
        )

