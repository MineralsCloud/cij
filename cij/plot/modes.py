import matplotlib
import matplotlib.pyplot as plt
import numpy

import cij.core.calculator
from cij.util.units import _from_ang3, _to_gpa, _to_ang3
from cij.plot.pvticker import PofVLocator, PofVFormatter
from cij.util.static_eos import get_static_p_of_v


class ModePlotter:
    '''The plotter for plotting interpolated phonon mode frequency and its derivatives (Gruneisen parameters, etc.) for checking
    the status of mode-frequency interpolation
    '''
    def __init__(self, calculator: cij.core.calculator.Calculator):
        self.calculator = calculator

    @property
    def qha_input(self):
        return self.calculator.qha_input

    @property
    def volumes(self):
        return _to_ang3(numpy.array([
            volume.volume
            for volume in self.qha_input.volumes
        ]))
    
    @property
    def v_array(self):
        return _to_ang3(self.calculator.v_array)

    def plot_modes(self, ax: matplotlib.axes.Axes, n: int = 0, iq: int = 0):
        '''Plotting the phonon modes and its derivatives

        :param ax:
        :param n: the order of derivatives to be plotted:
            :math:`n = 0` for :math:`\\omega_{qm}(V)`,
            :math:`n = 1` for :math:`\\gamma_{qm}(V)`,
            :math:`n = 2` for :math:`V\\frac{\partial\gamma_{qm}(V)}{\partial V}`
        :param iq: the index of :math:`q` point to be plotted
        '''

        if   n == 0:
            w_arrays = self.calculator.freq_array[:, iq, :]
        elif n == 1:
            w_arrays = self.calculator.mode_gamma[0][:, iq, :]
        elif n == 2:
            w_arrays = self.calculator.mode_gamma[1][:, iq, :]

        for k in range(self.calculator.np):
            if iq == 0 and k < 3: continue
            w_array = w_arrays[:, k]
            ax.plot(self.v_array, w_array)

        if n != 0: return

        for k in range(self.calculator.np):
            if iq == 0 and k < 3: continue
            freqs = numpy.array([
                volume.q_points[iq].modes[k]
                for volume in self.qha_input.volumes
            ])
            ax.scatter(self.volumes, freqs, s=10)

    def add_pressure_ticks(self, ax: matplotlib.axes.Axes, interval: float = 1.0, label: str = None):

        v_static, f_static = zip(*[(v_data.volume, v_data.energy) for v_data in self.calculator.qha_input.volumes])
        static_p_of_v = lambda v: _to_gpa(get_static_p_of_v(v_static, f_static)(_from_ang3(v)))

        ndec = max(-int(numpy.log10(interval)), 0)

        __ax = plt.gca()

        _ax = ax.twiny()
        _ax.set_xlim(*ax.get_xlim())
        _ax.xaxis.set_major_locator(PofVLocator(static_p_of_v, p_interval=interval))
        _ax.xaxis.set_minor_locator(PofVLocator(static_p_of_v, p_interval=interval / 5))
        _ax.xaxis.set_major_formatter(PofVFormatter(static_p_of_v, ndec=ndec))

        if label:
            _ax.set_xlabel(label)

        plt.sca(__ax)
        