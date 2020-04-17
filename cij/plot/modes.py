import cij.core.calculator
import matplotlib
import numpy
from cij.util import _to_ang3

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
            ax.scatter(self.volumes, freqs)