import math
import numpy
from scipy.interpolate import InterpolatedUnivariateSpline
from matplotlib.ticker import Locator, Formatter
from typing import Optional
from cij.util.units import _to_gpa, _from_ang3

class PofVLocator(Locator):

    def __init__(self, p_of_v: callable, p_interval: float = 2.0):

        super().__init__()

        self.p_of_v = p_of_v
        self.p_interval = p_interval

    def _get_v_of_p(self, vmin, vmax, N = 100):

        v_array = numpy.linspace(vmin, vmax, N + 1)
        p_array = self.p_of_v(v_array)

        return InterpolatedUnivariateSpline(p_array[::-1], v_array[::-1])
    
    def tick_values(self, vmin, vmax):

        vmin, vmax = sorted([vmin, vmax])

        pmin = self.p_of_v(vmax)
        pmax = self.p_of_v(vmin)

        imin = math.ceil(pmin / self.p_interval)
        imax = math.floor(pmax / self.p_interval)

        ni = imax - imin

        if ni < 0 or ni > 200: raise RuntimeError("Wrong PtoV conversion!")

        v_of_p = self._get_v_of_p(vmin, vmax)

        pticks = numpy.arange(imin, imax + 1) * self.p_interval
        vticks = v_of_p(pticks)

        return vticks

    def __call__(self):
        vmin, vmax = self.autoscale()
        return self.tick_values(vmin, vmax)

class PofVFormatter(Formatter):

    def __init__(self, p_of_v: callable, ndec: int = 1):

        super().__init__()

        self.p_of_v = p_of_v
        self.ndec = ndec

    def format_data(self, value):

        return ("%." + str(self.ndec) + "f") % numpy.round(self.p_of_v(value), self.ndec)

    def __call__(self, x, pos=None):

        return self.format_data(x)