'''Determine color based on a upper and lower bound
'''

import matplotlib
from typing import Callable

def color_x(
    vmin: float, vmax: float,
    cmap: matplotlib.colors.Colormap = None
) -> Callable[[float], tuple]:
    '''Generate a map from :math:`x` to color
    
    :param vmin: the lower bound of the color map
    :param vmax: the upper bound of the color map

    :returns: the function :math:`C(x)` that accepts :math:`x` and returns color :math:`C`
    '''
    
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    cmap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
    def get_color(v: float) -> tuple:
        return cmap.to_rgba(v)
    return get_color
