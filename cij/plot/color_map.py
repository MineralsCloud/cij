'''Determine color based on a upper and lower bound
'''

import matplotlib
from typing import Callable, Optional, Union

def color_x(
    vmin: float, vmax: float,
    cmap: Optional[Union[matplotlib.colors.Colormap, str]] = None
) -> Callable[[float], tuple]:
    '''Generate a map from :math:`x` to color
    
    :param vmin: the lower bound of the color map
    :param vmax: the upper bound of the color map

    :returns: the function :math:`C(x)` that accepts :math:`x` and returns color :math:`C`
    '''

    if cmap == None:
        cmap ="Blues"

    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    cmap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)

    def get_color(v: float) -> tuple:
        return cmap.to_rgba(v)

    return get_color
