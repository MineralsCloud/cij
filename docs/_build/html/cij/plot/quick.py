
import numpy
from matplotlib import pyplot as plt
from pathlib import Path
from .color_map import color_x
from typing import NamedTuple
from cij.util import C_, c_, units

def _guess_unit(fname: str):
    import re

    from cij.io.output.results_writer import DEFAULT_WRITER_RULES
    for rule in DEFAULT_WRITER_RULES:
        fname_pattern = rule["fname_pattern"]
        regex_pattern = re.compile(
            fname_pattern.format(
                ij=r"(?P<ij>\d+)",
                base=r"(?P<base>[vpt]{2,2})"
            )
        )
        res = regex_pattern.search(fname)
        if res: break
    
    if res is None:
        raise RuntimeError(f"File name '{fname}' has no match in rules!")

    key = None
    if 'ij' in dict(res.groups()).keys():
        key = c_(res.group('ij'))
    
    return rule, key, res.group('base')

def _guess_base(base: str):

    def _guess(c):
        return {
            "v": ("$V$", "Å^3"),
            "p": ("$P$", "GPa"),
            "t": ("$T$", "K")
        }[c]
    
    return [_guess(c) for c in base]

        
class PlotUnits(NamedTuple):

    x_label: str
    y_label: str
    z_label: str
    x_unit: str
    y_unit: str
    z_unit: str

    @classmethod
    def guess(cls, fname) -> 'PlotUnits':

        rule, key, base = _guess_unit(fname)
        _base = _guess_base(base)

        return cls(
            x_label=_base[1][0],
            y_label=_base[0][0],
            z_label=rule["description"],
            x_unit = _base[1][1],
            y_unit=_base[0][1],
            z_unit=rule["unit"]
        )

def read_table(fname: str, transpose: bool = False):
    with open(fname) as fp:
        first = next(fp).strip().split()
        x_array = [float(x) for x in first[1:]]
        arr = numpy.loadtxt(fp, dtype=float)
        y_array = arr[:, 0]
        z_data = arr[:, 1:]
    if not transpose:
        return x_array, y_array, z_data
    else:
        return y_array, x_array, z_data.T

def plot_table(fname: str, transpose: bool = False):

    plt.figure()

    x_array, y_array, z_data = read_table(fname, transpose)

    _units = PlotUnits.guess(fname)

    ymin = numpy.min(y_array)
    ymax = numpy.max(y_array)
    cmap = color_x(ymin, ymax)

    for y, z_array in zip(y_array, z_data):
        label = "{y:.0f} {unit}".format(y=y, unit=_units.y_unit)
        plt.plot(x_array, z_array, label=label, c=cmap(y))

    plt.xlabel("{label} ({unit})".format(label=_units.x_label, unit=_units.x_unit))
    plt.ylabel("{label} ({unit})".format(label=_units.z_label, unit=_units.z_unit))

    plt.legend(fontsize=8)
    plt.xlim(min(x_array), max(x_array))

    stem = Path(fname).stem

    plt.savefig(f"{stem}.png", dpi=300)
    plt.close()
