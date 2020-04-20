'''
Unit registery provided by `Pint <https://pint.readthedocs.io/>`_. Because the
entire program should share a same Pint unit register, ``units`` is the instance
that is created here and then shared upon.
'''

from typing import TypeVar
import pint

_T = TypeVar('_T') 

units = pint.UnitRegistry()

__all__ = [
    "units",
    "convert_unit",
    "_from_gpa", "_to_gpa",
    "_from_ang3", "_to_ang3"
]

def convert_unit(unit_from: units.Unit, unit_to: units.Unit, value: _T = None) -> _T:
    convert = lambda x: units.Quantity(x, unit_from).to(unit_to).magnitude
    if value is not None:
        return convert(value)
    else:
        return convert

def _to_gpa(value: _T) -> _T:
    return convert_unit(
        units.rydberg / units.bohr ** 3,
        units.GPa,
        value
    )

def _from_gpa(value: _T) -> _T:
    return convert_unit(
        units.GPa,
        units.rydberg / units.bohr ** 3,
        value
    )

def _to_ang3(value: _T) -> _T:
    return convert_unit(
        units.bohr ** 3,
        units.angstrom ** 3,
        value
    )

def _from_ang3(value: _T) -> _T:
    return convert_unit(
        units.angstrom ** 3,
        units.bohr ** 3,
        value
    )