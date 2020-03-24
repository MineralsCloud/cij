'''
Unit registery provided by `Pint <https://pint.readthedocs.io/>`_. Because the
entire program should share a same Pint unit register, ``units`` is the instance
that is created here and then shared upon.
'''

import pint

units = pint.UnitRegistry()

__all__ = ["units"]
