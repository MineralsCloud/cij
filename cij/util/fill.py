from pathlib import Path
from collections import OrderedDict
import itertools
import re

import pandas
import sympy
import numpy
from sympy.parsing.sympy_parser import parse_expr

from cij.data import get_data_fname


def fill_cij(
    elast: pandas.DataFrame,
    system: str = None,
    ignore_residuals: bool = False,
    ignore_rank: bool = False,
    drop_atol: float = 1e-8 
) -> pandas.DataFrame:

    if system is None:
        return elast

    # use V column as index

#     elast = elast.set_index("V")

    # symbols = c11, c22, c33, ...

    symbols = OrderedDict([
        (f"c{i}{j}", sympy.symbols(f"c{i}{j}"))
        for (i, j)
        in itertools.product(range(1, 7), range(1, 7))
        if i <= j
    ])
    nsym = len(symbols)

    # known variables

    a = []
    b = []

    for sym, col in elast.iteritems():
        if not re.search(r"c(\d)(\d)", sym): continue
        _a = numpy.zeros(nsym)
        _a[list(symbols.keys()).index(sym)] = 1
        a.append(_a)
        b.append(col.to_numpy())

    a = numpy.array(a)
    b = numpy.array(b)


    # constraints

    # ... try to find constraints file

    if not Path(system).exists():
        constraints = Path("constraints") / system
        constraints = get_data_fname(str(constraints))

    # ... apply constraints

    eqns = []

    with open(constraints) as fp:
        for line in fp:
            parts = [parse_expr(part) for part in line.split("=")]
            for part in parts[1:]:
                eqns.append(parts[0] - part)
    
    _a, _b = sympy.linear_eq_to_matrix(eqns, *symbols.values())
    _a = numpy.array(_a).astype(numpy.float64)
    _b = numpy.array(_b).astype(numpy.float64)
    _b = numpy.broadcast_to(_b, (_b.shape[0], b.shape[1]))

    a = numpy.concatenate((a, _a), axis=0)
    b = numpy.concatenate((b, _b), axis=0)

    # fitting

    x, residuals, rank, s = numpy.linalg.lstsq(a, b)

    # warning

    if rank < nsym and not ignore_rank:
        raise Warning(f"Rank of constraints {rank} is smaller than input {nsym}!")
    
    if numpy.any(residuals > 0.1) and not ignore_residuals:
        raise Warning(f"Residuals seems to be too large! -> (" + ", ".join(
            f"{sym}: {res:.3f}" for sym, res in zip(symbols, residuals)
        ) + ")")
    
    # write back to table

    for index, col in zip(symbols, x):
        elast.loc[:, index] = col
    
    # drop empty columns

    for index, col in list(elast.iteritems()):
        if numpy.allclose(col.to_numpy(), 0, atol=drop_atol):
            elast = elast.drop(index, axis=1)

#     elast = elast.reset_index()

    return elast