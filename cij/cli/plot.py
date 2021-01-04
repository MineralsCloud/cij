import sys
from glob import glob
import click


@click.command(help="Plot SAM-Cij calculation results.")
@click.argument("patterns", nargs=-1)
def main(patterns):

    from cij.plot.quick import plot_table

    for pattern in patterns:
        for fname in glob(pattern):
            plot_table(fname)