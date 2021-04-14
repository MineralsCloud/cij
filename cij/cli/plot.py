import sys
from glob import glob
import click


@click.command("plot", help="Plot SAM-Cij calculation results for files matching PATTERNS (e.g., *.txt, bm_*.txt, c*_*.txt or full path bm_tp_gpa.txt are all valid; multiple values can be seperated by spaces).")
@click.argument("patterns", nargs=-1)
def main(patterns):

    from cij.plot.quick import plot_table

    for pattern in patterns:
        for fname in glob(pattern):
            plot_table(fname)