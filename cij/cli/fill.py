from io import StringIO
import sys

import click
import pandas

from pathlib import Path

from cij.util.fill import fill_cij
from cij.data import get_data_fname


@click.command(help="Fill non-zero Cij terms based on symmetry.")
@click.argument("fname", type=click.Path(exists=True))
@click.option("-c", "--constraints")
@click.option("--ignore-residuals", is_flag=True)
@click.option("--ignore-rank", is_flag=True)
@click.option("--drop-atol", type=click.FLOAT)
def main(**kwargs):

    fname = kwargs.pop("fname")

    with open(fname) as fp:

        # Read header

        sys.stdout.write(fp.readline())

        line = fp.readline()
        N = int(line.strip().split()[1])
        sys.stdout.write(line)

        # Read elasticity table

        sio = StringIO()
        for i in range(N + 1):
            line = fp.readline()
            sio.write(line)
        sio.seek(0)

        # Convert the constraint name to its full path


        if kwargs["constraints"] is not None:
            constraints_fname = Path("constraints") / kwargs["constraints"]
            kwargs["constraints"] = get_data_fname(str(constraints_fname))

        # Write elasticity table

        elast = pandas.read_table(sio, header=0, index_col=None, sep="\s+")
        elast = fill_cij(elast, **kwargs)
        sys.stdout.write(elast.to_string(index=False) + "\n")

        # Write the rest of the file

        sys.stdout.write(fp.read())

if __name__ == "__main__":
    main()
