import click

@click.command("fill", help="Fill Cij tensor components table INPUT02 (elast.dat) based on the symmetry constraints of the crystal system.")
@click.argument("input02", type=click.Path(exists=True))
@click.option("-s", "--system", help="Name of the crystal system whose symmetry is applied to fill the missing elastic tensor components. Should be Should be one of: triclinic, monoclinic, hexagonal, trigonal6, trigonal7, orthorhombic, tetragonal6, tetragonal7, cubic.")
@click.option("--ignore-residuals", is_flag=True, help="Keep filling the missing elastic tensor components based on the crystal system even if the disagreements exceeds tolerance.")
@click.option("--ignore-rank", is_flag=True, help="Keep filling the missing elastic tensor components based on the crystal system even if not all necessary components are given")
@click.option("--drop-atol", type=click.FLOAT, default=1e-8, help="Drop the elastic tensor components if the magnitude of the elastic tenor components if this components in all volume points are less than this value.")
def main(**kwargs):

    import pandas
    import sys
    from io import StringIO
    from pathlib import Path
    
    from cij.util.fill import fill_cij

    fname = kwargs.pop("input02")

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

        # Write elasticity table

        elast = pandas.read_table(sio, header=0, index_col=None, sep="\s+")
        elast = fill_cij(elast, **kwargs)
        sys.stdout.write(elast.to_string(index=False) + "\n")

        # Write the rest of the file

        sys.stdout.write(fp.read())

if __name__ == "__main__":
    main()
