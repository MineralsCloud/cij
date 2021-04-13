import click

@click.command("fill", help="Fill non-zero Cij terms based on symmetry.")
@click.argument("input02", type=click.Path(exists=True))
@click.option("-s", "--system")
@click.option("--ignore-residuals", is_flag=True)
@click.option("--ignore-rank", is_flag=True)
@click.option("--drop-atol", type=click.FLOAT, default=1e-8)
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
