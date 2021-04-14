import click
from typing import List


def load_data(var):

    import pandas
    from glob import glob

    df = pandas.read_table(glob(f"{var}_tp_*")[0], sep=r"\s+", index_col=0)
    df.columns = [float(colname) for colname in df.columns]
    df.index = [float(idx) for idx in df.index]
    return df

def fit_data(df):

    from scipy.interpolate import RectBivariateSpline

    x = df.index.to_numpy()
    y = df.columns.to_numpy()
    z = df.to_numpy()
    return RectBivariateSpline(x, y, z)

@click.command("geotherm", help="Extract data from cij calculation results to table along a geotherm P, T and D (depth) map given in PATH.")
@click.option("-i", "-g", "--geotherm", required=True, help="The file name of geotherm P, D, T map.", type=click.Path(exists=True))
@click.option("--t-col", help="The name of geotherm pressure column.", default="P", show_default=True)
@click.option("--p-col", help="The name of geotherm temperature column", default="T", show_default=True)
@click.option("-v", "--variables", required=True, help="Variables to output, (e.g., 'c11s,c12s,bm,G'), values should be seperated with comma.")
@click.option("-h", "--hide-header", default=False, is_flag=True, help="Hide header or not.", show_default=True)
def main(variables: List[str], hide_header: bool, t_col: str = None, p_col: str = None, geotherm: str = None):

    import glob
    import pandas
    import numpy


    variables = variables.split(",")
    table = pandas.read_table(geotherm, sep=r"\s+", index_col=None, header=0)

    for var in variables:
        df = load_data(var)
        table[var] = fit_data(df)(table[p_col], table[t_col], grid=False)

    print(table.to_string(header=(not hide_header), index=False))

if __name__ == "__main__":
    main()