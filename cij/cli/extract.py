import click
from typing import List


def load_data(var):

    import pandas
    from glob import glob

    df = pandas.read_table(glob(f"{var}_tp_*")[0], sep="\s+", index_col=0)
    df.columns = [float(colname) for colname in df.columns]
    df.index = [float(idx) for idx in df.index]
    return df

@click.command("extract", help="Extract data from cij calculation results to a table at specific P or T (e.g., table with c11s, c12s, K, G vs. P at 300 K).")
@click.option("-v", "--variables", required=True, help="Variables to output, (e.g., 'c11s,c12s,bm,G'), values should be seperated with comma.")
@click.option("-T", "--temperature", type=click.FLOAT, help="Specify temperature to extract data at, value in unit K, result will be tabulated vs. P.")
@click.option("-P", "--pressure", type=click.FLOAT, help="Specify pressure to extract data at, value in unit GPa, result will be tabulated vs. T.")
@click.option("-h", "--hide-header", default=False, is_flag=True, help="Hide table header or not.", show_default=True)
def main(variables: List[str], hide_header: bool, temperature: float = None, pressure: float = None):

    import glob
    import pandas
    import numpy

    data = {}

    variables = variables.split(",")

    for var in variables:
        df = load_data(var)

        if temperature != None:
            y = temperature
        elif pressure != None:
            y = pressure
            df = df.T

        x_array = df.columns

        y_index = numpy.argmin(numpy.abs(df.index.to_numpy() - y))
        data[var] = df.iloc[y_index]
    
    table = pandas.DataFrame(columns=variables, index=x_array)
    for var in variables:
        table[var] = data[var]
    
    print(table.to_string(header=(not hide_header)))

if __name__ == "__main__":
    main()