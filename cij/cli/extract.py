import click
import pandas
from typing import List
import numpy
import glob

def load_data(var):
    df = pandas.read_table(glob.glob(f"{var}_tp_*")[0], sep="\s+", index_col=0)
    df.columns = [float(colname) for colname in df.columns]
    df.index = [float(idx) for idx in df.index]
    return df

@click.command(help="Extract the value and create table for multiple variables at specific P or T.")
@click.option("-T", "--temperature", type=click.FLOAT, help="Specify temperature in K.")
@click.option("-P", "--pressure", type=click.FLOAT, help="Specify tressure in GPa.")
@click.option("-v", "--variables", required=True, help="Variables to output.")
@click.option("-h", "--hide-header", default=False, is_flag=True, help="Hide header or not.")
def main(variables: List[str], hide_header: bool, temperature: float = None, pressure: float = None):

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