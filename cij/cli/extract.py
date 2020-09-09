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

@click.command()
@click.option("-T", "--temperature", type=click.FLOAT, required=True, help="Temperature in K")
@click.option("-v", "--variables", required=True, help="Variables for output")
@click.option("-h", "--hide-header", default=False, is_flag=True, help="Hide header or not")
def main(temperature: float, variables: List[str], hide_header: bool):

    data = {}

    variables = variables.split(",")
    for var in variables:
        df = load_data(var)
        t_index = numpy.argmin(numpy.abs(df.index.to_numpy() - temperature))
        data[var] = df.iloc[t_index]

    pressures = df.columns
    
    table = pandas.DataFrame(columns=variables, index=pressures)
    for var in variables:
        table[var] = data[var]
    
    print(table.to_string(header=(not hide_header)))

if __name__ == "__main__":
    main()