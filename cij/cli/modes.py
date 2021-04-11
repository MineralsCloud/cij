import click

@click.command("modes", help="Plot interpolated mode frequency vs volume.")
@click.argument("fname", type=click.Path(exists=True))
@click.option("-q", "--iq", type=click.INT, default=0, help="i-th q-point")
@click.option("-n", type=click.IntRange(0, 3), default=0, help="n-th derivative") 
@click.option("-o", "--output", help="save figure to", default="modes.png") 
@click.option("--y-max", help="y-max limit", default=None, type=click.FLOAT) 
@click.option("-p", "--interval", help="p-tick interval", default=0, type=click.FLOAT) 
def main(fname: str, iq: int, n: int, interval: float, output: str, y_max: float):

    import cij.core.calculator

    calculator = cij.core.calculator.Calculator(fname)

    from cij.plot import ModePlotter
    from matplotlib import pyplot as plt

    plotter = ModePlotter(calculator)

    plt.figure(figsize=(5,10))

    plotter.plot_modes(plt, n, iq)

    plt.xlabel("Volume (Å$^3$)")
    plt.ylabel("Frequency (cm$^{-1}$)")

    if interval:
        plotter.add_pressure_ticks(plt.gca(), label = "Pressure (GPa)", interval=interval)

    plt.xlim(max(plotter.volumes), min(plotter.volumes))
    plt.ylim(bottom=0)
    if y_max:
        plt.ylim(top=y_max)

    plt.savefig(output)

if __name__ == "__main__":
    main()