import click

@click.command("modes", help="Plot interpolated mode frequency vs volume for calculation details specified in SETTINGS.")
@click.argument("settings", type=click.Path(exists=True))
@click.option("-q", "--iq", type=click.INT, default=0, help="Plot for i-th q-point", show_default=True)
@click.option("-n", type=click.IntRange(0, 3), default=0, help="Plot for n-th derivative", show_default=True) 
@click.option("-o", "--output", help="File name for the output figure.", default="modes.png", type=click.Path(), show_default=True) 
@click.option("--y-max", help="y-max limit", default=None, type=click.FLOAT) 
@click.option("-p", "--interval", help="p-tick interval", default=0, type=click.FLOAT) 
def main(settings: str, iq: int, n: int, interval: float, output: str, y_max: float):

    import cij.core.calculator

    calculator = cij.core.calculator.Calculator(settings)

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