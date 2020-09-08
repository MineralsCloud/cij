import logging
import click
import cij

def run(config_fname: str):
    import cij.core.calculator
    calculator = cij.core.calculator.Calculator(config_fname)
    calculator.write_output()

@click.command()
@click.argument("settings_filename", type=click.Path(exists=True))
@click.version_option(version=cij.__version__, prog_name="Cij")
@click.option("--debug", default="INFO", type=click.Choice(logging._levelToName.values()), help="Logging level")
def main(settings_filename: str, debug: str):
    logging.basicConfig(level=debug)
    run(settings_filename)