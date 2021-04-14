import logging
import click
from pathlib import Path

with open(Path(__file__).parent / "../version.py") as fp: exec(fp.read())

def run(config_fname: str):
    import cij.core.calculator
    calculator = cij.core.calculator.Calculator(config_fname)
    calculator.write_output()
    
@click.command("run", help="Perform SAM-Cij calculation details specified in SETTINGS.")
@click.argument("settings", type=click.Path(exists=True))
@click.version_option(version=__version__, prog_name="Cij")     # pylint: disable=undefined-variable
@click.option("--debug", default="INFO", type=click.Choice(logging._levelToName.values()), help="Verbosity level of debug log emitted to the standard output.")
def main(settings: str, debug: str):

    logger = logging.getLogger("cij")
    logger.setLevel(debug)
    handler = logging.StreamHandler()
    logger.addHandler(handler)

    run(settings)