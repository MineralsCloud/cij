import click
from cij import __version__


@click.group()
@click.version_option(version=__version__, prog_name="Cij")
def main():
    pass

from cij.cli.main import main as _run
main.add_command(_run, "run")

from cij.cli.static import main as _static
main.add_command(_static, "run-static")

from cij.cli.fill import main as _fill
main.add_command(_fill, "fill")

from cij.cli.extract import main as _extract
main.add_command(_extract, "extract")

if __name__ == "__main__":
    main()