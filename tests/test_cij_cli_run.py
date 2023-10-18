import pytest
from pathlib import Path
from os import getcwd, chdir

settings = [
    "examples/diopside/settings.yaml",
    "examples/bridgmanite/settings.yaml",
    "examples/akimotoite/settings.yaml",
]

@pytest.mark.parametrize("settings", settings)
def test_cij_run_with_example(cli_runner, settings):
    import cij.cli.main
    settings = Path(settings)
    working_dir = str(settings.parent)
    settings_fname = settings.name
    cwd = getcwd()
    chdir(working_dir)
    result = cli_runner.invoke(cij.cli.main.main, [settings_fname])
    chdir(cwd)

    assert result.exit_code == 0
