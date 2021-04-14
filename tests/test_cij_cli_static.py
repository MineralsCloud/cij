import pytest


inputs = [
    ("examples/diopside/input01", "examples/diopside/input02", "monoclinic"),
    ("examples/bridgmanite/input01", "examples/bridgmanite/elast.dat", "orthorhombic"),
    ("examples/akimotoite/input01", "examples/akimotoite/input02", "trigonal7"),
]

@pytest.mark.parametrize("input01, input02, system", inputs)
def test_cij_static_with_example(cli_runner, input01, input02, system):
    from cij.cli.static import main

    result = cli_runner.invoke(main, [input01])
    assert result.exit_code == 0

    result = cli_runner.invoke(main, [input01, input02, "-s", system])
    assert result.exit_code == 0
