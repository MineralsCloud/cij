import pytest

input02s = [
    ("examples/diopside/input02", "monoclinic"),
    ("examples/bridgmanite/elast.dat", "orthorhombic"),
    ("examples/akimotoite/input02", "trigonal7"),
]

@pytest.mark.parametrize("input02, system", input02s)
def test_cij_fill_with_example(cli_runner, input02, system):
    import cij.cli.fill
    result = cli_runner.invoke(cij.cli.fill.main, ["-s", system, input02])
    assert result.exit_code == 0


@pytest.mark.parametrize("input02, system", input02s)
def test_cij_fill_with_example_no_symmetry(cli_runner, input02, system):
    import cij.cli.fill
    result = cli_runner.invoke(cij.cli.fill.main, [input02])
    assert result.exit_code == 0


