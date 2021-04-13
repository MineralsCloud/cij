import pytest
import numpy

from cij.util.fill import fill_cij


input02s = [
    ("examples/diopside/input02", "monoclinic"),
    ("examples/bridgmanite/elast.dat", "orthorhombic"),
    ("examples/akimotoite/input02", "trigonal7"),
]



@pytest.fixture
def input02(request):
    import pandas, io
    with open(request.param) as fp:
        for line in fp:
            if line.strip()[0].isdigit():
                break
        N = int(line.strip().split()[1])
        sio = io.StringIO("".join(fp.readlines()[:N]))
        df = pandas.read_table(sio, index_col=None, header=0, sep=r"\s+")
    return df

@pytest.fixture
def system(request):
    return request.param

input02s = [
    ("examples/diopside/input02", "monoclinic"),
    ("examples/bridgmanite/elast.dat", "orthorhombic"),
    ("examples/akimotoite/input02", "trigonal7"),
]


@pytest.mark.parametrize("input02, system", input02s, indirect=True)
def test_cij_fill_correctly_recognize_uppercase(cli_runner, input02, system):

    src = input02
    src.columns = [k.upper() for k in src.columns]
    dest = fill_cij(src, system=system)

    assert len(set(dest.columns)) == len(dest.columns)
    assert len(set(c.lower() for c in dest.columns)) == len(dest.columns)
    assert dest.notnull().values.any()


@pytest.mark.parametrize("input02, system", input02s, indirect=True)
def test_cij_fill_correctly_raise_error_with_wrong_system(cli_runner, input02, system):
    with pytest.raises(Warning):
        fill_cij(input02, system="cubic")


@pytest.mark.parametrize("input02, system", input02s, indirect=True)
def test_cij_fill_correctly_raise_error_with_missing_column(cli_runner, input02, system):
    with pytest.raises(Warning):
        fill_cij(input02.iloc[:, :3], system=system)


@pytest.mark.parametrize("input02, system", [("examples/akimotoite/input02", "trigonal7")], indirect=True)
@pytest.mark.parametrize("key", ["c11", "c12"])
def test_cij_fill_back_c22(cli_runner, input02, system, key):
    assert key in input02.columns
    broken = input02.drop(key, axis=1)
    assert key not in broken
    output = fill_cij(broken, system=system)
    assert key in output.columns
    assert numpy.allclose(output[key].values, input02[key].values)