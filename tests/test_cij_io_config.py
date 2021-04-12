import pytest
from glob import glob

from cij.io.config import read_config, validate_config

@pytest.yield_fixture
def all_config():
    return [read_config(fname) for fname in glob("examples/*/settings.yaml")]


def test_validaiton(all_config):
    for c in all_config:
        validate_config(c)