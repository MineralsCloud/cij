import pytest
from glob import glob

from cij.io.config import read_config, validate_config, update_config, apply_default_config

files = [*glob("examples/*/settings.yaml"), "cij/data/default/settings.yaml"]

@pytest.fixture
def config(request):
    return read_config(request.param)

@pytest.mark.parametrize("config", files, indirect=True)
def test_validate_config(config):
    validate_config(config)

@pytest.mark.parametrize("config", files, indirect=True)
def test_update_config(config):
    assert update_config(config, default_dict={}) == config
    assert update_config({}, default_dict=config) == config


@pytest.mark.parametrize("config", files, indirect=True)
def test_apply_default_config(config):
    assert apply_default_config(config) == config