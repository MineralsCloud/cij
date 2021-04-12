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

@pytest.mark.parametrize("config", [{}])
def test_apply_default_config_changes(config):
    assert apply_default_config(config) != config

@pytest.mark.parametrize( "a, b, ex", [
    ({"a": 1}, {"b": 2}, {"a": 1, "b": 2}),
    ({"a": 1}, {"a": 2}, {"a": 1}),
    ({"a": 1}, {"a": 1}, {"a": 1}),
    ({"a": {"b": 1}}, {"a": {"c": 2}}, {"a": {"b": 1, "c": 2}}),
    ({"b": 3}, {"a": {"c": 2}}, {"a": {"c": 2}, "b": 3})
])
def test_update_config(a, b, ex):
    assert update_config(a, b) == ex

from jsonschema.exceptions import ValidationError

@pytest.mark.parametrize("config", [
    {"qha": {}},
    {"elast": {}},
    {"elast": {"settings": {"symmetry": {"system": 2}}}, "qha": {}},
    {"elast": {"settings": {"symmetry": {"system": "x"}}}, "qha": {}},
    {"elast": {"settings": {"symmetry": {"system": "orthrohombic"}}}, "qha": {}},
    {"elast": {"settings": {"symmetry": {"system": "orthrohombic"}}}}
])
def test_validate_config2(config):
    with pytest.raises(ValidationError):
        validate_config(config)


@pytest.mark.parametrize("config", [
    { "elast": { "settings": { "symmetry": { "system": "orthorhombic" } } }, "qha": {} }
])
def test_validate_config3(config):
    validate_config(config)
