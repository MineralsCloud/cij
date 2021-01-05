import cij.data
import json
import jsonschema

__all__ = ["validate_config"]

def validate_config(config: dict) -> None:
    '''Validating the configuration object, raises exception when the
    configuration file is invalid.

    :param config: the configuration object to be validated
    :raise: ``jsonschema.exceptions.ValidationError``
    '''

    with open(cij.data.get_data_fname("schema/config.schema.json")) as fp:
        schema = json.load(fp)
    jsonschema.validate(instance=config, schema=schema)
