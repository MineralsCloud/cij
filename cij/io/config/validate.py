import cij.data
import json
import jsonschema

__all__ = "validate_config"

def validate_config(config):
    with open(cij.data.__path__ / "schema" / "config.schema.json") as fp:
        schema = json.load(fp)
    jsonschema.validate(instance=config, schema=schema)
