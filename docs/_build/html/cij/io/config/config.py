from pathlib import Path
from .validate import validate_config
from typing import Union

def read_config(fname: Union[str, Path], validate: bool = True) -> dict:
    '''Reading JSON and YAML configuration.

    :param fname: name or path of the configuration file; type of parser used
        will be determined by the extension of the configuration file
        (``.json``, ``.yml`` or ``.yaml``).
    :param validate: whether the input needs to be validated
    :return: the configuration object
    '''

    suffix = Path(fname).suffix
    with open(fname) as fp:
        if suffix in {".yml", ".yaml"}:
            import yaml
            config = yaml.load(fp, Loader=yaml.FullLoader)
        elif suffix in {".json"}:
            import json
            config = json.load(fp)
        else:
            raise RuntimeError(f"File type {suffix} is not supported")

    if validate:
        validate_config(config)

    return config
