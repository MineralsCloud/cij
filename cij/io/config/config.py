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


def update_config(input_dict: dict, default_dict: dict) -> dict:
    output_dict = {}
    for k in set([*input_dict.keys(), *default_dict.keys()]):
        if k not in input_dict.keys():
            output_dict[k] = default_dict[k]
        elif k not in default_dict.keys():
            output_dict[k] = input_dict[k]
        elif isinstance(input_dict[k], dict):
            output_dict[k] = update_config(input_dict[k], default_dict[k])
        else:
            output_dict[k] = input_dict[k]
    return output_dict


def apply_default_config(input_dict: dict) -> dict:
    import yaml
    import cij.data
    with open(cij.data.get_data_fname("default/settings.yaml")) as fp:
        default_dict = yaml.load(fp, Loader=yaml.FullLoader)
    return update_config(input_dict, default_dict)