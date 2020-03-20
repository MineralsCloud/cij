from pathlib import Path
from .validate import validate_config

def read_config(fname: str):
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
    validate_config(config)
    return config
