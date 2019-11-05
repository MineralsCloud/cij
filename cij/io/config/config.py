import yaml

def read_config(fname: str):
    with open(fname) as fp:
        return yaml.load(fp, Loader=yaml.FullLoader)