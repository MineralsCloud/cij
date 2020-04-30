from pathlib import Path
import pkg_resources

__path__ = Path(__file__).absolute().parent

def get_data_fname(fname) -> str:
    return pkg_resources.resource_filename(__name__, fname)

