from pathlib import Path
import sys

__path__ = Path(__file__).absolute().parent

def get_data_fname(fname: str) -> str:
    '''
    Get the full path of the data file from its relative path.
    '''
    if sys.version_info >= (3, 9):
        # Modern approach using importlib.resources
        from importlib.resources import files
        resource_path = files(__name__) / fname
        # For Python 3.9+, we need to use as_file context manager
        # but for simple path access, we can convert to string
        return str(resource_path)
    else:
        # Fallback for Python 3.8 using importlib_resources backport or traversable
        try:
            from importlib.resources import files
            resource_path = files(__name__) / fname
            return str(resource_path)
        except ImportError:
            # Ultimate fallback to simple path joining
            return str(__path__ / fname)

