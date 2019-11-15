from typing import List, Tuple, NamedTuple
from cij.util import c_
import re

class ElastVolumeData(NamedTuple):
    volume: float
    static_elastic_modulus: dict

class ElastData(NamedTuple):
    vref: float
    nv: int
    cellmass: float
    volumes: List[ElastVolumeData] = []

REGEX_MODULUS = r"^\D*(\d+)$"

def _find_modulus_key(key: str):
    res = re.search(REGEX_MODULUS, key)
    if res:
        return c_(res.group(1))
    else:
        return key

def read_elast_data(fname) -> ElastData:
    with open(fname, encoding="utf8") as fp:
        next(fp)
        fields = next(fp).strip().split()
        vref = float(fields[0])
        nv = int(fields[1])
        cellmass = float(fields[2])
        ret = ElastData(vref, nv, cellmass)

        keys = next(fp).strip().split()
        keys = [_find_modulus_key(x) for x in keys]

        for _ in range(nv):
            line = fp.readline()
            fields = tuple(map(float, line.strip().split()))
            ret.volumes.append(ElastVolumeData(fields[0], dict(zip(keys[1:], fields[1:]))))

    return ret