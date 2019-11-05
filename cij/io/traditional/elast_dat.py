from typing import List, Tuple, NamedTuple
import re

class ElastVolumeData(NamedTuple):
    volume: float
    static_elastic_coeficient: dict

class ElastData(NamedTuple):
    vref: float
    nv: int
    cellmass: float
    volumes: List[ElastVolumeData] = []

def read_elast_data(fname) -> ElastData:
    with open(fname, encoding="utf8") as fp:
        next(fp)
        fields = next(fp).strip().split()
        vref = float(fields[0])
        nv = int(fields[1])
        cellmass = float(fields[2])
        ret = ElastData(vref, nv, cellmass)

        keys = next(fp).strip().split()

        for line in fp:
            fields = tuple(map(float, next(fp).strip().split()))
            ret.volumes.append(ElastVolumeData(fields[0], dict(zip(keys[1:], fields[1:]))))

    return ret