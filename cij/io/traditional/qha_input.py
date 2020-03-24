from typing import List, Tuple, NamedTuple
import re


class QPointData(NamedTuple):
    coord: Tuple[float, float, float] = []
    modes: List[float] = []


class VolumeData(NamedTuple):
    pressure: float
    volume: float
    energy: float
    q_points: List[QPointData] = []


class QPointWeight(NamedTuple):
    coord: Tuple[float, float, float]
    weight: float


class QHAInputData(NamedTuple):
    nv: int
    nq: int
    np: int
    nm: int
    na: int
    weights: List[Tuple[Tuple[float, float, float], float]] = []
    volumes: List[VolumeData] = []


REGEX_INFO_START = r"^(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)$"
REGEX_Q_WEIGHT = r"^(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)$"
REGEX_PVE = r"P=\s+(-?\d*\.?\d*)\s+V=\s+(-?\d*\.?\d*)\s+E=\s+(-?\d*\.?\d*)"

def _read_weights(lines, nq):
    def _yield_weights():
        for _ in range(nq):
            line = next(lines)
            res = re.search(REGEX_Q_WEIGHT, line.strip())
            yield QPointWeight(tuple(map(float, res.groups()[:-1])), float(res.group(4)))
    return list(_yield_weights())

def _read_volume_data(lines, nv, nq, np):

    def _yield_mode_data():
        for _ in range(np):
            yield float(next(lines))

    def _yield_q_point_data():
        for _ in range(nq):
            coord = tuple(map(float, next(lines).strip().split()))
            yield QPointData(coord, list(_yield_mode_data()))

    def _yield_volume_data():
        for _ in range(nv):
            P, V, E = map(float, re.search(REGEX_PVE, next(lines)).groups())
            yield VolumeData(P, V, E, list(_yield_q_point_data()))

    return list(_yield_volume_data())

def read_energy(fname: str) -> QHAInputData:
    '''Read :math:`E`, :math:`V`, :math:`\omega` etc., from QHA data file

    :param fname: The path of the input file
    '''

    with open(fname, encoding="utf8") as fp:
        for line in fp:
            res = re.search(REGEX_INFO_START, line.strip())
            if res:
                (nv, nq, np, nm, na) = map(int, res.groups())
                break

        qha_input_data = QHAInputData(nv, nq, np, nm, na)
 
        for volume in _read_volume_data(fp, nv, nq, np):
            qha_input_data.volumes.append(volume)

        for line in fp:
            if line.strip() == "weight":
                break

        for weight in _read_weights(fp, qha_input_data.nq):
            qha_input_data.weights.append(weight)

    return qha_input_data

