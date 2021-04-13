from typing import List, Tuple, NamedTuple
import re


class QPointData(NamedTuple):
    coord: Tuple[float, float, float]
    modes: List[float]

class VolumeData(NamedTuple):
    pressure: float
    volume: float
    energy: float
    q_points: List[QPointData]


class QPointWeight(NamedTuple):
    coord: Tuple[float, float, float]
    weight: float


class QHAInputData(NamedTuple):
    nv: int
    nq: int
    np: int
    nm: int
    na: int
    weights: List[Tuple[Tuple[float, float, float], float]]
    volumes: List[VolumeData] = None


REGEX_INFO_START = r"^(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)$"
# REGEX_Q_WEIGHT = r"^(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)$"
REGEX_PVE = r"P=\s+(-?\d*\.?\d*)\s+V=\s+(-?\d*\.?\d*)\s+E=\s+(-?\d*\.?\d*)"
REGEX_PVE = r"\S=\s+(\S+)\s+\S=\s+(\S+)\s+\S=\s+(\S+)"

def _read_weights(lines, nq):
    def _yield_weights():
        for _ in range(nq):
            line = next(lines)
            words = line.strip().split()
            # res = re.search(REGEX_Q_WEIGHT, line.strip())
            yield QPointWeight(tuple(map(float, words[0:3])), float(words[3]))
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
            while True:
                line = next(lines)
                if line.strip() == "": continue
                break
            P, V, E = map(float, re.search(REGEX_PVE, line).groups())
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

        qha_input_data = QHAInputData(nv, nq, np, nm, na, [], [])
 
        for volume in _read_volume_data(fp, nv, nq, np):
            qha_input_data.volumes.append(volume)

        for line in fp:
            if line.strip() in ["weight", "weights"]:
                break

        for weight in _read_weights(fp, qha_input_data.nq):
            qha_input_data.weights.append(weight)

    return qha_input_data


def write_energy(fname: str, input_data: QHAInputData, comment: str = "QHA Input data"):
    '''Write QHA input data to file

    :param fname:
    :param input_data:
    :param comment: The comment line
    '''
    lines = []

    lines.append(comment)
    lines.append("")
    lines.append(" ".join("%s".rjust(4) % x for x in ('nv', 'nq', 'np', 'nm', 'na')))
    lines.append(" ".join("%4d" % x for x in (
        input_data.nv,
        input_data.nq,
        input_data.np,
        input_data.nm,
        input_data.na)))
    lines.append("")

    for p, v, e, v_data in input_data.volumes:
        lines.append(f"P= {p:12.6f} V= {v:12.6f} E= {e:12.6f}")
        for coords, modes in v_data:
            lines.append(" ".join("%10.4f" % c for c in coords))
            for cm_1 in modes:
                lines.append(f"{cm_1:12.6f}")

    lines.append("")
    lines.append("weight")
    for coords, weight in input_data.weights:
        lines.append("%10.6f %10.6f %10.6f %10.6f" % (*coords, weight))

    with open(fname, "w", encoding="utf8") as fp:
        fp.writelines([f"{line}\n" for line in lines])