import re
import io
import itertools
from logging import Logger

logger = Logger(__file__)

Q_COORDS_REGEX = re.compile(r"q\s*=\s*(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)")
MODE_INDEX_REGEX = re.compile(r"freq\s*\(\s*(\d+)\)\s*=\s*(-?\d+\.?\d*)\s*\[THz\]\s*=\s*(-?\d+\.?\d*)\s*\[cm\-1\]")

def _read_vecs(fp: io.StringIO, np: int):

    for m in range(np // 3):
        line = next(fp).strip()
        yield (
            float(line[ 2:12]) + float(line[13:23]) * 1j,
            float(line[26:36]) + float(line[37:47]) * 1j,
            float(line[50:60]) + float(line[61:71]) * 1j,
        )

def _read_modes(fp: io.StringIO, np: int):

    for l in range(np):

        line = next(fp).strip()
        res = MODE_INDEX_REGEX.search(line.strip())
        mode_id, thz, cm_1 = tuple((
            func(string) for func, string in 
            zip((int, float, float), res.groups())
        ))

        yield ((mode_id, thz, cm_1), tuple(itertools.chain.from_iterable(_read_vecs(fp, np))))

def _read_q_points(fp: io.StringIO, nq: int, np: int):

    for k in range(nq):

        for _ in range(2): next(fp)

        line = next(fp).strip()

        res = Q_COORDS_REGEX.search(line)
        q_coords = tuple([float(i) for i in res.groups()])

        next(fp)

        yield (q_coords, tuple(_read_modes(fp, np)))

        next(fp)


def evec_load(fname: str, nq: int, np: int) -> list:
    '''Load ``eig`` file from ``matdyn.x`` output.

    :param fname: The name of the ``eig`` file
    :param nq: The number of q points
    :param np: The number of modes
    '''
    with open(fname) as fp:
        return list(_read_q_points(fp, nq, np))

if __name__ == "__main__":
    import sys
    print(evec_load(sys.argv[1], int(sys.argv[2]), int(sys.argv[3])))