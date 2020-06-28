from cij.io.traditional.qha_input import *
from .evec_load import evec_load
from .evec_sort import evec_sort

import multiprocessing.pool
import numpy

import re

def _evec_load(params):
    return evec_load(*params)

def regen_freq(input_data: QHAInputData, eig_files: List[str]):
    '''Rebuild QHA input data with the reindexed mode frequencies

    :input_data: The QHA input data object.
    :eig_files: The name of the eigenvector files for each volumes in the QHA
        input_data object.
    '''


    if input_data.nv != len(eig_files): raise RuntimeError("Inconsistent nv")

    evecs = []

    p = multiprocessing.pool.Pool(12)
    evecs = p.map(_evec_load, [(fn, input_data.nq, input_data.np) for fn in eig_files])

    _evecs = dict()
    _freqs = dict()

    for i in range(input_data.nv):

        for k in range(input_data.nq): # q-points


            if i == 0:
                _sorted = evecs[i][k][1]
            else:
                _prev_evecs = _evecs[i-1,k]
                _curr_evecs = numpy.array([evecs[i][k][1][l][1] for l in range(input_data.np)])
                _sorted = evec_sort(evecs[i][k][1], _curr_evecs, _prev_evecs)

            _evecs[i,k] = [_sorted[l][1]    for l in range(input_data.np)]
            _freqs[i,k] = [_sorted[l][0][2] for l in range(input_data.np)]

    # Rebuild input data

    _input_data = QHAInputData(
        input_data.nv,
        input_data.nq,
        input_data.np,
        input_data.nm,
        input_data.na,
        input_data.weights,
        []
    )

    for i in range(input_data.nv):

        volume = VolumeData(
            input_data.volumes[i].pressure,
            input_data.volumes[i].volume,
            input_data.volumes[i].energy,
            []
        )

        _input_data.volumes.append(volume)

        for k in range(input_data.nq):

            volume.q_points.append(QPointData(
                evecs[i][k][0],
                _freqs[i,k]
            ))

    return _input_data