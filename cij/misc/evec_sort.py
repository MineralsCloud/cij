import numpy

def evec_sort(target_arr, target_evecs, base_evecs, filter = None, threshold: float = None):

    ndim = len(target_arr)
    sorted_arr = [None] * ndim

    s = set([len(target_evecs), len(base_evecs), *[len(i) for i in (target_evecs + base_evecs)]])
    if len(s) != 1 or ndim not in s:
        raise RuntimeError(f"Wrong input dimension of evec, they should be of {ndim} x {ndim}")

    m = numpy.conj(numpy.array(base_evecs)) @ numpy.array(target_evecs).T
    if filter: m = filter(m)
    for i in range(ndim):
        idx = numpy.unravel_index(numpy.argmax(m), m.shape)
        if threshold and m[idx] < threshold:
            raise RuntimeError("Low eval product")
        m[idx[0], :] = 0
        m[:, idx[1]] = 0
        sorted_arr[idx[0]] = target_arr[idx[1]]

    return sorted_arr
