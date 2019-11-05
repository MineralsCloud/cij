from ..io.traditional import models
import numpy
import scipy.interpolate, scipy.misc
import multiprocessing
import numba
from numpy.polynomial.polynomial import Polynomial

@numba.jit(parallel=True)
def interpolate_mode_spline(mode_volumes, mode_freqs, v_array, order=5):

    interp = scipy.interpolate.UnivariateSpline(
        numpy.flip(numpy.log(mode_volumes), axis=0),
        numpy.flip(numpy.log(mode_freqs), axis=0),
        k=order,
    )

    ln_v_array = numpy.log(v_array)

    return (
        numpy.exp(interp(ln_v_array)),
        interp(ln_v_array, nu=1),
        interp(ln_v_array, nu=2)
    )

@numba.jit(parallel=True)
def interpolate_mode_lagrange(mode_volumes, mode_freqs, v_array, order=6):

    # Lagrange interpolation is unstable when there is more than 6 nodes,
    # so pick at most 6 modes

    interval = int(numpy.ceil(mode_volumes.shape[0] / order))

    mode_volumes = mode_volumes[::interval]
    mode_freqs = mode_freqs[::interval]

    poly = scipy.interpolate.lagrange(
        numpy.flip(numpy.log(mode_volumes), axis=0),
        numpy.flip(numpy.log(mode_freqs), axis=0),
    )

    interp_lnfreq_lnv = poly
    interp_gamma_lnv =  numpy.polyder(poly, m=1)
    interp_vdr_dv_lnv = numpy.polyder(poly, m=2)

    ln_v_array = numpy.log(v_array)

    return (
        numpy.exp(interp_lnfreq_lnv(ln_v_array)),
        interp_gamma_lnv(ln_v_array),
        interp_vdr_dv_lnv(ln_v_array),
    )

@numba.jit(parallel=True)
def interpolate_mode_krogh(mode_volumes, mode_freqs, v_array, order=6):

    # Krogh interpolation is unstable when there is more than 12 nodes,
    # Krogh interpolation is unstable at ends for more than 8 nodes,
    # similar to lagrange interpolation at less than 6 nodes

    interval = int(numpy.ceil(mode_volumes.shape[0] / order))

    mode_volumes = mode_volumes[::interval]
    mode_freqs = mode_freqs[::interval]

    krogh = scipy.interpolate.KroghInterpolator(
        numpy.flip(numpy.log(mode_volumes), axis=0),
        numpy.flip(numpy.log(mode_freqs), axis=0),
    )

    ln_v_array = numpy.log(v_array)

    return (
        numpy.exp(krogh(ln_v_array)),
        krogh.derivative(ln_v_array, der=1),
        krogh.derivative(ln_v_array, der=2),
    )

@numba.jit(parallel=True)
def interpolate_mode_ppoly(mode_volumes, mode_freqs, v_array, method: str, order=6):

    # Krogh interpolation is unstable when there is more than 12 nodes,
    # Krogh interpolation is unstable at ends for more than 8 nodes,
    # similar to lagrange interpolation at less than 6 nodes

    interval = int(numpy.ceil(mode_volumes.shape[0] / order))

    mode_volumes = mode_volumes[::interval]
    mode_freqs = mode_freqs[::interval]

    if method == "pchip":
        Interpolator = scipy.interpolate.PchipInterpolator
    elif method == "akima":
        Interpolator = scipy.interpolate.Akima1DInterpolator
    elif method == "hermite":
        Interpolator = scipy.interpolate.CubicHermiteSpline

    interp = Interpolator(
        numpy.flip(numpy.log(mode_volumes), axis=0),
        numpy.flip(numpy.log(mode_freqs), axis=0),
    )

    ln_v_array = numpy.log(v_array)

    return (
        numpy.exp(interp(ln_v_array)),
        interp(ln_v_array, nu=1),
        interp(ln_v_array, nu=2),
    )

def interpolate_modes(
    qha_input: models.QHAInputData,
    v_array: numpy.ndarray,
    method: str = "spline",
    order = None
):

    nv = qha_input.nv
    nq = qha_input.nq
    np = qha_input.np

    ntv = v_array.shape[0]

    interp_freq = numpy.empty((ntv, nq, np))
    gamma_i= numpy.empty((ntv, nq, np))
    vdr_dv = numpy.empty((ntv, nq, np))

    mode_volumes = numpy.array([
        volume.volume
        for volume in qha_input.volumes
    ])

    for j in range(nq):
        for k in range(np):

            mode_freqs = numpy.array([
                volume.q_points[j].modes[k]
                for volume in qha_input.volumes
            ])

            if method == "lagrange":
                (   interp_freq[:, j, k],
                    gamma_i[:, j, k],
                    vdr_dv[:, j, k],
                ) = interpolate_mode_lagrange(
                    mode_volumes,
                    mode_freqs,
                    v_array,
                    order=order
                )
            elif method == "krogh":
                (   interp_freq[:, j, k],
                    gamma_i[:, j, k],
                    vdr_dv[:, j, k],
                ) = interpolate_mode_krogh(
                    mode_volumes,
                    mode_freqs,
                    v_array,
                    order=order
                )
            elif method in ["pchip", "hermite", "akima"]:
                (   interp_freq[:, j, k],
                    gamma_i[:, j, k],
                    vdr_dv[:, j, k],
                ) = interpolate_mode_ppoly(
                    mode_volumes,
                    mode_freqs,
                    v_array,
                    method=method,
                    order=order
                ) 
            elif method == "spline":
                (   interp_freq[:, j, k],
                    gamma_i[:, j, k],
                    vdr_dv[:, j, k],
                ) = interpolate_mode_spline(
                    mode_volumes,
                    mode_freqs,
                    v_array,
                    order=order
                )

    return interp_freq, gamma_i, vdr_dv
