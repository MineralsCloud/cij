import numpy

from qha.fitting import apply_finite_strain_fitting, polynomial_least_square_fitting
from qha.grid_interpolation import calculate_eulerian_strain, from_eulerian_strain
from scipy.interpolate import InterpolatedUnivariateSpline


def get_static_p_of_v(v_sparse, f_sparse, v0 = None, N = 100, order = 3, v_ratio = 1.2) -> callable:

    v0 = v0 if v0 != None else v_sparse[numpy.argmin(f_sparse)]
    x_sparse = calculate_eulerian_strain(v0, v_sparse)
    x_dense = numpy.linspace(
        calculate_eulerian_strain(v0, numpy.max(v_sparse) * v_ratio),
        calculate_eulerian_strain(v0, numpy.min(v_sparse) / v_ratio),
        N + 1)
    _, f_dense = polynomial_least_square_fitting(x_sparse, f_sparse, x_dense, order = order)
    v_dense = from_eulerian_strain(v0, x_dense)
    p_dense = -numpy.gradient(f_dense) / numpy.gradient(v_dense)

    def static_p_of_v(v_new: numpy.array) -> numpy.array:
        x_new = calculate_eulerian_strain(v0, v_new)
        return InterpolatedUnivariateSpline(x_dense, p_dense)(x_new)

    return static_p_of_v

if __name__ == "__main__":

    import sys

    from cij.io.traditional.qha_input import read_energy
    from cij.util.units import _to_gpa

    qha_input = read_energy(sys.argv[1])
    v_static, f_static = zip(*[(v_data.volume, v_data.energy) for v_data in qha_input.volumes])
    static_p_of_v = get_static_p_of_v(v_static, f_static)
    p_static = _to_gpa(static_p_of_v(v_static)) * 10

    print("{0:>15} {1:>15} {2:>15}".format("kbar", "bohr^3", "rydberg"))
    for p, v, e in zip(p_static, v_static, f_static):
        print(f"P= {p:12.6f} V= {v:12.6f} E= {e:12.6f}")

