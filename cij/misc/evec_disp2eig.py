import numpy
from numpy import newaxis as nax

def evec_disp2eig(a: numpy.ndarray, mass: list) -> numpy.ndarray:
    '''Convert phonon displacement vector to eigenvector.

    The displacement vectors (:math:`\\left|C_i\\right>`) are not orthogonal,
    whereas eigenvectors (:math:`\\left|C_i^m\\right>`) are orthogonal.
    The difference is the mass matrix (:math:`\\hat M`).
    
    :math:`\\left|C_i\\right>` and :math:`\\left|C_i^m\\right>` are related by

    .. math ::

        \\left|C_i^m\\right> = \\hat M^{1/2} \\left |C_i \\right>.

    This function returns normalized :math:`\\left|C_i^m\\right>`.


    :param a: 3N x 1 displacement vector or 3N x 3N displacement vector matrix
    :param mass: N x 1 atom mass vector, of any unit

    :returns: Eigenvector matrix, shape determined by the shape of ``a``:

        * If ``a`` is a displacement vector, return a 3N x 1 eigenvector
        * If ``a`` is a displacement matrix, return a 3N x 3N eigenvector matrix
    '''

    N = len(mass) # number of atoms

    m = numpy.repeat(mass, 3)   # Account for 3 polarization directions, N -> 3N
    a = numpy.copy(a)

    if a.shape == (3*N,):

        # Multiply by mass matrix

        a *= numpy.sqrt(m)

        # Renormalization

        norm = numpy.conj(a) @ a.T
        a /= numpy.sqrt(norm)

    elif a.shape == (3*N, 3*N):

        # Multiply by mass matrix

        a *= numpy.sqrt(m[nax, :])

        # Renormalization

        norm = numpy.diag(numpy.conj(a) @ a.T)
        a /= numpy.sqrt(norm)[:, nax]

    else:

        raise RuntimeError(
            f"Argument a ({a.shape}) or mass ({m.shape}) has wrong dimension!")


    return a