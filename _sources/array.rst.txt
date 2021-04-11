Conventions for using multi-demensional arrays
==============================================

Variables defined as a function of multiple variables, are stored as discrete
arrays in the form of ``numpy.ndarray`` internally. For the subscripts of such
variables, we follow the following conventions.

Conventions
------------

Variable on a grid of volume (pressure) and tempreature
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Index-conventions
   For an array of variable :math:`X` at
     - :math:`i`-th volume/pressure
     - :math:`j`-th tempreature

   we use
      X[i, j]

Examples
  - Helmholtz free energy :math:`F(V, T)`
  - Gibbs free energy :math:`G(P, T)`

Variable on a grid of volume, mode and q point
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Index-conventions
   For an array of variable :math:`X` at
     - :math:`i`-th volume
     - :math:`q`-th :math:`q`-point
     - :math:`m`-th phonon mode

   we use
        X[i, q, m]

Examples
  - Mode frequency :math:`\omega_{qm}(V)`
  - Grüneisen parameter :math:`\gamma_{qm}(V)`