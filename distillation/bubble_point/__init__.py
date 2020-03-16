"""The bubble-point temperature is the temperature of stage *j* at which
the liquid mixture begins to boil.

Here, the pressure *p_j* and the mole fractions of the liquid, *x_{i,j}*
will be known. We wish to find the temperature, *T_j*,
at which

.. math::
    \sum_i y_{i,j} = 1.0


The gas-phase mole fractions are calculated as

.. math::
    y_{i,j} = K_{i,j}(T_j) x_{i,j}
"""