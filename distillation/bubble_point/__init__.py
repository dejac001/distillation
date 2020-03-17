"""The bubble-point temperature is the temperature of stage *j* at which
the liquid mixture begins to boil.

Here, the pressure :math:`p` and the mole fractions of the liquid, :math:`x_i`
will be known. We wish to find the temperature, :math:`T`,
at which

.. math::

    \sum_i y_i = 1.0


The gas-phase mole fractions are calculated as

.. math::

    y_i = K_i(T) x_i

Where :math:`K_i` comes from the equilibrium relationship that has been chosen.

The calculation results in finding
the temperature that satisfies the stoichiometric equation

.. math::

    \sum_i y_i = \sum_i K_i (T) x_i = 1.0
"""