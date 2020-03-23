In Step 3, we calculate :math:`K_{i,j}`, the
:math:`K`-value for each component :math:`i` on stage :math:`j`.
We can do this now for our model by the following:

.. code-block:: python

    >>> model.update_K_values()

We can think of this step as doing something similar to the following:

.. code-block:: python

    for i in model.components:
        for j in range(model.N + 1):
            model.K[i][j] = K_func(i, T[j], model.P_feed)


where the :code:`K_func` function returns the :math:`K`-value for component
:math:`i` on stage :math:`j`.
The :code:`for` loops iteratively store the :math:`K`-values in each component of the amtrix.

In our example, we can checkout our new :math:`K`-values as follows:

.. code-block::

    >>> import pprint  # for pretty printing
    >>> pprint.pprint(model.K)
    {'n-Butane': array([1.61320297, 1.61320297, 1.61320297, 1.61320297]),
     'n-Pentane': array([0.49828848, 0.49828848, 0.49828848, 0.49828848])}

As expected, the values for butane are higher for the light component, butane.
The values are also the same for all stages because the temperature
of all stages are currently fixed at the feed temperature.

Behind the scenes, we then store the current stage temperatures
in an attribute :attr:`~distillation.amundson_1958.main.T_old` which
we will need later in :ref:`step-6` to determine whether the temperature
has converged.

The code for this step is depicted in
:meth:`distillation.amundson_1958.main.Model.update_K_values`.
