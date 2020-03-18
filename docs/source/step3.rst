In Step 3, we calculate :math:`K_{i,j}`, the
:math:`K`-value for each component :math:`i` on stage :math:`j`.
We can think of this step as doing something similar to the following:

.. code-block:: python

    for i in components:
        for j in stages:
            K[i,j] = K_func(i, T[j], p)


where the :code:`K_func` function returns the :math:`K`-value for component
:math:`i` on stage :math:`j`.
The :code:`for` loops iteratively store the :math:`K`-values in each component of the amtrix.

Behind the scenes, we then store the current stage temperatures
in an attribute :attr:`~distillation.amundson_1958.main.T_old` which
we will need later in :ref:`step-6` to determine whether the temperature
has converged.

The code for this step is depicted in
:meth:`distillation.amundson_1958.main.Model.update_K_values`.
