Before beginning the solution procedure,
we need to have an initial guess for the liquid flow rates,
vapor flow rates, and temperatures on each stage :math:`j`.
That is, we need to generate initial guesses for
:math:`L_j`, :math:`V_j`, and :math:`T_j`.

.. include:: calculate_feed_temperature.rst

Then, we generate the initial guesses
for these values by assuming constant molal overflow (the Lewis method).
Computationally, we can do this using the following

.. code-block:: python

    >>> model.generate_initial_guess()

We can then check the arrays of liquid flow rates (:math:`L`)
and vapor flow rates (:math:`V`) obtained by the calculations as below

.. code-block:: python

    >>> model.L
    array([ 400.,  400., 1400.,  600.])
    >>> model.V
    array([   0., 800., 800., 800.])

or, looking at the liquid flow rates by stage more specifically,
.. code-block:: python

    >>> for stage, L_j in enumerate(model.L):
    ...    print('L_j', stage, L_j)
    ...
    L_j 0 400.0
    L_j 1 400.0
    L_j 2 1400.0
    L_j 3 600.0

Here, we see that the liquid flow rates in the enriching
section are :math:`(L/D)\times D=1\times 400=400=L`.
Since the feed stage is stage 2, we notice that
the liquid leaving this stage is equal to :math:`\overline{L}=L + F=400 + 1000`.
Finally, we notice that the liquid leaving the bottom stage,
the partial reboiler, is the same as the bottoms
flow rate that would be calculated from an overall balance,
:math:`B = F - D = 1000 - 400 = 600`.

.. code-block:: python

    >>> for stage, V_j in enumerate(model.V):
    ...    print('V_j', stage, V_j)
    ...
    V_j 0 0.0
    V_j 1 800.0
    V_j 2 800.0
    V_j 3 800.0

We note that there is no vapor leaving stage 0 because the
column is equipped with a total condensor.
We note that :code:`V[1] = model.RR*model.D + model.D`,
which arises from the mass balance at the top of the column.
Similarly, we also realize that
:code:`V[3] = model.L[2] - model.B`.
