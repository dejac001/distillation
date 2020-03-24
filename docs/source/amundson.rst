Amundson_1958
=============
This code follows the approach that was published by Amundson and Pontinen [Amundson1958]_.

.. [Amundson1958] Amundson, NR and Pontinen, AJ. Multicomponent Distillation Calculations on a Large Digital Computer. *Ind. Eng. Chem.* 1958;50:730--736.

This approach takes advantage of the following mathematical trick.

.. _tri-matrix:

Mathematical Trick
------------------

A system of equations such as

.. math::

    Ax=b

where :math:`A` is a matrix of constants,
:math:`x` is a vector of unknowns,
and :math:`b` is a vector of constants,
can be solved very efficiently if the matrix
:math:`A` can be made into the following form:

.. math::

    A = \begin{bmatrix}
        d_0 & u_0 & 0 & 0 & 0  \\
        l_0 & d_1 & u_1 & 0 & 0 \\
        0 & l_1 & d_2 & u_2 & 0  \\
        0 & 0 & l_2 & d_3 & u_3  \\
        0 & 0 &  0 & l_3 & d_4
    \end{bmatrix}

This special type of matrix
is called a tridiagonal_matrix_.

.. _tridiagonal_matrix: https://en.wikipedia.org/wiki/Tridiagonal_matrix

It can be solved, for example, if we define the following vectors


.. math::

    \begin{align}
        l &= [l_0, l_1, l_2, l_3] \\
        d &= [d_0, d_1, d_2, d_3, d_4] \\
        u &= [u_0, u_1, u_2, u_3] \\
        b &= [b_0, b_1, b_2, b_3]
    \end{align}

or, with python code, like

.. code-block:: python

   l = [l_0, l_1, l_2, l_3]
   d = [d_0, d_1, d_2, d_3, d_4]
   u = [u_0, u_1, u_2, u_3]
   b = [b_0, b_1, b_2, b_3]


the :math:`5\times5` matrix can be solved
efficiently using the following code:

.. code-block:: python

    from scipy.sparse import linalg, diags
    A = diags(
        diagonals=[l, d, u],
        offsets=[-1, 0, 1],
        shape=(5, 5),
    )
    x = linalg.spsolve(A, b)

This mathematical trick is used throughout
the solution algorithm.

The algorithm
-------------

The complete algorithm is shown in the diagram below.

.. graphviz::

    digraph {
        input -> initial_guess -> calc_K -> mass_balance -> bubble_pt -> converg_1;
        converg_1 -> calc_K[label="No",color=red];
        converg_1 -> energy_balance[label="Yes",color=green];
        energy_balance -> converg_2;
        converg_2 -> finish[label="Yes",color=green];
        converg_2 -> mass_balance[label="No",color=red];

        input [shape=polygon,side=4,label="1. Input equilibrium and enthalpy data.\n Input design variables"];
        initial_guess [shape=polygon,side=4,label="2. Pick L, V, T\n on every stage"];
        calc_K [shape=polygon,side=4,label="3. Calculate K for each component\n on each stage"];
        mass_balance [shape=polygon,side=4,label="4. Solve component mass\n balances in matrix form",color=lightblue,style=filled];
        bubble_pt [shape=polygon,side=4,label="5. Calculate T on each stage\n using bubble point calculation"];
        converg_1 [shape=diamond,label="6. Determine whether\n T is converged\n for all stages"];
        energy_balance [shape=polygon,side=4,label="7. Use energy balance to calculate\n L, V on every stage",color=lightblue,style=filled];
        converg_2 [shape=diamond,label="8. Determine whether\n L, V are converged\n for all stages"];
        finish [shape=ellipse,label="9. Finished"];
    }


Tutorial
--------
Step 1, Input
*************

In Step 1, the property data and design variables are input.
The property data is described in :ref:`equil-data`.
A model with this property data and input design variables is created
by the following code:

.. code-block:: python

    >>> from distillation.amundson_1958 import Model

    >>> model = Model(
            components=['n-Butane', 'n-Pentane'],
            F=1000.,              # feed flow rate [kmol/h]
            P=101325.*2,          # pressure (constant) [Pa]
            z_feed=[0.45, 0.55],  # mole fraction n-Butane = 0.45
            RR=1.,                # reflux ratio [L/D]
            D=400.,               # distillate flow rate [kmol/h]
            N=3,                  # number of equilibrium contacts
            feed_stage=2,         # stage at which feed goes in
    )


Currently, the input design flow rates must be
specified in kmol/h and the pressure must be specified in Pa.
The method called in this process is
:py:meth:`distillation.amundson_1958.main.Model.__init__`
(see :ref:`amund-code`).


.. _step2:

Step 2, Make Initial Guess
**************************
.. include:: step2.rst

Step 3, Update :math:`K`-values
*******************************

.. include:: step3.rst

.. _stage-4:

Step 4, Solve Component Mass Balances
*************************************

The two stages in which tridiagonal matrices are used
for computational efficiency are :ref:`stage-4` and :ref:`step-7`,
which are depicted in colored blue boxes.

In this step, the matrices are tridiagonal for the complete
mass balances for *each* component (as opposed to all components together).
The variables that are calculated are :math:`l_{ij}` for each component :math:`i`
and each stage :math:`j`.
This new variable is introduced where :math:`l_{ij}=x_iL_j`,
and can be interpreted as a *component* flow rate.
By formulating the mass balance equations into function of :math:`l_{ij}`
the system of equations can be converted in to a tridiagonal matrix.

.. include:: component-mass-balances.rst

Computationally, we can perform this step
using our model by the following:

.. code-block:: python

    >>> for i in model.components:
    ...    model.solve_component_mass_bal(i)
    ...
    >>>


First, the tridiagonal matrix is generated.
Then, the efficient approach described in :ref:`tri-matrix` is used.
to solve the equations.

We can see the results of the calculation by probing
the model attributes as follows

.. code-block:: python

    >>> import pprint # pretty printing
    >>> pprint.pprint(model.l)
    {'n-Butane': array([288.88918652, 179.07801532, 507.65007098, 161.11081348]),
     'n-Pentane': array([ 74.88288594, 150.28018773, 790.77762525, 475.11711406])}

Since :math:`\sum_i l_{i,j} = L_j`, it is worthwhile to compare
our calculations to what was calculated for :math:`L_j` with the
Lewis method (i.e. assuming Constant Molal Overflow)

.. code-block:: python

    >>> for stage in model.stages:
    ...     print(stage, model.L[stage], model.l['n-Butane'][stage]+model.l['n-Pentane'][stage])
    ...
    0 400.0 363.7720724598715
    1 400.0 329.3582030449366
    2 1400.0 1298.4276962288486
    3 600.0 636.2279275401286
    >>>

From this analysis, we realize that our initial guess for the liquid and vapor flow rates was not that bad.

.. _step-5:

Step 5
******

In :ref:`step2`, we performed a :ref:`bubble` calculation
to determine the temperature of the feed.
We do the same thing here, except we do it multiple times
(i.e., for each stage).
However, we have not calculated the liquid phase mole fractions
explicitly yet.
We can use the liquid component mass balances calculated in :ref:`stage-4`
to calculate the mole fractions as

.. math::

    x_{ij} = \frac{l_{ij}}{\sum_k l_{kj}}

and then perform the :ref:`bubble` calculations
with these mole fractions.

Computationally, we can perform this step with our
model using the following

.. code-block:: python

    >>> model.T
    array([306.37018411, 306.37018411, 306.37018411, 306.37018411])
    >>> model.update_T_values()
    >>> model.T
    array([295.34137166 302.94861315 308.72858323 314.97039634])

And we can see how the temperatures change from the initial guess (the feed temperature)
to the results of the first calculation.

.. _step-6:

Step 6, :math:`T_j` Convergence Check
*************************************

.. include:: temp-converge.rst

We can see from :ref:`step-5` that the temperature
has clearly changed by more than 0.01~K between
the first iteration. When we perform this step,
we find that the temperature has not converged.

.. code-block:: python

    >>> model.T_is_converged()
    False

Following the algorithmic diagram,
we go back to Step 3 and perform a few more iterations.

Its convenient to perform these iterations with a while loop, as follows

.. code-block:: python

    while not model.T_is_converged():
    ...     model.update_K_values()
    ...     for i in model.components:
    ...         model.solve_component_mass_bal(i)
    ...     model.update_T_values()
    ...     print(model.T)
    ...
    [295.50131711 303.41982848 309.35773209 315.96224522]
    [295.6166139  303.62102821 309.52271272 316.15309479]
    [295.65094925 303.67756389 309.56508571 316.19923696]
    [295.66006297 303.69238402 309.57598554 316.21097278]
    [295.66242518 303.69621481 309.57879188 316.21398756]

The temperature is converged after 5 more iterations,
as we can see from the output above.
Since the temperature has converged,
we can now proceed to the next step.

.. _step-7:

Step 7, Solve Energy Balances
*****************************
In Step 7, we solve the energy balances.
Here, all the balances can again be combined into one banded matrix.

.. include:: energy-bal.rst

A helper function for creation and solution of these
matrices is provided in the model.
In other words, the energy balances can be solved
by the following code:

.. code-block:: python

    >>> model.L
    array([ 400.,  400., 1400.,  600.])
    >>> model.V
    array([  0., 800., 800., 800.])
    >>> model.solve_energy_balances()
    >>> model.L
    array([ 400.        ,  312.3605968 , 1300.10206406,  600.        ])
    >>> model.V
    array([  0.        , 730.64500782, 712.3605968 , 700.10206406])
    >>>

Where we have added extra steps so that we can check the
values of :code:`model.L` and :code:`model.V` before
and after solving the energy balances.

Step 8, :math:`V_j`, :math:`L_j` Convergence Check
**************************************************
In this step, we determine if the simulation has converged.
If the following holds true for all stages :math:`j`

.. math::

    \sqrt{\left(\frac{X_{j,\mathrm{new}} - X_{j,\mathrm{old}}}{X_{j,\mathrm{new}}}\right)^2} < \epsilon

for each variable :math:`X=V` and :math:`X=L`.
The flow rate tolerance, :math:`\epsilon`,
is found in the attribute :attr:`distillation.amundson_1958.main.Model.flow_rate_tol`.
The code for this step is found in :meth:`distillation.amundson_1958.main.Model.flow_rates_converged`

Obviously, our flow rates have not converged after just one iteration,
as we have seen that several values have changed by :math:`\\approx 100` kmol/h.
We will need to come up with a looping algorithm to
get the simulation to converge, as described in the next section.

Putting it all together
-----------------------
So far, we can combine our code from the tutorial into
the following

.. code-block:: python

    >>> from distillation.amundson_1958 import Model

    >>> model = Model(
            components=['n-Butane', 'n-Pentane'],
            F=1000.,              # feed flow rate [kmol/h]
            P=101325.*2,          # pressure (constant) [Pa]
            z_feed=[0.45, 0.55],  # mole fraction n-Butane = 0.45
            RR=1.,                # reflux ratio [L/D]
            D=400.,               # distillate flow rate [kmol/h]
            N=3,                  # number of equilibrium contacts
            feed_stage=2,         # stage at which feed goes in
    )
    >>> model.generate_initial_guess()
    >>> while not model.T_is_converged():
    ...     model.update_K_values()
    ...     for i in model.components:
    ...         model.solve_component_mass_bal(i)
    ...     model.update_T_values()
    ...
    >>> model.solve_energy_balances()

But we need to figure out how to do the last iteration to automatically try to converge.
It turns out, we can do this with a nested :code:`while` loop,
as shown below

.. code-block:: python

    >>> while not model.flow_rates_converged():
    ...     for i in self.components:
    ...         model.solve_component_mass_bal(i)
    ...     model.update_T_values()
    ...     while not model.T_is_converged():
    ...         model.update_K_values()
    ...         for i in model.components:
    ...             model.solve_component_mass_bal(i)
    ...         model.update_T_values()
    ...     model.solve_energy_balances()
    ...

This is what happens when the user invokes :code:`model.run()`
For our example, it turns out that it takes about 4 loops
through the step 4-8 loop, and then we will be converged.
The final results for the flow rates can be calculated as

.. code-block:: python

    >>> model.L
    array([ 400.        ,  353.80944637, 1340.92773179,  600.        ])
    >>> model.V
    array([  0.        , 733.53885641, 753.80944637, 740.92773179])

which is still quite similar to what we calculated by the Lewis method.
We can also take a look at the stage temperatures and mole fractions
as below

.. code-block:: python

    >>> model.T
    array([295.60405689, 303.58875055, 309.37291267, 315.8254746 ])
    >>> list(model.x_ij_expr('n-Butane', i) for i in model.stages)
    [0.7648435805973592, 0.5587863474554224, 0.38180690272177853, 0.24010427960176078]
    >>> list(model.y_ij_expr('n-Butane', i) for i in model.stages)
    [0.902879383384724, 0.8342045394354831, 0.6681268674620511, 0.4965317369291433]

As expected, the distillate is greatly enriched in the light component (n-Butane).

.. _amund-code:

Class Method Reference
----------------------

.. autoclass:: distillation.amundson_1958.main.Model
    :members:
