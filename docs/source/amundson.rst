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

    from distillation.model.amundson_1958.main import Model

    model = Model(
        components=['n-Butane', 'n-Pentane'],
        F=1000.,              # feed flow rate [kmol/h]
        P=101325.*2,          # pressure (constant) [Pa]
        z_feed=[0.45, 0.55],  # mole fraction n-Butane = 0.45
        RR=1.,                # reflux ratio [L/D]
        D=500.,               # distillate flow rate [kmol/h]
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
In this case, the approach becomes something of the following:

.. code-block:: python

    for i in components:
        A = generate_A_matrix(i, *args)
        x = solve(A, b)

First, the tridiagonal matrix is generated.
Then, the efficient approach described in :ref:`tri-matrix` is used.
to solve the equations.

Behind the scenes, we save the current values of :math:`L_j`
and :math:`V_j` into the :code:`numpy` arrays
:meth:`distillation.amundson_1958.main.Model.L_old` and
:meth:`distillation.amundson_1958.main.Model.V_old`,
respectively, because we will need them
to determine the convergence of the flow rates
in :ref:`step-7`.

Step 5
******

In :ref:`step2`, we performed a :ref:`bubble` calculation
to determine the temperature of the feed.
We do the same thing here, except we do it multiple times
(i.e., for each stage).

.. _step-6:

Step 6, :math:`T_j` Convergence Check
*************************************
In Step 6, we determine whether the temperatures on all stages :math:`j`,
:math:`T_j`, have converged.
This is performed in :meth:`distillation.amundson_1958.main.Model.T_is_converged`

.. _step-7:

Step 7, Solve Energy Balances
*****************************
In Step 7, we solve the energy balances.
Here, all the balances can again be combined into one banded matrix.

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

.. _amund-code:

Class Method Reference
----------------------

.. autoclass:: distillation.amundson_1958.main.Model
    :members:
