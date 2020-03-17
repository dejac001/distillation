.. image:: column_diagram.png
    :width: 600

.. image:: stage_balance_diagram.png
    :width: 400


Mathematical Tricks
===================
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

For example, this :math:`5\times5` matrix can be solved
efficiently using the following code:

.. code-block:: python

    from scipy.sparse import linalg, diags
    A = diags(
        diagonals=[l, d, u],
        offsets=[-1, 0, 1],
        shape=(5, 5),
    )
    x = linalg.spsolve(A, b)


where

.. math::

    \begin{align}
        l &= [l_0, l_1, l_2, l_3] \\
        d &= [d_0, d_1, d_2, d_3, d_4] \\
        u &= [u_0, u_1, u_2, u_3] \\
        b &= [b_0, b_1, b_2, b_3]
    \end{align}
