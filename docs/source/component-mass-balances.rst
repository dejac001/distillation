For a general stage :math:`j`,
the mass balance on component :math:`i` is

.. math::

    V_jy_{i,j} + L_j x_{i,j} - V_{j+1}y_{i,j+1} - L_{j-1}x_{i,j-1} = F_jz_{i,j}

Using the following relationships

.. math::
    \begin{align}
        y_{i,j} &= K_jx_{i,j} \\
        y_{i,j+1} &= K_{j+1}x_{i,j+1} \\
    \end{align}

and

.. math::
    \begin{align}
        x_{i,j} &= \frac{l_{i,j}}{L_j}\\
        x_{i,j-1} &= \frac{l_{i,j-1}}{L_{j-1}}\\
    \end{align}

we obtain

.. math::

    -l_{i,j-1} + \left(1 + \frac{V_jK_j}{L_j}\right)l_{i,j} - \left(\frac{V_{j+1}K_{j+1}}{L_{j+1}}\right)l_{i,j+1}=F_jz_{i,j}

