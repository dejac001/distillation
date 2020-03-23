In this step, we determine whether the temperatures on all stages :math:`j`,
:math:`T_j`, have converged.
Mathematically, we require the following

.. math::

    \sqrt{\left(T_{j,\mathrm{new}} - T_{j,\mathrm{old}}\right)^2} < \epsilon

The temperature tolerance :math:`\epsilon`
is attribute :attr:`distillation.amundson_1958.main.Model.temperature_tol`;
It is currently set to 0.01.
