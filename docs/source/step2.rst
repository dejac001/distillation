
Before beginning the solution procedure,
we need to have an initial guess for the liquid flow rates,
vapor flow rates, and temperatures on each stage :math:`j`.
That is, we need to generate initial guesses for
:math:`L_j`, :math:`V_j`, and :math:`T_j`.

First we calculate the feed temperature
(which is a saturated liquid) using
a bubble point calculation (see :ref:`bubble`).
We set the temperatures :math:`T_j`
to be the same as the feed temperature.

Then, we generate the initial guesses
for these values by assuming constant molal overflow.

The code for this step
is depicted in :meth:`distillation.amundson_1958.main.Model.generate_initial_guess`.
