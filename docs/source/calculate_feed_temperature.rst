First we calculate the feed temperature
(which is a saturated liquid) using
a bubble point calculation (see :ref:`bubble`).

.. code-block:: python:

    >>> model.calculate_T_feed()


We set the temperatures :math:`T_j`
to be the same as the feed temperature.
This step is performed automatically when the model
is initialized with input parameters.
We can check to see the feed temperature calculated
for our model by the following

.. code-block:: python

    >>> model.T_feed
    306.37018410667076

