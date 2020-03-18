.. _equil-data:

Equilibrium Properties
======================

The equilibrium behavior
for each component :math:`i` are specified using the following :math:`K`-value relationships

.. math::

    y_i = K_i(T,p)x_i

where :math:`y_i` is the vapor phase mole fraction
and :math:`x_i` is the liquid phase mole fraction.
Here, we have assumed for simplicity that :math:`K_i` does not depend on composition.

In order to calculate the enthalpies of the liquid and vapor streams,
heat capacities and heats of vaporization are required.

The heats of vaporization and heat capacities of the liquid
are obtained from [Perry]_.
The equilibrium properties are taken from
the DePriester charts in [Wankat]_.
The heat capacities of the vapor are currently estimated as a polyatomic ideal gas.

.. [Perry] Green, Don W., and Robert H. Perry. Perry's Chemical Engineers' Handbook (8th Edition), McGraw-Hill Professional Publishing, 2007. ProQuest Ebook Central,

.. [Wankat] Wankat, P C. Separation Process Engineering (3rd Edition), Prentice Hall, 2012, p. 33
