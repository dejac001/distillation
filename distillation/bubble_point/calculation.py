def residual(T, x, K, p):
    return sum(x_i*K_i(T, p) for x_i, K_i in zip(x, K)) - 1.


def bubble_point(x, K, p, T_guess):
    """
    Find the temperature that satisfies the stoichiometric equation

    .. math::
        \sum_i y_i = \sum_i K_i (T) x_i = 1.0

    :param x: mole fractions in liquid
    :param p: total pressure
    :param K: functions calculating K for each component
    :param T_guess: guess temperature
    :return: Temperature at which the liquid mixture begins to boil.
    """
    from scipy.optimize import root_scalar
    sol = root_scalar(residual, args=(x, K, p), x0=T_guess-10, x1=T_guess+10)  #, method='newton')
    return sol.root
