import numpy as np


def make_ABC(V: np.array, L: np.array, K: np.array, F: np.array, z: np.array,
             Distillate: float, Bottoms: float, N: int):
    """
    Distillation column with partial reboiler and total condenser

    .. note::
        K_j is assumed to depend on *T* and *p*, but not composition

    :param V: vapor molar flow rate out of stage 0 to *N*
    :param L: liquid molar flow rate out of stage 0 to *N*
    :param K: equilibrium expressions for stage 0 to *N*
    :param F: feed flow rate into stage for stage 0 to *N*
    :param z: feed composition into stage for stage 0 to *N*
    :param Distillate: distillate flow rate
    :param Bottoms: bottoms flow rate
    :param N: number of equilibrium stages
    :return: A, B, C, D
    """
    B = np.zeros(N + 1)  # diagonal
    A = -1 * np.ones(N)  # lower diagonal
    C = np.zeros(N)  # upper diagonal
    D = np.zeros(N + 1)

    assert abs(V[0]) < 1e-8, 'Vapor flow rate out of total condenser is non-zero!'
    # total condenser
    B[0] = 1. + Distillate / L[0]
    C[0] = -V[1] * K[1] / L[1]
    D[0] = F[0] * z[0]
    # reboiler
    B[N] = 1 + V[N] * K[N] / Bottoms
    D[N] = F[N] * z[N]

    D[1:N] = F[1:N] * z[1:N]
    B[1:N] = 1 + V[1:N] * K[1:N] / L[1:N]
    C[1:N] = -V[2:(N + 1)] * K[2:(N + 1)] / L[2:(N + 1)]
    return A, B, C, D


def solve_diagonal(lower, diag, upper, b):
    """Solve matrix Ax=b when A is diagonal"""
    from scipy.sparse import linalg, diags
    N = len(diag)
    A = diags(
        diagonals=[lower, diag, upper],
        offsets=[-1, 0, 1],
        shape=(N, N),
        format='csr'
    )
    return linalg.spsolve(A, b)


def solve_component_mass_balances(*args):
    """
    Distillation column with partial reboiler and total condenser

    .. note::
        K_j is assumed to depend on *T* and *p*, but not composition

    :param V: vapor molar flow rate out of stage 0 to *N*
    :param L: liquid molar flow rate out of stage 0 to *N*
    :param K: equilibrium expressions for stage 0 to *N*
    :param F: feed flow rate into stage for stage 0 to *N*
    :param z: feed composition into stage for stage 0 to *N*
    :param Distillate: distillate flow rate
    :param Bottoms: bottoms flow rate
    :param N: number of equilibrium stages
    :return: l
    """
    A, B, C, D = make_ABC(*args)
    return solve_diagonal(A, B, C, D)


def main(components: list, F: float, P: float, z_feed: list, RR: float, D: float, N: int, feed_stage: int):
    """Distillation column with partial reboiler and total condenser.
    Feed is saturated liquid.

    :param components: list of component names
    :param F: feed molar flow rate
    :param P: pressure (constant throughout column), [Pa]
    :param z_feed: mole fractions of each component, ordered
    :param RR: reflux ratio (L/D)
    :param D: distillate molar flow rate
    :param N: number of equilibrium contacts
    :param feed_stage: stage where feed is input
    :return: x_ij, y_ij, T_j for each stage j
    """
    # input data
    from distillation.equilibrium_data.depriester_charts import DePriester
    import os
    from distillation.bubble_point.calculation import bubble_point
    file = os.path.join('..', '..', 'distillation', 'equilibrium_data', 'depriester.csv')
    import numpy as np
    K = {
        key: DePriester(key, file) for key in components
    }

    # bubble point calculation on the feed
    K_vals = [K[key].eval_SI for key in components]
    T_feed = bubble_point(z_feed, K_vals, P, 300.)

    # make matrix assuming CMO
    L = np.zeros(N + 1)
    V = np.zeros(N + 1)
    L[:feed_stage] = RR * D
    L[feed_stage:N] = RR * D + F
    Bottoms_flow = F - D
    L[N] = Bottoms_flow
    V[0] = 0.  # total condenser
    V[1:] = RR * D + D

    F_stages = np.zeros(N + 1)
    F_stages[feed_stage] = F

    l = {}
    for i, component in enumerate(components):
        z = np.zeros(N + 1)
        z[feed_stage] = z_feed[i]
        K_stages = K[component].eval_SI(T_feed, P) * np.ones(N + 1)
        l[component] = solve_component_mass_balances(
            V, L, K_stages, F_stages, z, D, Bottoms_flow, N
        )

    x = {}
    for component in components:
        x[component] = np.zeros(N + 1)
        for i in range(N + 1):
            x[component][i] = l[component][i]/sum(l[c][i] for c in components)

    # calculate temperature on each stage with a bubble point calculation
    T_stages = np.zeros(N + 1)
    for i in range(N + 1):
        x_vals = [x[c][i] for c in components]
        K_vals = [K[i].eval_SI for i in components]
        T = bubble_point(x_vals, K_vals, P, T_feed)
        T_stages[i] = T

    # calculate mole fractions
    y = {
        key: np.zeros(N + 1) for key in components
    }
    for i in range(N + 1):
        for component in components:
            y[component][i] = K[component].eval_SI(T_stages[i], P)*x[component][i]

    return x, y, T_stages
