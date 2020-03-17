import numpy as np
from distillation.solvers import solve_diagonal


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


def one_iteration(*args):
    pass

