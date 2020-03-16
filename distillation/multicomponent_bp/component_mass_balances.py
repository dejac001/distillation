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


def example_6_1():
    """A distillation column with a partial reboiler and a total condenser
    is separating nC4 nC5 and nC8. The column has two equilibrium stages
    (a total of three equilibrium contacts), and the feed is a saturated liquid fed into
    the bottom stage of the column. The column operates at 2 atm. Feed rate is
    1000 kmol/h. zC4 = 0.20, zC5 = 0.35, zC8 = 0.45 (mole fraction).

    The reflux is a saturated liquid, and L/D=1.5.
    The distillate rate is D = 550 kmol/h.
    """
    # input data
    from distillation.equilibrium_data.depriester_charts import DePriester
    components = ['n-Butane', 'n-Pentane', 'n-Octane']
    z_feed = {
        'n-Butane': 0.2,
        'n-Pentane': 0.35,
        'n-Octane': 0.45
    }
    K = {
        key: DePriester(key) for key in components
    }
    p_feed = 101325. * 2  # Pressure in Pa
    Feed_flow = 1000.  # kmol/h
    Distillate_flow = 550.  # kmol/h
    reflux_ratio = 1.5  # L/D

    # bubble point calculation on the feed
    from distillation.bubble_point.calculation import bubble_point
    z_vals = [z_feed[key] for key in components]
    K_vals = [K[key].eval_SI for key in components]
    T_feed = bubble_point(z_vals, K_vals, p_feed, 400.)
    T_diff = T_feed - 273. - 60.
    print('Difference in bubble point temperature with Wankat is %e deg K' % T_diff)

    # make matrix assuming CMO
    N = 3  # number of equilibrium contacts
    L = np.zeros(N + 1)
    V = np.zeros(N + 1)
    L[0] = reflux_ratio * Distillate_flow
    L[1] = reflux_ratio * Distillate_flow
    L[2] = reflux_ratio * Distillate_flow + Feed_flow
    Bottoms_flow = Feed_flow - Distillate_flow
    L[3] = Bottoms_flow
    V[0] = 0.  # total condenser
    V[1:] = reflux_ratio * Distillate_flow + Distillate_flow

    # n-Butane-------------------
    F_stages = np.array([0., 0., Feed_flow, 0.])
    z = np.array([0., 0., z_feed['n-Butane'], 0.])
    K_stages = K['n-Butane'].eval_SI(T_feed, p_feed) * np.ones(N + 1)
    A, B, C, D = make_ABC(V, L, K_stages, F_stages, z, Distillate_flow, Bottoms_flow, N)

    # this matrix is for n-Butane
    A_wankat = np.array([-1, -1, -1])
    B_wankat = np.array([1.6, 6.00, 3.26, 10.17])
    C_wankat = np.array([-5.00, -2.26, -9.17])
    D_wankat = np.array([0, 0, 200, 0])

    print(np.linalg.norm(A - A_wankat))
    print(np.linalg.norm(B - B_wankat))
    print(np.linalg.norm(C - C_wankat))
    print(np.linalg.norm(D - D_wankat))

    l_C4 = solve_diagonal(A, B, C, D)
    l_C4_Wankat = np.array([280.217, 93.593, 124.492, 12.241])
    print(np.linalg.norm(l_C4 - l_C4_Wankat))

    # n-Pentane-------------------
    F_stages = np.array([0., 0., Feed_flow, 0.])
    z = np.array([0., 0., z_feed['n-Pentane'], 0.])
    K_stages = K['n-Pentane'].eval_SI(T_feed, p_feed) * np.ones(N + 1)
    A, B, C, D = make_ABC(V, L, K_stages, F_stages, z, Distillate_flow, Bottoms_flow, N)
    l_C5 = solve_diagonal(A, B, C, D)
    l_C5_Wankat = np.array([303.037, 289.141, 623.042, 148.284])
    print(np.linalg.norm(l_C5 - l_C5_Wankat))

    # n-Octane--------------------
    F_stages = np.array([0., 0., Feed_flow, 0.])
    z = np.array([0., 0., z_feed['n-Octane'], 0.])
    K_stages = K['n-Octane'].eval_SI(T_feed, p_feed) * np.ones(N + 1)
    A, B, C, D = make_ABC(V, L, K_stages, F_stages, z, Distillate_flow, Bottoms_flow, N)
    l_C8 = solve_diagonal(A, B, C, D)
    l_C8_Wankat = np.array([1.993, 27.741, 549.114, 450.273])
    print(np.linalg.norm(l_C8 - l_C8_Wankat))

    # calculate mole fractions
    x_C4 = l_C4 / (l_C4 + l_C5 + l_C8)
    x_C5 = l_C5 / (l_C4 + l_C5 + l_C8)
    x_C8 = l_C8 / (l_C4 + l_C5 + l_C8)

    # calculate temperature on each stage with a bubble point calculation
    # start with reboiler, N = 3
    T_wankat = [32.0, 44.8, 67.2, 103.4]
    T_stages = np.zeros(N+1)
    y_C4 = np.zeros(N+1)
    y_C5 = np.zeros(N+1)
    y_C8 = np.zeros(N+1)
    for i in range(N+1):
        x = [x_C4[i], x_C5[i], x_C8[i]]
        K_vals = [K['n-Butane'].eval_SI, K['n-Pentane'].eval_SI, K['n-Octane'].eval_SI]
        T = bubble_point(x, K_vals, p_feed, 300.)
        T_diff = T - 273. - T_wankat[i]
        print('Difference in T on stage %i with Wankat is %e deg K' % (i, T_diff))
        T_stages[i] = T

        y_C4[i] = K['n-Butane'].eval_SI(T, p_feed)*x_C4[i]
        y_C5[i] = K['n-Pentane'].eval_SI(T, p_feed)*x_C5[i]
        y_C8[i] = K['n-Octane'].eval_SI(T, p_feed)*x_C8[i]

    y_C4_Wankat = np.array([0.747, 0.496, 0.344, 0.134])
    y_C5_Wankat = np.array([0.243, 0.501, 0.620, 0.659])
    y_C8_Wankat = np.array([6.8e-5, 0.002, 0.036, 0.207])
    print(np.linalg.norm(y_C4-y_C4_Wankat))
    print(np.linalg.norm(y_C5-y_C5_Wankat))
    print(np.linalg.norm(y_C8-y_C8_Wankat))


if __name__ == '__main__':
    example_6_1()
