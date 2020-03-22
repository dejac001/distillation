from distillation.bubble_point.calculation import bubble_point
from distillation.solvers import solve_diagonal
from scipy.sparse import linalg, diags
import numpy as np


class Model:
    def __init__(self, components: list = None, F: float = 0., P: float = 101325.,
                 z_feed: list = None, RR: float = 1, D: float = 0, N: int = 1, feed_stage: int = 0,
                 T_feed_guess: float = 300.):
        """Distillation column with partial reboiler and total condenser.
        Feed is saturated liquid.

        .. todo::
            implement component liquid flow rates (:attr:`Model.l`) throughout instead of mole fractions

        :param components: list of component names
        :param F: feed molar flow rate
        :param P: pressure (constant throughout column), [Pa]
        :param z_feed: mole fractions of each component, ordered
        :param RR: reflux ratio (L/D)
        :param D: distillate molar flow rate
        :param N: number of equilibrium contacts
        :param feed_stage: stage where feed is input
        :param T_feed_guess: guess temperature for feed stage
        """
        self.flow_rate_tol = 1.e-4
        self.temperature_tol = 1.e-2
        from distillation.equilibrium_data.depriester_charts import DePriester
        from distillation.equilibrium_data.heat_capacity_liquid import CpL
        from distillation.equilibrium_data.heat_capacity_vapor import CpV
        from distillation.equilibrium_data.heats_of_vaporization import dH_vap
        self.components = components
        self.F_feed = F
        self.P_feed = P
        self.z_feed = {key: val for key, val in zip(components, z_feed)}
        self.RR = RR
        self.D = D
        self.B = F - D
        self.N = N
        self.feed_stage = feed_stage
        self.T_feed_guess = T_feed_guess
        self.K_func = {
            key: DePriester(key) for key in self.components
        }
        self.CpL_func = {
            key: CpL(key) for key in self.components
        }
        self.CpV_func = {
            key: CpV() for key in self.components
        }
        self.dH_func = {
            key: dH_vap(key) for key in self.components
        }
        self.T_ref = {
            key: val.T_ref for key, val in self.dH_func.items()
        }

        # create matrices for variables
        self.num_stages = self.N + 1
        self.stages = range(self.num_stages)
        self.L = np.zeros(self.num_stages)
        self.V = np.zeros(self.num_stages)
        self.L_old = np.zeros(self.num_stages)
        self.V_old = np.zeros(self.num_stages)
        self.L[N] = self.B
        self.V[0] = 0.  # total condenser
        self.F = np.zeros(self.num_stages)
        self.F[self.feed_stage] = self.F_feed
        self.z = {
            key: np.zeros(self.num_stages) for key in components
        }
        # todo: implement component flow rates throughout
        self.l = {
            key: np.zeros(self.num_stages) for key in components
        }
        for component in self.components:
            self.z[component][feed_stage] = self.z_feed[component]

        self.T_feed = self.T_feed_guess
        self.T = self.T_feed_guess * np.ones(self.num_stages)
        self.T_old = np.zeros(self.num_stages)
        self.K = {
            key: self.K_func[key].eval_SI(self.T_feed, self.P_feed) * np.ones(self.num_stages) for key in
            self.components
        }

        # solver parameters
        self.df = 1.  # Dampening factor to prevent excessive oscillation of temperatures

    def h_pure_rule(self, c, T):
        """rule for liquid enthalpy of pure component"""
        return self.CpL_func[c].integral_dT(self.T_ref[c], T)

    def h_j_rule(self, stage):
        """Enthalpy of liquid on stage *j*.
        Calculated for ideal mixture

        .. math::

            h_j = \\sum_i x_{ij}h^*_i(T_j)

        where the asterisk indicates the pure component enthalpy

        :return: :math:`h_j` [J/kmol]
        """
        return sum(
            self.x_ij_expr(c, stage) * self.h_pure_rule(c, self.T[stage]) for c in self.components
        )

    def x_ij_expr(self, i, j):
        """

        :param i: component name
        :param j: stage number
        :return: mole fraction on stage
        """
        return self.l[i][j] / self.L[j]

    def h_feed_rule(self, stage):
        """Enthalpy of liquid in feed mixture
        Calculated for ideal mixture

        .. math::

            h = \\sum_i x_{ij}h^*_i(T_j)

        where the asterisk indicates the pure component enthalpy

        :return: :math:`h` [J/kmol]
        """
        return sum(
            self.z[c][stage] * self.h_pure_rule(c, self.T_feed) for c in self.components
        )

    def H_pure_rule(self, c, T):
        """Rule for vapor enthalpy of pure component"""
        return self.CpV_func[c].integral_dT(self.T_ref[c], T) + self.dH_func[c].eval()

    def H_j_rule(self, stage):
        """Enthalpy of vapor on stage *j*.
        Calculated for ideal mixture

        .. math::
            H_j = \\sum_i y_{ij}H^*_i(T_j)

        where the asterisk indicates the pure component enthalpy

        .. todo::
            convert y mole fractions to dynamic expression

        :return: :math:`H_j` [J/kmol]
        """
        return sum(
            self.y_ij_expr(c, stage) * self.H_pure_rule(c, self.T[stage]) for c in self.components
        )

    def y_ij_expr(self, i, j):
        """

        :param i: component name
        :param j: stage number
        :return: gas-phase mole fraction on stage
        """
        return self.K_func[i].eval_SI(self.T[j], self.P_feed) * self.x_ij_expr(i, j)

    def Q_condenser_rule(self):
        """Condenser requirement can be determined from balances around total condenser"""
        return self.D * (1. + self.RR) * (self.h_j_rule(0) - self.H_j_rule(1))

    def Q_reboiler_rule(self):
        return self.D * self.h_j_rule(0) + self.B * self.h_j_rule(self.N) \
               - self.F_feed * self.h_feed_rule(self.feed_stage) - self.Q_condenser_rule()

    def step_3_to_step_6(self):
        num_iter = 0
        while not self.T_is_converged():
            self.update_K_values()
            for i in self.components:
                self.solve_component_mass_bal(i)
            self.update_T_values()
            num_iter += 1
        print('while loop exits with %i iterations' % num_iter)

    def run(self):
        self.generate_initial_guess()
        self.step_3_to_step_6()
        self.solve_energy_balances()
        main_loop = 0
        while not self.flow_rates_converged():
            for i in self.components:
                self.solve_component_mass_bal(i)
            self.update_T_values()
            self.step_3_to_step_6()
            self.solve_energy_balances()
            main_loop += 1
            print(main_loop)

    def update_K_values(self):
        """
        .. include:: step3.rst

        """
        for c in self.components:
            self.K[c][:] = self.K_func[c].eval_SI(self.T[:], self.P_feed)

        self.T_old[:] = self.T[:]

    def update_T_values(self):
        """Update temperatures in all stages
        by performing bubble point calculation

        .. todo::
            vectorize with matrix multiplication

        """
        # update from old calculations
        for i in range(self.num_stages):
            # calculate stage temperature now that all liquid-phase mole fractions are known
            K_vals = [self.K_func[c].eval_SI for c in self.components]
            l_total = sum(self.l[c][i] for c in self.components)
            x_vals = [self.l[c][i] / l_total for c in self.components]
            self.T[i] = self.T_old[i] + self.df * (
                    bubble_point(x_vals, K_vals, self.P_feed, self.T_old[i]) - self.T_old[i]
            )

    def generate_initial_guess(self):
        """
        .. include:: step2.rst

        """
        # initialize temperatures as T (feed)
        self.T_feed = bubble_point(
            [self.z_feed[i] for i in self.components],
            [self.K_func[i].eval_SI for i in self.components], self.P_feed, self.T_feed_guess
        )
        self.T[:] = self.T_feed

        # initialize L, V with CMO
        self.L[:self.feed_stage] = self.RR * self.D
        self.L[self.feed_stage:self.N] = self.RR * self.D + self.F_feed
        self.V[1:] = self.RR * self.D + self.D

    def T_is_converged(self):
        """
        .. include:: temp-converge.rst


        :return: True if T is converged, else False
        """
        eps = np.abs(self.T - self.T_old)
        return eps.max() < self.temperature_tol

    def solve_component_mass_bal(self, component):
        """Solve component mass balances

        .. todo:
            dont need to calculate D and A here as they are constant

        """
        A, B, C, D = make_ABC(
            self.V, self.L, self.K[component], self.F, self.z[component], self.D, self.B, self.N
        )
        self.l[component][:] = solve_diagonal(A, B, C, D)

    def update_flow_rates(self):
        for i in self.stages:
            self.L[i] = sum(self.l[c][i] for c in self.components)

        self.V[0] = 0.  # total condenser
        self.V[1] = (self.RR + 1.) * self.D
        for i in range(2, self.num_stages):
            self.V[i] = self.L[i - 1] + self.D - sum(self.F[k] for k in range(i))

    def solve_energy_balances(self):
        """Solve energy balances"""

        self.L_old[:] = self.L[:]
        self.V_old[:] = self.V[:]

        BE = np.zeros(self.num_stages)
        CE = np.zeros(self.num_stages)
        DE = np.zeros(self.num_stages)

        # total condenser
        BE[0] = 0.
        CE[0] = self.h_j_rule(0) - self.H_j_rule(1)
        DE[0] = self.F[0] * self.h_feed_rule(0) + self.Q_condenser_rule()

        # stages 1 to N-1
        for j in range(1, self.N):
            BE[j] = self.H_j_rule(j) - self.h_j_rule(j - 1)
            CE[j] = self.h_j_rule(j) - self.H_j_rule(j + 1)
            DE[j] = self.F[j] * self.h_feed_rule(j) - self.D * (self.h_j_rule(j - 1) - self.h_j_rule(j)) \
                    - sum(self.F[k] for k in range(j + 1)) * self.h_j_rule(j) \
                    + sum(self.F[k] for k in range(j)) * self.h_j_rule(j - 1)

        # partial reboiler
        BE[self.N] = self.H_j_rule(self.N) - self.h_j_rule(self.N - 1)
        DE[self.N] = self.F[self.N] * self.h_feed_rule(self.N) + self.Q_reboiler_rule() \
                     - self.B * (self.h_j_rule(self.N - 1) - self.h_j_rule(self.N)) \
                     - self.F[self.N - 1] * self.h_j_rule(self.N - 1)

        A = diags(
            diagonals=[BE[1:], CE[1:-1]],
            offsets=[0, 1],
            shape=(self.N, self.N),
            format='csr'
        )
        self.V[1:] = linalg.spsolve(A, DE[1:])
        self.L[0] = self.RR * self.D
        for i in range(1, self.N):
            self.L[i] = self.V[i + 1] - self.D + sum(self.F[k] for k in range(i + 1))
        self.L[self.N] = self.B

    def flow_rates_converged(self):
        """Determine if flow rates are converged

        Use the mathematical criterion in :meth:`Model.is_below_relative_error`

        """
        return self.is_below_relative_error(self.L, self.L_old) and self.is_below_relative_error(self.V[1:],
                                                                                                 self.V_old[1:])

    def is_below_relative_error(self, new, old):
        """Determine relative error between two vectors

        .. math::

            \\sqrt{\\left(\\frac{X_{\\mathrm{new}} - X_{\\mathrm{old}}}{X_{\\mathrm{new}}}\\right)^2} < \\epsilon

        The flow rate tolerance, :math:`\epsilon`,
        is found in the attribute :attr:`distillation.amundson_1958.main.Model.flow_rate_tol`

        :param new:
        :param old:
        :rtype: bool

        """
        return np.abs((new - old) / new).max() < self.flow_rate_tol


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


if __name__ == '__main__':
    cls = Model(
        ['n-Butane', 'n-Pentane', 'n-Octane'],
        1000., 101325. * 2,
        [0.2, 0.35, 0.45],
        1.5, 550., 3, 2
    )
    cls.run(num_iter=100)
