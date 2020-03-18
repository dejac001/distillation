from distillation.amundson_1958.component_mass_balances import make_ABC
from distillation.bubble_point.calculation import bubble_point
from distillation.solvers import solve_diagonal
import os
from scipy.sparse import linalg, diags
import numpy as np


class Model:

    def __init__(self, components: list=None, F: float=0., P: float=101325.,
                 z_feed: list=None, RR: float=1, D: float=0, N: int=1, feed_stage: int=0):
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
        file = os.path.join('..', '..', 'distillation', 'equilibrium_data', 'depriester.csv')
        self.K_func = {
            key: DePriester(key, file) for key in self.components
        }
        file = os.path.join('..', '..', 'distillation', 'equilibrium_data', 'heat_capacity_liquid.csv')
        self.CpL_func = {
            key: CpL(key, file) for key in self.components
        }
        self.CpV_func = {
            key: CpV() for key in self.components
        }
        file = os.path.join('..', '..', 'distillation', 'equilibrium_data', 'heats_of_vaporization.csv')
        self.dH_func = {
            key: dH_vap(key, file) for key in self.components
        }
        self.T_ref = {
            key: val.T_ref for key, val in self.dH_func.items()
        }

        # create matrices for variables
        self.L = np.zeros(self.N + 1)
        self.V = np.zeros(self.N + 1)
        self.L[N] = self.B
        self.V[0] = 0.  # total condenser
        self.F = np.zeros(N + 1)
        self.F[self.feed_stage] = self.F_feed
        self.z = {
            key: np.zeros(N + 1) for key in components
        }
        self.x = {}
        self.y = {}
        for component in self.components:
            self.z[component][feed_stage] = self.z_feed[component]
            self.x[component] = self.z_feed[component]*np.ones(self.N + 1)
            self.y[component] = self.z_feed[component]*np.ones(self.N + 1)

        # perform bubble point calculation on the feed
        self.T_feed = self.calculate_T_feed()

        # initialize T and K for all stages at feed temperature
        self.T_old = self.T_feed * np.ones(self.N + 1)
        self.T_new = self.T_feed * np.ones(self.N + 1)
        self.K = {
            key: self.K_func[key].eval_SI(self.T_feed, self.P_feed)*np.ones(self.N + 1) for key in self.components
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
            h_j = \sum_i x_{ij}\overbar{h}_i(T_j)

        where the overbar indicates the pure component enthalpy
        units of J/kmol
        """
        return sum(
            self.x[c][stage]*self.h_pure_rule(c, self.T_new[stage]) for c in self.components
        )

    def h_feed_rule(self, stage):
        """Enthalpy of liquid in feed mixture
        Calculated for ideal mixture

        .. math::
            h_j = \sum_i x_{ij}\overbar{h}_i(T_j)

        where the overbar indicates the pure component enthalpy
        units of J/kmol
        """
        return sum(
           self.z[c][stage]*self.h_pure_rule(c, self.T_feed) for c in self.components
        )

    def H_pure_rule(self, c, T):
        """Rule for vapor enthalpy of pure component"""
        return self.CpV_func[c].integral_dT(self.T_ref[c], T) + self.dH_func[c].eval()

    def H_j_rule(self, stage):
        """Enthalpy of vapor on stage *j*.
        Calculated for ideal mixture

        .. math::
            H_j = \sum_i y_{ij}\overbar{H}_i(T_j)

        where the overbar indicates the pure component enthalpy
        units of J/kmol
        """
        return sum(
            self.y[c][stage] * self.H_pure_rule(c, self.T_new[stage]) for c in self.components
        )

    def Q_condenser_rule(self):
        """Condenser requirement can be determined from balances around total condenser"""
        return self.D * (1. + self.RR) * (self.h_j_rule(0) - self.H_j_rule(1))

    def Q_reboiler_rule(self):
        """Condenser requirement can be determined from balances around total condenser"""
        return self.D*self.h_j_rule(0) + self.B*self.h_j_rule(self.N) \
               - self.F_feed*self.h_feed_rule(self.feed_stage) - self.Q_condenser_rule()

    def run(self, num_iter=1):
        self.initialize_CMO()
        for i in range(1, num_iter + 1):
            self.solve_component_mass_balances()
            if self.temperature_is_converged():
                print('temperature converged in %i iterations' % i)
                for j in range(1, num_iter + 1):
                    converged = self.solve_energy_balances()
                    if converged:
                        print('flow rates converged in %i iterations' % i)

                print('flow ratesdid not converge in %i iterations' % num_iter)
                return

        print('temperature did not converge in %i iterations' % num_iter)

    def temperature_is_converged(self):
        eps = np.abs(self.T_new - self.T_old)
        if eps.max() < self.temperature_tol:
            return True

        self.T_old = self.T_new[:]
        return False

    def calculate_T_feed(self):
        """Feed mixture is saturated liquid. Calculate T_feed with bubble point calculation"""
        return bubble_point(
            [self.z_feed[i] for i in self.components],
            [self.K_func[i].eval_SI for i in self.components], self.P_feed, 300.
        )

    def initialize_CMO(self):
        """Make initial guesses for flow rates and temperatures"""

        # make matrix assuming CMO
        self.L[:self.feed_stage] = self.RR * self.D
        self.L[self.feed_stage:self.N] = self.RR * self.D + self.F_feed
        self.V[1:] = self.RR * self.D + self.D

    def solve_component_mass_balances(self):
        """Calculate component liquid flow rates."""
        l = {}
        for i, component in enumerate(self.components):
            # todo: dont need to calculate D and A as these are constant
            A, B, C, D = make_ABC(
                self.V, self.L, self.K[component], self.F, self.z[component], self.D, self.B, self.N
            )
            l[component] = solve_diagonal(A, B, C, D)

        # todo: vectorize all this with matrix multiplication
        # update from old calculations
        for i in range(self.N + 1):

            # update mole fractions
            for c in self.components:
                self.x[c][i] = l[c][i]/sum(l[c][i] for c in self.components)

            # calculate stage temperature now that all liquid-phase mole fractions are known
            K_vals = [self.K_func[c].eval_SI for c in self.components]
            x_vals = [self.x[c][i] for c in self.components]
            self.T_new[i] = self.T_old[i] + self.df * (
                    bubble_point(x_vals, K_vals, self.P_feed, self.T_old[i]) - self.T_old[i]
            )
            for c in self.components:
                # update K
                self.K[c][i] = self.K_func[c].eval_SI(self.T_old[i], self.P_feed)
                # update y
                self.y[c][i] = self.K[c][i]*self.x[c][i]

        # update L
        # for i in range(self.N + 1):
        #     self.L[i] = sum(l[c][i] for c in self.components)
        # calculate V from overall mass balance at each equilibrium stage
        # self.V[self.N] = self.L[self.N-1] - self.B
        # for i in range(self.N-1, 1):
        #     self.V[i] = self.L[i-1] + self.V[i + 1] + self.F[i] - self.L[i]
        # self.V[1] = self.B + self.L[0]
        # self.V[0] = 0. # total condenser

    def solve_energy_balances(self):
        """Solve energy balances"""
        V_old = self.V[:]
        L_old = self.L[:]

        BE = np.zeros(self.N + 1)
        CE = np.zeros(self.N)
        DE = np.zeros(self.N + 1)

        # total condenser
        BE[0] = 0.
        CE[0] = self.h_j_rule(0) - self.H_j_rule(1)
        DE[0] = self.F[0]*self.h_feed_rule(0) + self.Q_condenser_rule()

        # partial reboiler
        BE[self.N] = self.H_j_rule(self.N) - self.H_j_rule(self.N-1)
        DE[self.N] = self.F[self.N]*self.h_feed_rule(self.N) + self.Q_reboiler_rule() \
                     + self.B*(self.h_j_rule(self.N-1)-self.h_j_rule(self.N)) \
                     - self.F[self.N-1]*self.h_j_rule(self.N-1)

        # stages 1 to N-1
        BE[1:self.N] = list(self.H_j_rule(j) - self.h_j_rule(j-1) for j in range(1, self.N))
        CE[1:self.N] = list(self.h_j_rule(j) - self.H_j_rule(j+1) for j in range(1, self.N))
        DE[1:self.N] = list(self.F[j]*self.h_feed_rule(j) + self.D*(self.h_j_rule(j-1) - self.h_j_rule(j))
                        - sum(self.F[k] for k in range(j+1))*self.h_j_rule(j)
                        + sum(self.F[k] for k in range(j))*self.h_j_rule(j-1)
                        for j in range(1, self.N))
        A = diags(
            diagonals=[BE[1:], CE[1:]],
            offsets=[0, 1],
            shape=(self.N, self.N),
            format='csr'
        )
        self.V[1:] = linalg.spsolve(A, DE[1:])
        for i in range(self.N):
            self.L[i] = self.V[i+1] - self.D + sum(self.F[k] for k in range(i + 1))
        self.L[self.N] = self.B
        return self.L_is_converged(L_old) and self.V_is_converged(V_old)

    def L_is_converged(self, L_old):
        eps_L = self.relative_error(self.L, L_old)
        if eps_L.max() < self.flow_rate_tol:
            return True

        return False

    def relative_error(self, new, old):
        return np.abs((new - old)/new)

    def V_is_converged(self, V_old):
        eps_V = self.relative_error(self.V[1:], V_old[1:])
        if eps_V.max() < self.flow_rate_tol:
            return True

        return False


if __name__ == '__main__':
    cls = Model(
        ['n-Butane', 'n-Pentane', 'n-Octane'],
        1000., 101325. * 2,
        [0.2, 0.35, 0.45],
        1.5, 550., 3, 2
    )
    cls.run(num_iter=100)
