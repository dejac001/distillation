from distillation.multicomponent_bp.component_mass_balances import make_ABC
from distillation.bubble_point.calculation import bubble_point
from distillation.solvers import solve_diagonal
import numpy as np


class Model:
    def __init__(self, components: list, F: float, P: float,
                 z_feed: list, RR: float, D: float, N: int, feed_stage: int):
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
        from distillation.equilibrium_data.depriester_charts import DePriester
        import os
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

    def run(self, num_iter=1):
        self.initialize_CMO()
        for i in range(1, num_iter + 1):
            self.solve()
            if self.temperature_is_converged():
                print('temperature converged in %i iterations' % i)
                return

        print('temperature did not converge in %i iterations' % num_iter)

    def temperature_is_converged(self):
        tol = 1e-2
        eps = np.abs(self.T_new - self.T_old)
        if eps.max() < tol:
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

    def solve(self):
        """Calculate component liquid flow rates."""
        l = {}
        for i, component in enumerate(self.components):
            A, B, C, D = make_ABC(
                self.V, self.L, self.K[component], self.F, self.z[component], self.D, self.B, self.N
            )
            l[component] = solve_diagonal(A, B, C, D)

        # update from old calculations
        for i in range(self.N + 1):
            # update L
            # self.L[i] = sum(l[c][i] for c in self.components)

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

        # calculate V from overall mass balance at each equilibrium stage
        # self.V[self.N] = self.L[self.N-1] - self.B
        # for i in range(self.N-1, 1):
        #     self.V[i] = self.L[i-1] + self.V[i + 1] + self.F[i] - self.L[i]
        # self.V[1] = self.B + self.L[0]
        # self.V[0] = 0. # total condenser


if __name__ == '__main__':
    cls = Model(
        ['n-Butane', 'n-Pentane', 'n-Octane'],
        1000., 101325. * 2,
        [0.2, 0.35, 0.45],
        1.5, 550., 3, 2
    )
    cls.run(num_iter=10)
