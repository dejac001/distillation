import unittest


class MyTestCase(unittest.TestCase):
    def test_example_6_1_Wankat(self):
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
        import os
        from distillation.bubble_point.calculation import bubble_point
        from distillation.multicomponent_bp.component_mass_balances import make_ABC, solve_diagonal
        import numpy as np
        components = ['n-Butane', 'n-Pentane', 'n-Octane']
        z_feed = {
            'n-Butane': 0.2,
            'n-Pentane': 0.35,
            'n-Octane': 0.45
        }
        file = os.path.join('..', '..', 'distillation', 'equilibrium_data', 'depriester.csv')
        K = {
            key: DePriester(key, file) for key in components
        }
        p_feed = 101325. * 2  # Pressure in Pa
        Feed_flow = 1000.  # kmol/h
        Distillate_flow = 550.  # kmol/h
        reflux_ratio = 1.5  # L/D

        # bubble point calculation on the feed
        z_vals = [z_feed[key] for key in components]
        K_vals = [K[key].eval_SI for key in components]
        T_feed = bubble_point(z_vals, K_vals, p_feed, 400.)
        self.assertAlmostEqual(T_feed/100, (60.+273)/100, places=1)

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

        self.assertLess(np.linalg.norm(A - A_wankat), 1e-8)
        self.assertLess(np.linalg.norm(B - B_wankat), 0.1)
        self.assertLess(np.linalg.norm(C - C_wankat), 0.1)
        self.assertLess(np.linalg.norm(D - D_wankat), 1e-8)

        l_C4 = solve_diagonal(A, B, C, D)
        l_C4_Wankat = np.array([280.217, 93.593, 124.492, 12.241])
        self.assertLess(np.linalg.norm(l_C4 - l_C4_Wankat), 2.)

        # n-Pentane-------------------
        F_stages = np.array([0., 0., Feed_flow, 0.])
        z = np.array([0., 0., z_feed['n-Pentane'], 0.])
        K_stages = K['n-Pentane'].eval_SI(T_feed, p_feed) * np.ones(N + 1)
        A, B, C, D = make_ABC(V, L, K_stages, F_stages, z, Distillate_flow, Bottoms_flow, N)
        l_C5 = solve_diagonal(A, B, C, D)
        l_C5_Wankat = np.array([303.037, 289.141, 623.042, 148.284])
        self.assertLess(np.linalg.norm(l_C5 - l_C5_Wankat), 3.)

        # n-Octane--------------------
        F_stages = np.array([0., 0., Feed_flow, 0.])
        z = np.array([0., 0., z_feed['n-Octane'], 0.])
        K_stages = K['n-Octane'].eval_SI(T_feed, p_feed) * np.ones(N + 1)
        A, B, C, D = make_ABC(V, L, K_stages, F_stages, z, Distillate_flow, Bottoms_flow, N)
        l_C8 = solve_diagonal(A, B, C, D)
        l_C8_Wankat = np.array([1.993, 27.741, 549.114, 450.273])
        self.assertLess(np.linalg.norm(l_C8 - l_C8_Wankat), 16.)

        # calculate mole fractions
        x_C4 = l_C4 / (l_C4 + l_C5 + l_C8)
        x_C5 = l_C5 / (l_C4 + l_C5 + l_C8)
        x_C8 = l_C8 / (l_C4 + l_C5 + l_C8)

        # calculate temperature on each stage with a bubble point calculation
        # start with reboiler, N = 3
        T_wankat = [32.0, 44.8, 67.2, 103.4]
        T_stages = np.zeros(N + 1)
        y_C4 = np.zeros(N + 1)
        y_C5 = np.zeros(N + 1)
        y_C8 = np.zeros(N + 1)
        for i in range(N + 1):
            x = [x_C4[i], x_C5[i], x_C8[i]]
            K_vals = [K['n-Butane'].eval_SI, K['n-Pentane'].eval_SI, K['n-Octane'].eval_SI]
            T = bubble_point(x, K_vals, p_feed, 300.)
            self.assertAlmostEqual(T / 100, (T_wankat[i] + 273) / 100, places=1)
            T_stages[i] = T

            y_C4[i] = K['n-Butane'].eval_SI(T, p_feed) * x_C4[i]
            y_C5[i] = K['n-Pentane'].eval_SI(T, p_feed) * x_C5[i]
            y_C8[i] = K['n-Octane'].eval_SI(T, p_feed) * x_C8[i]

        y_C4_Wankat = np.array([0.747, 0.496, 0.344, 0.134])
        y_C5_Wankat = np.array([0.243, 0.501, 0.620, 0.659])
        y_C8_Wankat = np.array([6.8e-5, 0.002, 0.036, 0.207])
        self.assertLess(np.linalg.norm(y_C4 - y_C4_Wankat), 0.01)
        self.assertLess(np.linalg.norm(y_C5 - y_C5_Wankat), 0.01)
        self.assertLess(np.linalg.norm(y_C8 - y_C8_Wankat), 0.01)


if __name__ == '__main__':
    unittest.main()
