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
        from distillation.amundson_1958.main import Model
        cls = Model(
            ['n-Butane', 'n-Pentane', 'n-Octane'],
            1000., 101325. * 2,
            [0.2, 0.35, 0.45],
            1.5, 550., 3, 2
        )
        cls.run()
        import numpy as np
        x_Wankat = {
            'n-Butane': np.array([0.497, 0.228, 0.096, 0.020]),
            'n-Pentane': np.array([0.508, 0.704, 0.481, 0.243]),
            'n-Octane': np.array([0.0034, 0.068, 0.423, 0.737]),
        }
        y_Wankat = {
            'n-Butane': np.array([0.747, 0.496, 0.344, 0.134]),
            'n-Pentane': np.array([0.243, 0.501, 0.620, 0.659]),
            'n-Octane': np.array([6.8e-5, 0.002, 0.036, 0.207]),
        }
        T_Wankat = np.array([32. + 273., 44.8 + 273., 67.2 + 273., 103.4 + 273.])

        for key in cls.x.keys():
            self.assertLess(np.linalg.norm(cls.y[key]-y_Wankat[key]), 1.5)
            self.assertLess(np.linalg.norm(cls.x[key]-x_Wankat[key]), 1.52)

        self.assertLess(np.linalg.norm(T_Wankat - cls.T_old), 1.)


if __name__ == '__main__':
    unittest.main()
