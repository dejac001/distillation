import unittest


class MyTestCase(unittest.TestCase):
    def test_ABC_construction(self):
        """A distillation column with a partial reboiler and a total condenser
        is separating nC4 nC5 and nC8. The column has two equilibrium stages
        (a total of three equilibrium contacts), and the feed is a saturated liquid fed into
        the bottom stage of the column. The column operates at 2 atm. Feed rate is
        1000 kmol/h. zC4 = 0.20, zC5 = 0.35, zC8 = 0.45 (mole fraction).

        The reflux is a saturated liquid, and L/D=1.5.
        The distillate rate is D = 550 kmol/h.
        """

        # bubble point calculation on the feed


        import numpy as np
        A_wankat = np.array([-1, -1, -1])
        B_wankat = np.array([1.6, 6.00, 3.26, 10.17])
        C_wankat = np.array([-5.00, -2.26, -9.17])
        D_wankat = np.array([0, 0, 200, 0])


if __name__ == '__main__':
    unittest.main()
