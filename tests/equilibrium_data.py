import unittest
import os


class MyTestCase(unittest.TestCase):
    def test_heat_capacities(self):
        from distillation.equilibrium_data.heat_capacity_liquid import CpL
        components = ['n-Butane', 'n-Pentane', 'n-Octane']
        file = os.path.join('..', 'distillation', 'equilibrium_data', 'heat_capacity_liquid.csv')
        for key in components:
            cls = CpL(key, file)
            self.assertAlmostEqual(1., cls.eval(cls.T_min)/cls.value_T_min, places=4)


if __name__ == '__main__':
    unittest.main()
