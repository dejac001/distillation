def read_chart(f_name):
    data = {}
    with open(f_name, 'r') as f:
        header = next(f).split(',')
        for line in f:
            compound, *vals = line.split(',')
            data[compound] = {
                key: val for key, val in zip(header[1:], vals)
            }

    return data


class DePriester:
    def __init__(self, compound):
        from distillation import ROOT_DIR, os
        f_name = os.path.join(ROOT_DIR, 'equilibrium_data', 'depriester.csv')
        data = read_chart(f_name)
        assert compound in data.keys(), 'Compound not found!'
        my_data = data[compound]

        self.a_T1 = float(my_data.pop('a_T1'))
        self.a_T2 = float(my_data.pop('a_T2'))
        self.a_T6 = float(my_data.pop('a_T6'))
        self.a_p1 = float(my_data.pop('a_p1'))
        self.a_p2 = float(my_data.pop('a_p2'))
        self.a_p3 = float(my_data.pop('a_p3'))

    def eval(self, T, p):
        """

        :param T: temperature in Rankine
        :param p: pressure in psia
        :return: K-value for component at specific *T* and *p*
        """
        import numpy as np
        return np.exp(
            self.a_T1 / T / T + self.a_T2 / T + self.a_T6
            + self.a_p1 * np.log(p) + self.a_p2 / p / p + self.a_p3 / p
        )

    def eval_SI(self, T, p):
        """

        :param T: temperature in K
        :param p: pressure in Pa
        :return: K-value for component at specific *T* and *p*
        """
        from distillation.unit_conversions.temperature import Kelvin_to_Rankine
        from distillation.unit_conversions.pressure import Paa_to_psia
        return self.eval(
            Kelvin_to_Rankine(T), Paa_to_psia(p)
        )
