class CpL:
    """Heat capacities of liquids"""
    def __init__(self, compound_name):
        from distillation import os, ROOT_DIR
        file = os.path.join(ROOT_DIR, 'equilibrium_data', 'heat_capacity_liquid.csv')
        self.num_constants = 5
        self.units = 'J/kmol/K'
        with open(file, 'r') as f:
            header = next(f).rstrip('\n').split(',')
            for line in f:
                vals = line.rstrip('\n').split(',')
                if vals[0] == compound_name:
                    self.constants = {
                        key: float(vals[header.index('C%i' % key)]) for key in range(1, self.num_constants + 1)
                    }
                    self.value_T_min = float(vals[header.index('Val(Tmin)')])
                    self.T_min = float(vals[header.index('Tmin [K]')])
                    self.T_max = float(vals[header.index('Tmax [K]')])

    def eval(self, T):
        """return heat capacity liquid in J/kmol/K"""
        return sum(self.constants[i]*pow(T, i-1) for i in range(1, self.num_constants + 1))

    def integral(self, T):
        return sum(self.constants[i]*pow(T, i)/i for i in range(1, self.num_constants + 1))

    def integral_dT(self, T_ref, T):
        """
        .. math::

            \\int_{Tref}^T CpL dT

        """
        return self.integral(T) - self.integral(T_ref)
