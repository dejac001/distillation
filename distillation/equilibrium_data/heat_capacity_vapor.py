class CpV:
    """Heat capacities of vapor"""
    def __init__(self, compound_name, verbose=False):
        self.name = compound_name
        self.value = 4.*8.314*1000.
        self.units = 'J/kmol/K'

        if verbose:
            print('Setting heat capacity for %s to be 4R (ideal polyatomic gas)' % compound_name)

    def eval(self):
        return self.value

    def integral_dT(self, T_ref, T):
        return self.value*(T - T_ref)

