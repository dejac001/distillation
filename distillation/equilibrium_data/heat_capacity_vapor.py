class CpV:
    """Heat capacities of vapor"""
    def __init__(self):
        self.value = 4.*8.314*1000.
        self.units = 'J/kmol/K'

    def eval(self):
        return self.value

    def integral_dT(self, T_ref, T):
        return self.value*(T - T_ref)

