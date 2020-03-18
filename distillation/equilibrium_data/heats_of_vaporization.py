class dH_vap:
    """Heat of vaporization"""
    def __init__(self, compound_name):
        from distillation import os, ROOT_DIR
        file = os.path.join(ROOT_DIR, 'equilibrium_data', 'heats_of_vaporization.csv')
        self.units = 'J/kmol'
        with open(file, 'r') as f:
            header = next(f).rstrip('\n').split(',')
            for line in f:
                vals = line.rstrip('\n').split(',')
                if vals[0] == compound_name:
                    self.value = float(vals[header.index('Value')]) * 1e6   # convert from kJ/mol to J/kmol
                    self.T_ref = float(vals[header.index('T_ref')])

    def eval(self):
        return self.value
