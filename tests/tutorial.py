from distillation.amundson_1958 import Model

model = Model(
             components=['n-Butane', 'n-Pentane'],
             F=1000.,              # feed flow rate [kmol/h]
             P=101325.*2,          # pressure (constant) [Pa]
             z_feed=[0.45, 0.55],  # mole fraction n-Butane = 0.45
             RR=1.,                # reflux ratio [L/D]
             D=400.,               # distillate flow rate [kmol/h]
             N=3,                  # number of equilibrium contacts
             feed_stage=2,         # stage at which feed goes in
)

model.generate_initial_guess()
model.update_K_values()
import pprint
pprint.pprint(model.K)
for i in model.components:
    model.solve_component_mass_bal(i)
pprint.pprint(model.l)
# print(model.L, 'L before update flow rates')
# print(model.V, 'V before update flow rates')
# model.update_flow_rates()
# print(model.L, 'L after update flow rates')
# print(model.V, 'V after update flow rates')
print(model.T, 'T before update')
model.update_T_values()
print(model.T, 'T after update')

print(model.T_is_converged(), 'T_is_converged after first call')
while not model.T_is_converged():
    model.update_K_values()
    for i in model.components:
        model.solve_component_mass_bal(i)
    model.update_T_values()
    print(model.T, '<= temperature at end of while loop')

# print(model.L, 'L before solve energy balances')
# print(model.V, 'V before solve energy balances')
# model.solve_energy_balances()
# print(model.L, 'L after solve energy balances')
# print(model.V, 'V after solve energy balances')
# print(model.flow_rates_converged(), 'Flow rate converged value after first iter')
#
