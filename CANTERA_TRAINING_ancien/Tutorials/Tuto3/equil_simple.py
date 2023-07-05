# homogeneous equilibrium of a gas

import cantera as ct

# create an object representing the gas phase
gas = ct.Solution('gri30.cti')

compo = {}
compo['CH4'] = 0.5
compo['O2'] = 1
compo['N2'] = 3.76

# set initial state
gas.TPX = (300.0, 1.0e05, compo)

gas.equilibrate('TP')

print(gas())
