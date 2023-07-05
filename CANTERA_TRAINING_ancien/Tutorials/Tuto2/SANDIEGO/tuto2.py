import cantera as ct

gas = ct.Solution('h2_sandiego.cti')

print(gas.reaction_type(15))
