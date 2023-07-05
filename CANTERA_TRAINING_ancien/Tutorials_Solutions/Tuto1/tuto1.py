import cantera as ct

gas = ct.Solution('gri30.xml')

print(gas())

print(gas.species_names)


