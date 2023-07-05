#!/usr/bin/env python
# coding: utf-8


import cantera as ct
import numpy as np
import csv
from matplotlib import *
import matplotlib.pyplot as plt
import sys



gas = ct.Solution('gri30.cti')           # create an object representing the gas phase

gas.TPX = 300, 100000,{'CH4':1, 'O2':2, 'N2':7.52}      # set initial state

gas.equilibrate('TP')                 # equilibrate using Temperature (T) and Pressure (P)

print(gas())



pressure = 1.0e5                # pressure
temperature = 400.0             # unburned gas temperature
comp = 'CH4:0.5, O2:1, N2:3.76' # premixed gas composition

gas = ct.Solution('gri30.xml')
gas.TPX = temperature, pressure, comp


# Initial state of the gas


print("******************************************************** ")
print("    Initial state :")
print("******************************************************** ")
print("P =  ", "%10.4e  " % (gas.P) + "    Pa")
print("T =  ", "%10.4e  " % (gas.T) + "    K")
print("V =  ", "%10.4e  " % (gas.volume_mass) + "    m3/kg")
print("U =  ", "%10.4e  " % (gas.int_energy_mass) + "    J/kg")
print("H =  ", "%10.4e  " % (gas.enthalpy_mass) + "    J/kg")
print("S =  ", "%10.4e  " % (gas.entropy_mass) + "    J/kg/K")
print("")
print("")


# Comparing chemical potentials and element potentials

chemeq = np.zeros(gas.n_species)
chemeq = gas.chemical_potentials

mu_H2 = chemeq[gas.species_index("H2")]
mu_OH = chemeq[gas.species_index("OH")]
mu_H2O = chemeq[gas.species_index("H2O")]
mu_O2 = chemeq[gas.species_index("O2")]
lambda_H = chemeq[gas.species_index("H")]
lambda_O = chemeq[gas.species_index("O")]

print()
print("Comparison between Chem potentials and element potentials:")
print()
s_mu_H2 = "%11.4e" % mu_H2
s_lam_mu_H2 = "%11.4e" % (2.0 * lambda_H)
print("mu_H2   : ", s_mu_H2, ",    ", s_lam_mu_H2)

s_mu_O2 = "%11.4e" % mu_O2
s_lam_mu_O2 = "%11.4e" % (2.0 * lambda_O)
print("mu_O2   : ", s_mu_O2, ",    ", s_lam_mu_O2)

s_mu_OH = "%11.4e" % mu_OH
s_lam_mu_OH = "%11.4e" % (lambda_H + lambda_O)
print("mu_OH   : ", s_mu_OH, ",    ", s_lam_mu_OH)

s_mu_H2O = "%11.4e" % mu_H2O
s_lam_mu_H2O = "%11.4e" % (2.0 * lambda_H + lambda_O)
print("mu_H2O  : ", s_mu_H2O, ",    ", s_lam_mu_H2O)


# Program equilibrate

try:
    # print("0")
    gas.equilibrate('TP', solver='element_potential')  # use the ChemEquil solver
except:
    print("")
    print("ChemEquil solver failed! Trying the vcs solver...")
    gas.equilibrate('TP', solver='vcs', maxsteps=1500)
    # gas.equilibrate('TP', solver = 'gibbs')              # the gibbs solver works also


# Compare the results with the initial values


print("")
print("******************************************************** ")
print("    Final state :")
print("    Tadiabatique = " + str(gas.T) + " K")
print("******************************************************** ")
print("P =  ", "%10.4e  " % (gas.P) + "    Pa")
print("T =  ", "%10.4e  " % (gas.T) + "    K")
print("V =  ", "%10.4e  " % (gas.volume_mass) + "    m3/kg")
print("U =  ", "%10.4e  " % (gas.int_energy_mass) + "    J/kg")
print("H =  ", "%10.4e  " % (gas.enthalpy_mass) + "    J/kg")
print("S =  ", "%10.4e  " % (gas.entropy_mass) + "    J/kg/K")
print("")
print("")


# Comparing chemical and element equilibrium for the equilibrate mixture

chemeq = gas.chemical_potentials
mu_H2 = chemeq[gas.species_index("H2")]
mu_OH = chemeq[gas.species_index("OH")]
mu_H2O = chemeq[gas.species_index("H2O")]
mu_O2 = chemeq[gas.species_index("O2")]
lambda_H = chemeq[gas.species_index("H")]
lambda_O = chemeq[gas.species_index("O")]

print()
print("Comparison between Chem potentials and element potentials:")
print()
s_mu_H2 = "%11.4e" % mu_H2
s_lam_mu_H2 = "%11.4e" % (2.0 * lambda_H)
print("mu_H2   : ", s_mu_H2, ",    ", s_lam_mu_H2)

s_mu_O2 = "%11.4e" % mu_O2
s_lam_mu_O2 = "%11.4e" % (2.0 * lambda_O)
print("mu_O2   : ", s_mu_O2, ",    ", s_lam_mu_O2)

s_mu_OH = "%11.4e" % mu_OH
s_lam_mu_OH = "%11.4e" % (lambda_H + lambda_O)
print("mu_OH   : ", s_mu_OH, ",    ", s_lam_mu_OH)

s_mu_H2O = "%11.4e" % mu_H2O
s_lam_mu_H2O = "%11.4e" % (2.0 * lambda_H + lambda_O)
print("mu_H2O  : ", s_mu_H2O, ",    ", s_lam_mu_H2O)


# Saving the results

csv_file = '2-Output/all_mole_fractions.csv'
with open(csv_file, 'w') as outfile:
    writer = csv.writer(outfile)
    writer.writerow(['phi', 'T (K)'] + gas.species_names)
    writer.writerow(['1', gas.T] + list(gas.X))
print(('Output written to {0}'.format(csv_file)))





