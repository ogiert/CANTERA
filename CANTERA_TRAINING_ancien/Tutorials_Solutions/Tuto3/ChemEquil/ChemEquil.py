###############################################################
#
# Equilibrium calculations
# with different solvers 
#              
###############################################################

# Cantera has 3 different equilibrium solvers, 2 of them are worth 
# mentionning: The 'ChemEquil' solver uses the element potential 
# method for homogeneous equilibrium in gas
# mixtures. It is fast, but sometimes doesn't converge. 
# The 'VCS' solver uses the VCS algorithm (Gibbs
# minimization), which is slower but more robust. 
# It can also handle multiple phases. Here we'll solve a
# problem for which the ChemEquil solver fails, but the
# VCS solver has no problem.

# import :

import cantera as ct
import numpy as np
import csv

#################################################################
# Prepare your run
#################################################################
# Parameter values :

# General
pressure = 1.0e5  # pressure
temperature = 400.0  # unburned gas temperature
comp = 'CH4:0.5, O2:1, N2:3.76'  # premixed gas composition

print("")
print("Computing Equilibirum at Phi = 1, T = " + str(temperature) + " K, P = " + str(pressure) + " Pa")
print("Equilibrate holding TP constants")
print("Using different solvers")
print("")

# Use GRI-Mech 3.0 for the methane/air mixture, and set its initial state
gas = ct.Solution('gri30.xml')

#################
# Assembling objects :

# Set gas state to that of the unburned gas
gas.TPX = temperature, pressure, comp

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

# Print the chemical potentials for further comparison
# (The element potentials are the chemical potentials of the atomic vapors.)
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

#################################################################
# Program starts here
#################################################################
# Equilibrium:
try:
    gas.equilibrate('TP', solver='element_potential')  # use the ChemEquil solver
except:
    print("")
    print("ChemEquil solver failed! Trying the vcs solver...")
    gas.equilibrate('TP', solver='vcs', maxsteps=1500)
    # gas.equilibrate('TP', solver = 'gibbs')              # the gibbs solver works also

#################################################################
# Print results
#################################################################
# On screen    
tad = gas.T
xeq = np.zeros(gas.n_species)
xeq[:] = gas.X

print("")
print("******************************************************** ")
print("    Final state :")
print("    Tadiabatique = " + str(tad) + " K")
print("******************************************************** ")
print("P =  ", "%10.4e  " % (gas.P) + "    Pa")
print("T =  ", "%10.4e  " % (gas.T) + "    K")
print("V =  ", "%10.4e  " % (gas.volume_mass) + "    m3/kg")
print("U =  ", "%10.4e  " % (gas.int_energy_mass) + "    J/kg")
print("H =  ", "%10.4e  " % (gas.enthalpy_mass) + "    J/kg")
print("S =  ", "%10.4e  " % (gas.entropy_mass) + "    J/kg/K")
print("")
print("")

# To check that this is an equilibrium state, verify that the chemical
# potentials may be computed by summing the element potentials for each atom.
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

#################################################################
# Save results
#################################################################
# Save the state after the simulation
csv_file = 'all_mole_fractions.csv'
with open(csv_file, 'w') as outfile:
    writer = csv.writer(outfile)
    writer.writerow(['phi', 'T (K)'] + gas.species_names)
    writer.writerow(['1', tad] + list(xeq[:]))
print(('Output written to {0}'.format(csv_file)))
