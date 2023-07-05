###############################################################
#
# ADIABATIC_FLAME - A freely-propagating, premixed flat flame 
#              
###############################################################

# import :
import cantera as ct

#################################################################
# Prepare your run
#################################################################

# Import gas phases with mixture transport model
gas = ct.Solution('gri30.cti')

# Parameter values :
# General
p = 1e5  # pressure
tin = 600.0  # unburned gas temperature
phi = 0.8

fuel = {'CH4': 1}
oxidizer = {'O2': 1, 'N2': 3.76}

# Set gas state to that of the unburned gas
gas.TP = tin, p
gas.set_equivalence_ratio(phi, fuel, oxidizer)

# Create the free laminar premixed flame
f = ct.FreeFlame(gas, width=0.02)

#################################################################
# Program starts here
#################################################################

f.restore('ch4_adiabatic.xml', 'energy continuation')

# Set gas state to that of the unburned gas
gas.TP = tin, p
gas.set_equivalence_ratio(phi, fuel, oxidizer)

# set inlet conditions
f.inlet.X = gas.X
f.inlet.T = gas.T

f.solve()

points = f.flame.n_points
print('mixture-averaged flamespeed continuation = {0:7f} m/s'.format(f.u[0]))  # m/s
print('mixture-averaged final T = {0:7f} K'.format(f.T[points - 1]))  # K

#################################################################
# Save your results if needed
#################################################################
# Write the velocity, temperature, density, and mole fractions to a CSV file
f.write_csv('ch4_adiabatic.csv', quiet=False)
