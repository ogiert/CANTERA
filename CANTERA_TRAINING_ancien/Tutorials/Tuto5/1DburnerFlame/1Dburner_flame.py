###############################################################
#
# ADIABATIC_FLAME - A freely-propagating, premixed flat flame 
#              
###############################################################

# import :
import cantera as ct
import numpy as np

#################################################################
# Prepare your run
#################################################################

# Import gas phases with mixture transport model
gas = ct.Solution('gri30.cti')

# Parameter values :
# General
p = 1e5  # pressure
tin = 300.0  # unburned gas temperature
phi = 1

fuel = {'CH4': 1}
oxidizer = {'O2': 1, 'N2': 3.76}

initial_grid = np.linspace(0, 0.02, 300)

# Set gas state to that of the unburned gas
gas.TP = tin, p
gas.set_equivalence_ratio(phi, fuel, oxidizer)

# Create the free laminar premixed flame
f = ct.BurnerFlame(gas,initial_grid)

# set inlet conditions
f.inlet.X = gas.X
f.inlet.T = gas.T

#################################################################
# Program starts here
#################################################################
# First flame:

# No energy for starters
f.energy_enabled = False

# Tolerance properties
tol_ss = [1.0e-5, 1.0e-9]  # [rtol atol] for steady-state problem
tol_ts = [1.0e-5, 1.0e-9]  # [rtol atol] for time stepping

f.flame.set_steady_tolerances(default=tol_ss)
f.flame.set_transient_tolerances(default=tol_ts)

# Max number of times the Jacobian will be used before it must be re-evaluated
f.set_max_jac_age(50, 50)

# Set time steps whenever Newton convergence fails
f.set_time_step(1.0e-5, [2, 5, 10, 20, 80])  # s

# Refinement criteria
f.set_refine_criteria(ratio=10.0, slope=1, curve=1)

# Calculation
loglevel = 1  # amount of diagnostic output (0
# to 5)

refine_grid = False  # True to enable refinement, False to
# disable

f.solve(loglevel, refine_grid)

f.save('ch4_adiabatic.xml', 'no energy',
       'solution with no energy')

################
# Second flame:

# Energy equation enabled
f.energy_enabled = True

# Calculation and save of the results
refine_grid = True

# Refinement criteria when energy equation is enabled
f.set_refine_criteria(ratio=5.0, slope=0.5, curve=0.5)

f.solve(loglevel, refine_grid)

f.save('ch4_adiabatic.xml', 'energy',
       'solution with the energy equation enabled')

# See the sl to get an idea of whether or not you should continue
print('mixture-averaged flamespeed = {0:7f} m/s'.format(f.u[0]))  # m/s

#################
# Third flame and so on ...:

# Refinement criteria should be changed ...
f.set_refine_criteria(ratio=2.0, slope=0.05, curve=0.05)

# Calculation and saving of the results should always be done
f.solve(loglevel, refine_grid)

f.save('ch4_adiabatic.xml', 'energy continuation',
       'solution with the energy equation enabled continuation')

points = f.flame.n_points
print('mixture-averaged flamespeed continuation = {0:7f} m/s'.format(f.u[0]))  # m/s
print('mixture-averaged final T = {0:7f} K'.format(f.T[points - 1]))  # K

#################
# Fourth flame and so on ...:

# Switch transport model
f.transport_model = 'Multi'

f.solve(loglevel, refine_grid)

f.save('ch4_adiabatic.xml', 'energy_multi',
       'solution with the multicomponent transport and energy equation enabled')

points = f.flame.n_points
print('multicomponent flamespeed = {0:7f} m/s'.format(f.u[0]))  # m/s
print('multicomponent final T = {0:7f} K'.format(f.T[points - 1]))  # K

#################################################################
# Save your results if needed
#################################################################
# Write the velocity, temperature, density, and mole fractions to a CSV file
f.write_csv('ch4_adiabatic.csv', quiet=False)
