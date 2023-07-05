"""
An opposed-flow ethane/air diffusion flame
"""

import cantera as ct
import numpy as np

##################################################################
# Create the gas object 
gas = ct.Solution('gri30.xml', 'gri30_mix')

# Parameter values :
# General
p = 1e5  # pressure

# Different input parameters
tin_f = 300.0  # fuel inlet temperature
tin_o = 300.0  # oxidizer inlet temperature
mdot_f = 0.24  # fuel inlet kg/m^2/s
mdot_o = 0.72  # oxidizer inlet kg/m^2/s

comp_f = 'C2H6:1'  # fuel composition
comp_o = 'O2:0.21, N2:0.78, AR:0.01'  # air composition

# Set gas state for the flame domain
# gas.TP = tin_o, p
# Distance between inlets is 2 cm.
# Start with an evenly-spaced 6-point grid.
initial_grid = np.linspace(0, 0.02, 6)

################ generate initial strain info #######
rho_o = p / (8.314 / 0.029 * tin_o)
rho_f = p / (8.314 / 0.030 * tin_f)

vel_o = mdot_o / rho_o
vel_f = mdot_f / rho_f

a = (vel_o + vel_f) / 0.02  # s-1

################ generate initial strain info #######

# Create an object representing the counterflow flame configuration,
# which consists of a fuel inlet on the left, the flow in the middle,
# and the oxidizer inlet on the right.
f = ct.CounterflowDiffusionFlame(gas, initial_grid)

# Set the state of the two inlets
f.fuel_inlet.mdot = mdot_f
f.fuel_inlet.X = comp_f
f.fuel_inlet.T = tin_f

f.oxidizer_inlet.mdot = mdot_o
f.oxidizer_inlet.X = comp_o
f.oxidizer_inlet.T = tin_o

# construct the initial solution estimate. To do so, it is necessary
# to specify the fuel species. If a fuel mixture is being used,
# specify a representative species here for the purpose of
# constructing an initial guess.
f.set_initial_guess(fuel='C2H6')

#################################################################
# Program starts here
#################################################################
# First flame:
# disable the energy equation
f.energy_enabled = False

# Set error tolerances
tol_ss = [1.0e-5, 1.0e-11]  # [rtol, atol] for steady-state problem
tol_ts = [1.0e-5, 1.0e-11]  # [rtol, atol] for time stepping
f.flame.set_steady_tolerances(default=tol_ss)
f.flame.set_transient_tolerances(default=tol_ts)

# and solve the problem without refining the grid
f.solve(loglevel=1, refine_grid=False)

#################
# Second flame:
# specify grid refinement criteria, turn on the energy equation,
f.energy_enabled = True

f.set_refine_criteria(ratio=3, slope=0.8, curve=0.8)

# and solve the problem again
f.solve(loglevel=1, refine_grid=True)

# save your results
f.save('c2h6_diffusion.xml', 'energy')

#################
# Third flame:
# specify grid refinement criteria, turn on the energy equation,
f.energy_enabled = True

f.set_refine_criteria(ratio=2, slope=0.2, curve=0.2, prune=0.04)

# and solve the problem again
f.solve(loglevel=1, refine_grid=True)

# save your results
f.save('c2h6_diffusion.xml', 'energy continuation')

#################################################################
# Save your results if needed
#################################################################
# write the velocity, temperature, and mole fractions to a CSV file
f.write_csv('c2h6_diffusion.csv', quiet=False)
