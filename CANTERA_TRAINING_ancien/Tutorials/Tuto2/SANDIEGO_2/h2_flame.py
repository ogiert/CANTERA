###############################################################
#
# ADIABATIC_FLAME - A freely-propagating, premixed flat flame 
#              
###############################################################

# import :

import sys
import cantera as ct
import matplotlib.pyplot as plt
import numpy as np

#################################################################
# Prepare your run
#################################################################
# Parameter values :

# General
p = 1.0e05  # pressure
tin = 300.0  # unburned gas temperature
# comp       =  'H2:2., O2:1., N2:3.76' # composition
phi_min = 0.6
phi_max = 2.6

# Composition parameters and storage
npoints = 21
phi = np.zeros(npoints, 'd')
sl_cas = np.zeros(npoints, 'd')

# Several choices for the initial grid:

# The advised choice, only the domain width is specified :
width = 0.02

# Initial grids, chosen to be 0.02cm long :
# - Refined grid at inlet and outlet, 6 points in x-direction :
initial_grid = width * np.array([0.0, 0.1, 1, 2, 2.9, 3], 'd') / 3  # m

# - Uniform grid, 6 points in x-direction (import numpy):
# initial_grid = 0.02*array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0],'d') # m

# - Uniform grid of 300 points using numpy :
# initial_grid = numpy.linspace(0,0.02 , 300)

# Tolerance properties
tol_ss = [1.0e-6, 1.0e-13]       # [rtol atol] for steady-state problem
tol_ts = [1.0e-5, 1.0e-9]        # [rtol atol] for time stepping

loglevel = 0  # amount of diagnostic output (0 to 5)

refine_grid = "refine"  # True to enable refinement, False to disable

# Import gas phases with mixture transport model
gas = ct.Solution('h2_sandiego.cti', 'gas')
# gas = Solution('gri30.cti','gri30_mix')

# You need to specify the composition of you fresh gases
# You can either specify a composition by hand :

# Initial composition of the mixture:
# gaseous fuel species
fuel = 'H2'

# air composition
air = {'O2': 1, 'N2': 3.76}

phi[0] = phi_min

# Everything can be done by hand :

# stoichio O2
stoich_O2 = 0.25 * gas.n_atoms(fuel, 'H')  # For hydrogen only

# nb species in case
m = gas.n_species
nsp_cas = m

# find fuel, nitrogen, and oxygen indices for this case
ifuel = gas.species_index(fuel)
io2 = gas.species_index('O2')
in2 = gas.species_index('N2')

# You can either specify an array
x = np.zeros(nsp_cas, 'd')
x[ifuel] = phi[0]
x[io2] = stoich_O2
x[in2] = stoich_O2 * air['N2']

# Or a dictionary
# Being for concise
x = {fuel: phi[0], 'O2': stoich_O2, 'N2': stoich_O2 * air['N2']}

gas.TPX = tin, p, x

# # You can also use the built-in function :
# gas.TP = tin, p
# gas.set_equivalence_ratio(phi[0], fuel, air)


#################################################################
# FIRST SIMULATION :
#################################################################

phi[0] = phi_min
gas.TP = tin, p
gas.set_equivalence_ratio(phi[0], fuel, air)

# Create the free laminar premixed flame
f = ct.FreeFlame(gas, initial_grid)

f.flame.set_steady_tolerances(default=tol_ss)
f.flame.set_transient_tolerances(default=tol_ts)

f.inlet.X = x
f.inlet.T = tin

#################################################################
# Program starts here
#################################################################
# First flame:

# No energy for starters
f.energy_enabled = False

# Refinement criteria
f.set_refine_criteria(ratio=10.0, slope=1, curve=1)

# Max number of times the Jacobian will be used before it must be re-evaluated
f.set_max_jac_age(50, 50)

# Set time steps whenever Newton convergence fails
f.set_time_step(1.0e-6, [2, 5, 10, 20, 80])  # s

# Calculation
f.solve(loglevel, refine_grid)

#################
# Second flame:

# Energy equation enabled
f.energy_enabled = True

# Refinement criteria when energy equation is enabled
f.set_refine_criteria(ratio=5.0, slope=1, curve=1)

# Calculation and save of the results
f.solve(loglevel, refine_grid)

#################
# Third flame :

# Refinement criteria should be changed ...
f.set_refine_criteria(ratio=3.0, slope=0.8, curve=0.9)

# Calculation
f.solve(loglevel, refine_grid)

#################
# Fourth flame and so on ...:

# Refinement criteria should be changed ...
f.set_refine_criteria(ratio=2.0, slope=0.5, curve=0.5)

# Calculation ...
f.solve(loglevel, refine_grid)

# Refinement criteria should be changed ...
f.set_refine_criteria(ratio=2.0, slope=0.1, curve=0.2, prune=0.05)

# Calculation ...
f.solve(loglevel, refine_grid)

sl_cas[0] = f.u[0]

# ...and saving of the results
f.save('h2-air_adiabatic.xml', 'energy_0',
       'solution with the energy equation enabled')

#################################################################
# Rest of the simulation :
#################################################################

# Start the loop on phis
for j in range(1, npoints):
    phi[j] = phi_min + (phi_max - phi_min) * j / (npoints - 1)

    #################
    # Assembling objects :

    # Set gas state to that of the unburned gas
    gas.TP = tin, p
    gas.set_equivalence_ratio(phi[j], fuel, air)

    # Create the free laminar premixed flame
    f = ct.FreeFlame(gas, initial_grid)
    f.restore('h2-air_adiabatic.xml', 'energy_' + str(j - 1))

    # tol_ss    = [1.0e-8, 1.0e-13]
    f.flame.set_steady_tolerances(default=tol_ss)
    f.flame.set_transient_tolerances(default=tol_ts)
    # f.set_refine_criteria(ratio = 2.0, slope = 0.1, curve = 0.5)

    f.inlet.X = gas.X
    f.inlet.T = tin
    f.P = p

    f.solve(loglevel, refine_grid)
    sl_cas[j] = f.u[0]

    # ...and saving of the results
    f.save('h2-air_adiabatic.xml', 'energy_' + str(j),
           'solution with the energy equation enabled')

    #################################################################
    # Save your results if needed
    #################################################################
    # Write the velocity, temperature, density, and mole fractions to a CSV file
    f.write_csv('h2-air_adiabatic_phi' + str(phi[j]) + '.csv', quiet=False)

#################################################################
# Plot your results
#################################################################
# import matplotlib.pylab as plt
# create plot
# fig=figure(1)

# create first subplot - adiabatic flame temperature
# a=fig.add_subplot(111)
if '--plot' in sys.argv:
    plt.plot(phi, sl_cas, '-')
    # hold(True)

    plt.title(r'Sl vs. $\Phi$ for $H_{2}/Air$ flames, at P = ' + str(p) + 'Pa, Tin = ' + str(tin) + 'K')
    plt.xlabel(r'$\Phi$', fontsize=20)
    plt.ylabel("Laminar flame speed [m/s]")
    # a.axis([0.5,1.7,0.0,3.5])
    # ax = gca()
    # ax.set_autoscale_on(False)
    # a.xaxis.set_major_locator(MaxNLocator(13)) # this controls the number of tick marks on the axis
    # a.yaxis.set_major_locator(MaxNLocator(20)) # this controls the number of tick marks on the axis
    # hold(False)
    # legend(bbox_to_anchor=(0.9, 0.9,1,1),loc=8)

    # fig.text(0.5,0.95,r'Laminar flame speed for $H_{2}$ + Air, at P = '+ str(p)+'Pa, Tin = '+str(tin)+'K.',
    # fontsize=22, horizontalalignment='center')
    plt.grid()
    # subplots_adjust(left=0.08, right=0.96, wspace=0.25)

    # show()
    plt.savefig('plot_flamespeed-' + str(tin) + '-' + str(p) + '.png', bbox_inches='tight')
