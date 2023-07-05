#!/usr/bin/env python
# coding: utf-8

import cantera as ct
import numpy as np
from matplotlib.pylab import *

gas = ct.Solution('gri30.cti')                  # Import gas phases with mixture transport model

# General
p = 1e5                                         # pressure
tin = 600.0                                     # unburned gas temperature
phi = 0.8                                       # equivalence ratio

fuel = {'CH4': 1}                               # Methane composition
oxidizer = {'O2': 1, 'N2': 3.76}                # Oxygen composition

gas.TP = tin, p
gas.set_equivalence_ratio(phi, fuel, oxidizer)

f = ct.FreeFlame(gas, width=0.02)   # Create the free laminar premixed flame specifying the width of the grid
f.inlet.X = gas.X                   # Inlet condition on mass fraction
f.inlet.T = gas.T                   # Inlet condition on temperature

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

f.save('4-Output/ch4_adiabatic.xml', 'no energy',
       'solution with no energy')

################
# Second flame:

# Energy equation enabled
f.energy_enabled = True

# Calculation and save of the results
refine_grid = 'refine'

# Refinement criteria when energy equation is enabled
f.set_refine_criteria(ratio=5.0, slope=0.5, curve=0.5)

f.solve(loglevel, refine_grid)

f.save('4-Output/ch4_adiabatic.xml', 'energy',
       'solution with the energy equation enabled')

# See the sl to get an idea of whether or not you should continue
print('mixture-averaged flamespeed = {0:7f} m/s'.format(f.u[0]))  # m/s

#################
# Third flame and so on ...:

# Refinement criteria should be changed ...
f.set_refine_criteria(ratio=2.0, slope=0.05, curve=0.05)

# Calculation and saving of the results should always be done
f.solve(loglevel, refine_grid)

f.save('4-Output/ch4_adiabatic.xml', 'energy continuation',
       'solution with the energy equation enabled continuation')

points = f.flame.n_points
print('mixture-averaged flamespeed continuation = {0:7f} m/s'.format(f.u[0]))  # m/s
print('mixture-averaged final T = {0:7f} K'.format(f.T[points - 1]))  # K

#################
# Fourth flame and so on ...:

# Switch transport model
f.transport_model = 'Multi'

f.solve(loglevel, refine_grid)

f.save('4-Output/ch4_adiabatic.xml', 'energy_multi',
       'solution with the multicomponent transport and energy equation enabled')

points = f.flame.n_points
print('multicomponent flamespeed = {0:7f} m/s'.format(f.u[0]))  # m/s
print('multicomponent final T = {0:7f} K'.format(f.T[points - 1]))  # K

f.write_csv('4-Output/ch4_adiabatic.csv', quiet=False)  # Write the velocity, temperature, density, and mole fractions 

rcParams['figure.figsize'] = (14, 10)

# Get the different arrays
z = f.flame.grid
T = f.T
u = f.u
ifuel = gas.species_index('CH4')

fig=figure(1)

# create first subplot - adiabatic flame temperature
a=fig.add_subplot(221)
a.plot(z,T)
title(r'$T_{adiabatic}$ vs. Position', fontsize=15)
xlabel(r'Position [m]')
ylabel("Adiabatic Flame Temperature [K]")
a.xaxis.set_major_locator(MaxNLocator(10)) # this controls the number of tick marks on the axis

# create second subplot - velocity
b=fig.add_subplot(222)
b.plot(z,u)
title(r'Velocity vs. Position', fontsize=15)
xlabel(r'Position [m]')
ylabel("velocity [m/s]")
b.xaxis.set_major_locator(MaxNLocator(10)) 

# create third subplot - rho
c=fig.add_subplot(223)
p = zeros(f.flame.n_points,'d')
for n in range(f.flame.n_points):
    f.set_gas_state(n)
    p[n]= gas.density_mass
c.plot(z,p)
title(r'Rho vs. Position', fontsize=15)
xlabel(r'Position [m]')
ylabel("Rho [$kg/m^3$]")
c.xaxis.set_major_locator(MaxNLocator(10)) 


# create fourth subplot - specie CH4
d=fig.add_subplot(224)
ch4 = zeros(f.flame.n_points,'d')
for n in range(f.flame.n_points):
    f.set_gas_state(n)
    ch4[n]= gas.Y[ifuel]
d.plot(z,ch4)
title(r'CH4 vs. Position', fontsize=15)
xlabel(r'Position [m]')
ylabel("CH4 Mole Fraction")
d.xaxis.set_major_locator(MaxNLocator(10))

# Set title
fig.text(0.5,0.95,r'Adiabatic $CH_{4}$ + Air Free Flame at Phi = 0.8 Ti = 600K and P = 1atm',fontsize=22,horizontalalignment='center')

subplots_adjust(left=0.08, right=0.96, wspace=0.25, hspace=0.25)
show()

