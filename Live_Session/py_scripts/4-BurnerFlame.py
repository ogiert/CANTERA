#!/usr/bin/env python
# coding: utf-8


import cantera as ct
import numpy as np
from matplotlib.pylab import *

# Import gas phases with mixture transport model
gas = ct.Solution('gri30.cti')

# Parameter values :
# General
p = 1e5  # pressure
tin = 373  # unburned gas temperature
phi = 1.3

fuel = 'CH4: 1'
oxidizer = 'O2:1.0, N2:3.76'

# Set gas state to that of the unburned gas
gas.TP = tin, p
gas.set_equivalence_ratio(phi, fuel, oxidizer)

f = ct.BurnerFlame(gas, width=0.2)

f.burner.T =  gas.T
f.burner.X =  gas.X              # Conditions

mdot = 0.04
f.burner.mdot = mdot

#################
#f.energy_enabled = False

tol_ss = [1.0e-5, 1.0e-9]  # [rtol atol] for steady-state
tol_ts = [1.0e-5, 1.0e-4]  # [rtol atol] for time stepping

f.flame.set_steady_tolerances(default=tol_ss)
f.flame.set_transient_tolerances(default=tol_ts)

f.set_refine_criteria(ratio=10.0, slope=1, curve=1)

f.set_max_jac_age(50, 50)

f.set_time_step(1.0e-5, [1, 2, 5, 10, 20])

loglevel = 1  # amount of diagnostic output (0 to 5)

#refine_grid = 'refine'  # True to enable refinement, False to

#f.solve(loglevel, refine_grid)

#f.save('4-Output/ch4_burner_flame.xml', 'no_energy',
#       'solution with the energy equation disabled')

#################
f.energy_enabled = True

f.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.2)

#f.solve(loglevel, refine_grid='refine')

#f.save('4-Output/ch4_burner_flame.xml', 'energy',
#       'solution with the energy equation enabled')

#################
f.transport_model = 'Multi'

#f.solve(loglevel, refine_grid='refine')

#f.save('4-Output/ch4_burner_flame.xml', 'energy_multi',
#       'solution with the energy equation enabled and multicomponent transport')

#################
f.soret_enabled = True

f.solve(loglevel, refine_grid='refine')

f.save('4-Output/ch4_burner_flame.xml', 'energy_soret',
       'solution with the energy equation enabled and multicomponent transport')

import matplotlib.pyplot as plt
plt.plot(f.grid, f.u)
plt.xlabel('grid (m)',fontsize=15)
plt.ylabel('Laminar flame speed (m/s)',fontsize=15)
plt.title('Laminar flame speed vs x-axis for an equivalence ratio of 1.3', fontsize=20)
plt.show()

f.write_csv('4-Output/ch4_burner_flame.csv', quiet=False)

