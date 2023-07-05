#!/usr/bin/env python
# coding: utf-8


import cantera as ct
import numpy as np
from matplotlib.pylab import *

gas = ct.Solution('gri30.xml', 'gri30_mix')

p = 1e5  # pressure

comp_f = 'C2H6:1'                       # fuel composition
tin_f = 300.0                           # fuel inlet temperature

rho_f = p / (8.314 / 0.030 * tin_f)     # fuel inlet density
mdot_f = 0.24                           # fuel inlet mass flow rate (kg/m^2/s)
vel_f = mdot_f / rho_f                  # fuel inlet velocity
print('Velocity of the fuel : ' + str(vel_f))

comp_o = 'O2:0.21, N2:0.78, AR:0.01'    # oxidizer composition
tin_o = 300.0                           # oxidizer inlet temperature

rho_o = p / (8.314 / 0.029 * tin_o)     # oxidizer inlet density
mdot_o = 0.72                           # oxidizer inlet mass flow rate (kg/m^2/s)
vel_o = mdot_o / rho_o                  # oxidizer inlet velocity
print('Velocity of the oxidizer : ' + str(vel_o))

gas.TP = tin_o, p

width = 0.02
a = (vel_o + vel_f) / width              # calculation of the strain rate (s^-1)
print('Strain rate of the diffusion flame : ' + str(a))


# This will instanciate the counterflow diffusion flame. As you might notice, and this is valid for all the other cases, you can instanciate the initial grid yourself. This allows you to refine the grid when you want. This can be a good trick when your simulation is crashing and you want to compute it again.<br><br>
# As far as the inlets are concerned, the composition, temperature and mass flow rate needs to be specified.

# In[22]:


initial_grid = np.linspace(0, width, 6)
f = ct.CounterflowDiffusionFlame(gas, initial_grid)

# Set the state of the two inlets
f.fuel_inlet.mdot = mdot_f
f.fuel_inlet.X = comp_f
f.fuel_inlet.T = tin_f

f.oxidizer_inlet.mdot = mdot_o
f.oxidizer_inlet.X = comp_o
f.oxidizer_inlet.T = tin_o


# This starts the calculation of the counterflow diffusion flame. By experience, you will see that, out of this wonderful tutorial case, diffusion flames are hard to converge. You can either play with the size of the domain and every parameter said all along the tutorial (indeed, sometimes it converges better if the domain is smaller), or you can try a method that is developed here in CERFACS that resolves **the diffusion flame in the z-space**.<br>
# If you want some more information about that, come and ask the chemistry team of cerfacs.

# In[23]:


#%%capture
# First flame:
# disable the energy equation
f.energy_enabled = False

# Set error tolerances
tol_ss = [1.0e-5, 1.0e-11]  # [rtol, atol] for steady-state problem
tol_ts = [1.0e-5, 1.0e-11]  # [rtol, atol] for time stepping
f.flame.set_steady_tolerances(default=tol_ss)
f.flame.set_transient_tolerances(default=tol_ts)

# and solve the problem without refining the grid
f.solve(loglevel=1, refine_grid='disabled')

#################
# Second flame:
# specify grid refinement criteria, turn on the energy equation,
f.energy_enabled = True

f.set_refine_criteria(ratio=3, slope=0.8, curve=0.8)

# and solve the problem again
f.solve(loglevel=1, refine_grid='refine')

# save your results
f.save('4-Output/c2h6_diffusion.xml', 'energy')

#################
# Third flame:
# specify grid refinement criteria, turn on the energy equation,
f.energy_enabled = True

f.set_refine_criteria(ratio=2, slope=0.2, curve=0.2, prune=0.04)

# and solve the problem again
f.solve(loglevel=1, refine_grid='refine')

# save your results
f.save('4-Output/c2h6_diffusion.xml', 'energy continuation')

#################################################################
# Save your results if needed
#################################################################
# write the velocity, temperature, and mole fractions to a CSV file
f.write_csv('4-Output/c2h6_diffusion.csv', quiet=False)


# These are interesting quantities to plot, such as the different mass fractions and the temperature of the flame. This is all done in the following graph (with the values being adimensionalised).

# In[24]:


# Get interesting values
z = f.flame.grid
T = f.T
u = f.u

# Get interesting indices for computation of species
fuel_species = 'C2H6'
ifuel = gas.species_index(fuel_species)
io2 = gas.species_index('O2')
in2 = gas.species_index('N2')

# Initiate interesting vectors
c2h6 = np.zeros(f.flame.n_points,'d')
o2 = np.zeros(f.flame.n_points,'d')
hr = np.zeros(f.flame.n_points,'d')

# Computes interesting quantities for analyzing a counter-flow flame
for n in range(f.flame.n_points):
    f.set_gas_state(n)
    c2h6[n]= gas.Y[ifuel]
    o2[n]= gas.Y[io2]
    hr[n] = - np.dot(gas.net_production_rates, gas.partial_molar_enthalpies)


# In[25]:


fig=figure(1)

a=fig.add_subplot(111)
a.plot(z,T/np.max(T),z,c2h6/np.max(c2h6),z,o2/np.max(o2))
plt.title(r'$T_{adiabatic}$ vs. Position',fontsize=25)
plt.xlabel(r'Position [m]', fontsize=15)
plt.ylabel('Normalized values of different quantities',fontsize=15)
plt.legend(['Temperature','$Y_{C_2H_6}$', '$Y_{O_2}$'],fontsize=15)
show()

