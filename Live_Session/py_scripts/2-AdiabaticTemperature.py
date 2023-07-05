#!/usr/bin/env python
# coding: utf-8
#
import cantera as ct
import numpy as np
import csv
from matplotlib import *
import matplotlib.pyplot as plt
import sys

gas = ct.Solution('Mechanisms/mech.cti')

T = 300.0                     # Temperature
P = 101325.0                  # Pressure

phi_min = 0.3                 # Minimal equivalence ratio
phi_max = 3.5                 # Maximal equivalence ratio
npoints = 50                  # Point in-between the two preceeding values

fuel_species = 'C2H4'         # fuel species
air_N2_O2_molar_ratio = 3.76  # ratio representing the air


phi = np.zeros(npoints)                  # 1D array
tad = np.zeros(npoints)                  # 1D array

xeq = np.zeros((gas.n_species, npoints)) # 2D array

for i in range(npoints):

    gas.TP = T, P
    
    phi[i] = phi_min + (phi_max - phi_min) * i / (npoints - 1)
    gas.set_equivalence_ratio(phi[i], {fuel_species: 1}, {'O2': 1, 'N2': air_N2_O2_molar_ratio})

    gas.equilibrate('HP') # Equilibrate the mixture adiabatically at constant P with the solver vcs
    
    xeq[:, i] = gas.X
    tad[i] = gas.T
    print("At phi = ","%10.4f"% (phi[i])+ "  Tad = ","%10.4f"% (tad[i]))


csv_file = '2-Output/yourfile.csv'
with open(csv_file, 'w') as outfile:
    writer = csv.writer(outfile)
    writer.writerow(['Phi', 'T (K)'] + gas.species_names)
    for i in range(npoints):
        writer.writerow([phi[i], tad[i]] + list(xeq[:, i]))
print("Output written to", "%s" % csv_file)

rcParams['figure.figsize'] = (14, 10)

plt.plot(phi, tad, '-')

plt.title(r'Tad vs. $\Phi$ for $C_{2}H_{4}/Air$ flames')
plt.xlabel(r'$\Phi$', fontsize=20)
plt.ylabel("Adiabatic flame temperature (K)")

plt.grid()

plt.show()
plt.savefig('2-Output/plot_tad.png', bbox_inches='tight')

for i, cas in enumerate(gas.species_names):
    if cas in ['O2','CO2','CO']:
        plt.plot(phi,xeq[i,:], label = cas)
        plt.xlabel('Equivalence ratio')
        plt.ylabel('Mole fractions')
        plt.legend(loc='best')
        plt.show()

