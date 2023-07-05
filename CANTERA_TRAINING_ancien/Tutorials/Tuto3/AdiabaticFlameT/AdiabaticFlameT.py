"""
Adiabatic flame temperature and equilibrium composition for a fuel/air mixture
as a function of equivalence ratio.
"""

import cantera as ct
import matplotlib.pyplot as plt
import numpy as np
import sys
import csv

##############################################################################
# Edit these parameters to change the initial temperature, the pressure, and
# the phases in the mixture.

# Import the gas :
gas = ct.Solution('mech.cti')

# Choose the equivalence ratio range :
phi_min = 0.3
phi_max = 3.5
npoints = 50

# Set the gas composition :
T = 300.0
P = 101325.0

# find fuel, nitrogen, and oxygen indices
# to set the composition of the gas
fuel_species = 'C2H4'

##############################################################################

# Create some arrays to hold the data with numpy
# 1D arrays : phi, tad
phi = np.zeros(npoints)
tad = np.zeros(npoints)

# 2D arrays :
xeq = np.zeros((gas.n_species, npoints))

# Start the loop on Phi :
for i in range(npoints):
    # Start with setting the composition of the gas
    air_N2_O2_molar_ratio = 3.76

    phi[i] = phi_min + (phi_max - phi_min) * i / (npoints - 1)

    gas.set_equivalence_ratio(phi[i], {fuel_species: 1}, {'O2': 1, 'N2': air_N2_O2_molar_ratio})

    gas.TP = T, P

    # Equilibrate the mixture adiabatically at constant P
    # with the solver vcs
    gas.equilibrate('HP')

    # Save the adiabatic temperature each time
    tad[i] = gas.T
    # you can even print it on screen
    print("At phi = ","%10.4f"% (phi[i])+ "  Tad = ","%10.4f"% (tad[i]))

    # You could also save the equilibrium mass fractions at each
    # phi, as long as you initialized an array :
    xeq[:, i] = gas.X

# Save your results in a CSV file
csv_file = 'yourfile.csv'
with open(csv_file, 'w') as outfile:
    writer = csv.writer(outfile)
    writer.writerow(['Phi', 'T (K)'] + gas.species_names)
    for i in range(npoints):
        writer.writerow([phi[i], tad[i]] + list(xeq[:, i]))
print("Output written to", "%s" % csv_file)

# Add a plot option
if '--plot' in sys.argv:
    plt.figure(1)
    plt.plot(phi, tad, '-')

    plt.title(r'Tad vs. $\Phi$ for $C_{2}H_{4}/Air$ flames')
    plt.xlabel(r'$\Phi$', fontsize=20)
    plt.ylabel("Adiabatic flame temperature (K)")

    plt.grid()
   
    plt.savefig('plot_tad.png', bbox_inches='tight')
 
    plt.figure(2)
    for i, cas in enumerate(gas.species_names):
        if cas in ['O2','CO2','CO']:
            plt.plot(phi,xeq[i,:], label = cas)
    plt.xlabel('Equivalence ratio')
    plt.ylabel('Mole fractions')
    plt.legend(loc='best')
    plt.title(r'Mole fraction vs. $\Phi$ for $C_{2}H_{4}/Air$ flames')
    plt.grid()

    plt.savefig('plot_Mole_Frac.png', bbox_inches='tight')

    plt.show()
