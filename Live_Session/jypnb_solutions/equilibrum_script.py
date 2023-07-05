import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

NAME = ['Mechanisms/GRI-Mech2.11.cti', 'gri30.cti']
T = 400.0                     # Temperature
P = 2.0e5                     # Pressure

phi_min = 0.3                 # Minimal equivalence ratio
phi_max = 10.0                 # Maximal equivalence ratio
npoints = 50                  # Point in-between the two preceeding values

fuel_species = 'CH4'         # fuel species
air_N2_O2_molar_ratio = 3.76  # ratio representing the air
phi = np.zeros(npoints)                  # 1D array
tad = np.zeros(npoints)                  # 1D array
for name in NAME:
    gas = ct.Solution(name)
    gas.TP = T,P
    for i in range(npoints):
        phi[i] = phi_min + (phi_max - phi_min) * i / (npoints - 1)
        gas.set_equivalence_ratio(phi[i], {fuel_species: 1}, {'O2': 1, 'N2': air_N2_O2_molar_ratio})
        gas.equilibrate('HP') # Equilibrate the mixture adiabatically at constant P with the solver vcs
        tad[i] = gas.T
        print("At phi = ","%10.4f"% (phi[i])+ "  Tad = ","%10.4f"% (tad[i]))

    plt.plot(phi, tad, label=name)
plt.title(r'Tad vs. $\Phi$ for $CH_{4}/Air$ flames')
plt.xlabel(r'$\Phi$', fontsize=20)
plt.ylabel("Adiabatic flame temperature (K)")
plt.legend()
plt.grid()
plt.show()
