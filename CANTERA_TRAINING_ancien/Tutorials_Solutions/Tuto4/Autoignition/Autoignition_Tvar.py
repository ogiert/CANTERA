###############################################################
#
# Autoignition of a methane air mixture at stoichiometry
#             and atmospheric pressure, for different
#             initial temperature
#
###############################################################

import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

#################################################################

# Mechanism used for the process
gas = ct.Solution('gri30.cti')

# Initial temperatures
Temperature_range = list(range(800, 1700, 100))

# Specify the number of time steps and the time step size
nt = 100000
dt = 1.e-4  # s

# Storing auto ignitions
auto_ignitions = []

for index, Temperature in enumerate(Temperature_range):
    #################################################################

    # Initial temperature, Pressure and stoichiometry
    gas.TPX = Temperature, ct.one_atm, 'CH4:0.5, O2:1, N2:3.76'

    # Create the batch reactor
    r = ct.IdealGasReactor(gas)
    # Now create a reactor network consisting of the single batch reactor
    sim = ct.ReactorNet([r])

    # Storage space
    mfrac = []
    # ...
    time = []
    temperature = []
    HR = []

    # Run the simulation
    # Initial simulation time
    current_time = 0.0

    # Loop for nt time steps of dt seconds.
    for n in range(nt):
        current_time += dt
        sim.advance(current_time)

        time.append(current_time)
        temperature.append(r.T)
        mfrac.append(r.thermo.Y)
        HR.append(- np.dot(gas.net_production_rates, gas.partial_molar_enthalpies))

    #################################################################
    # Catch the autoignition timing
    #################################################################

    # Get the ignition delay time by the maximum value of the Heat Release rate
    auto_ignition = time[HR.index(max(HR))]
    print('For T = ' + str(Temperature) + ', Autoignition time = ' + str(auto_ignition) + ' s')
    # Posterity
    FinalTemp = temperature[nt - 1]

    auto_ignitions.append(auto_ignition)

    # #################################################################
    # # Save results
    # #################################################################
    # # write output CSV file for importing into Excel
    # csv_file = 'Phi-1_P-1_T-' + str(Temperature) + '_UV.csv'
    # with open(csv_file, 'w') as outfile:
    #     writer = csv.writer(outfile)
    #     writer.writerow(['Auto ignition time [s]', 'Final Temperature [K]'] + gas.species_names)
    #     writer.writerow([auto_ignition, FinalTemp] + list(mfrac[:]))
    # print('output written to ' + csv_file)

T_invert = [1000 / Temperature for Temperature in Temperature_range]
#################################################################
# Plot results
#################################################################
# create plot
plt.plot(Temperature_range, auto_ignitions, 'b-o')
plt.xlabel(r'Temperature [K]', fontsize=20)
plt.ylabel("Auto ignition [s]", fontsize=20)
plt.yscale('log')
plt.title(r'Autoignition of $CH_{4}$ + Air mixture at $\Phi$ = 1, and P = 1 bar',
          fontsize=22, horizontalalignment='center')
plt.axis(fontsize=20)
plt.grid()
plt.savefig('Phi-1_P-1_Trange_UV.png', bbox_inches='tight')
plt.show()
