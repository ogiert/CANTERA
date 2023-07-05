###############################################################
#
# Autoignition of a methane air mixture at stoichiometry
#             and atmospheric pressure, for different 
#             initial temperature
#              
###############################################################

import cantera as ct
import csv
import numpy as np

#################################################################

# Mechanism used for the process
gas = ct.Solution('gri30.cti')

# Initial temperature, Pressure and stoichiometry
gas.TPX = 1250, ct.one_atm, 'CH4:0.5, O2:1, N2:3.76'

# Set gas state

# Specify the number of time steps and the time step size
nt = 1000
dt = 1.e-4  # s

# Storage space
mfrac = []
time = []
temperature = []
HR = []

#################################################################

# Create the batch reactor
r = ct.IdealGasReactor(gas)
# Now create a reactor network consisting of the single batch reactor
sim = ct.ReactorNet([r])

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
Autoignition = time[HR.index(max(HR))]
print('Autoignition time = ' + str(Autoignition))
# Posterity
Autoignition = Autoignition * 1000  # ms
FinalTemp = temperature[nt - 1]

#################################################################
# Save results
#################################################################
# write output CSV file for importing into Excel
csv_file = 'Phi-1_P-1_T-1250_UV.csv'
with open(csv_file, 'w') as outfile:
    writer = csv.writer(outfile)
    writer.writerow(['Auto ignition time [ms]', 'Final Temperature [K]'] + gas.species_names)
    writer.writerow([Autoignition, FinalTemp] + list(mfrac[:]))
print('output written to ' + csv_file)
