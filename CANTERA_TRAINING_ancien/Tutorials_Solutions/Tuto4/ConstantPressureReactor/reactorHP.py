"""
Constant-pressure, adiabatic kinetics simulation.
"""

import sys
import numpy as np
import cantera as ct
import matplotlib.pyplot as plt

#################################################################################

# Mechanisms used for the process
gri3 = ct.Solution('gri30.xml')
air = ct.Solution('air.xml')

# Gas state
gri3.TPX = 1000.0, ct.one_atm, 'CH4:0.5,O2:1,N2:3.76'

# Create Reactors used for the process and initialize them
r = ct.IdealGasReactor(gri3)
env = ct.Reservoir(air)

# Define a wall between the reactor and the environment, and
# make it flexible, so that the pressure in the reactor is held
# at the environment pressure.
w = ct.Wall(r, env)
w.expansion_rate_coeff = 1.0e6  # set expansion parameter. dV/dt = KA(P_1 - P_2)
w.area = 1.0

# Prepare the simulation with a ReactorNet object
sim = ct.ReactorNet([r])
time = 4.e-1

# Arrays to hold the datas
times = np.zeros(200)
data = np.zeros((200, 4))

# Advance the simulation in time
print(('%10s %10s %10s %14s' % ('t [s]', 'T [K]', 'P [Pa]', 'h [J/kg]')))
for n in range(200):
    time += 5.e-3
    sim.advance(time)
    times[n] = time  # time in s
    data[n, 0] = r.T
    data[n, 1:] = r.thermo['O2', 'H', 'CH4'].X
    print(('%10.3e %10.3f %10.3f %14.6e' % (sim.time, r.T,
                                            r.thermo.P, r.thermo.h)))

# Plot the results if matplotlib is installed.
if '--plot' in sys.argv[1:]:

    plt.clf()
    plt.subplot(2, 2, 1)
    plt.plot(times, data[:, 0])
    plt.xlabel('Time (ms)')
    plt.ylabel('Temperature (K)')
    plt.subplot(2, 2, 2)
    plt.plot(times, data[:, 1])
    plt.xlabel('Time (ms)')
    plt.ylabel('O2 Mole Fraction')
    plt.subplot(2, 2, 3)
    plt.plot(times, data[:, 2])
    plt.xlabel('Time (ms)')
    plt.ylabel('H Mole Fraction')
    plt.subplot(2, 2, 4)
    plt.plot(times, data[:, 3])
    plt.xlabel('Time (ms)')
    plt.ylabel('CH4 Mole Fraction')
    plt.tight_layout()
    plt.show()
else:
    print("To view a plot of these results, run this script with the option --plot")
