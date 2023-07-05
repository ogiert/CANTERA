#!/usr/bin/env python
# coding: utf-8

import sys
import numpy as np
import cantera as ct
import matplotlib.pyplot as plt
from matplotlib import *


gas = ct.Solution('gri30.cti')
gas.TPX = 1000.0, ct.one_atm, 'CH4:0.5,O2:1,N2:3.76'

r = ct.Reactor(gas)
sim = ct.ReactorNet([r])
time = 4.e-1

n_points = 400
times = np.zeros(n_points)
data = np.zeros((n_points, 4))

print(('%10s %10s %10s %14s' % ('t [s]', 'T [K]', 'vol [m3]', 'u [J/kg]')))
for n in range(n_points):
    time += 5.e-3
    sim.advance(time)
    times[n] = time  # time in s
    data[n, 0] = r.T                               # set the temperature in the first column
    data[n, 1:] = r.thermo['O2', 'CO2', 'CH4'].X     # set the molar fractions in the other column
    print(('%10.3e %10.3f %10.3f %14.6e' % (sim.time, r.T,
                                            r.thermo.v, r.thermo.u)))

rcParams['figure.figsize'] = (14, 10)

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
plt.ylabel('CO2 Mole Fraction')

plt.subplot(2, 2, 4)
plt.plot(times, data[:, 3])
plt.xlabel('Time (ms)')
plt.ylabel('CH4 Mole Fraction')

plt.show()

