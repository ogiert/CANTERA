#!/usr/bin/env python
# coding: utf-8

import sys
import numpy as np
import cantera as ct
import matplotlib.pyplot as plt
from matplotlib import *

gas_a = ct.Solution('air.xml')
gas_a.TPX = 300, ct.one_atm, 'O2:0.21, N2:0.78, AR:0.01'
mw_a = gas_a.mean_molecular_weight/1000 #kg/mol

gas_b = ct.Solution('gri30.xml')
gas_b.TPX = 300.0, ct.one_atm, 'CH4:1'
mw_b = gas_b.mean_molecular_weight/1000 #kg/mol

res_a = ct.Reservoir(gas_a)
res_b = ct.Reservoir(gas_b)
downstream = ct.Reservoir(gas_b)

gas_b.TPX = 300, ct.one_atm, 'O2:1., N2:3.78, CH4:0.5'
mixer = ct.IdealGasReactor(gas_b, energy='on')

mfca = ct.MassFlowController(res_a, mixer, mdot=mw_a*2./0.21)
mfcb = ct.MassFlowController(res_b, mixer, mdot=mw_b*1.0)

outlet = ct.Valve(mixer, downstream, K=10.0)

sim = ct.ReactorNet([mixer])

print('{0:>14s} {1:>14s} {2:>14s} {3:>14s} {4:>14s}'.format('t [s]', 'T [K]', 'h [J/kg]', 'P [Pa]', 'X_CH4'))
t = 0.0
for n in range(100):
    tres = mixer.mass/(mfca.mdot(t) + mfcb.mdot(t))
    t += 2*tres
    sim.advance(t)
    print('{0:14.5g} {1:14.5g} {2:14.5g} {3:14.5g} {4:14.5g}'.format(t, mixer.T, mixer.thermo.h, mixer.thermo.P, mixer.thermo['CH4'].X[0]))

print(mixer.thermo.report())

