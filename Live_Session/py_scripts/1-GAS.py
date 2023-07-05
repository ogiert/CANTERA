#!/usr/bin/env python
# coding: utf-8

# Import
import cantera as ct

# Solution
gas = ct.Solution("gri30.xml")

help(gas)

# Set gas properties (version 1)

gas.TPX = 500, 101325, {'CH4':1, 'O2':2, 'N2':7.52}

print(gas())

# Set gas properties (version 2)

gas.set_equivalence_ratio(1, 'CH4: 1', 'O2:1.0, N2:3.76')

print(gas())

# Interesting features

print(gas.species_names)

print(gas.species_index('CO2'))

print(gas.net_production_rates)

#print(gas.binary_diff_coeffs)

# ck2cti function

import subprocess
result = subprocess.check_output('ck2cti --input=CK2CTI/mech.inp --thermo=CK2CTI/therm.dat --transport=CK2CTI/tran.dat', shell=True)

