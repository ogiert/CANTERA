import cantera as ct
import numpy as np
from matplotlib.pylab import *

names = ['gri30.cti', 'Lu.cti']
fuels = ['CH4','CH4']
fig=figure(1)
for i,n in enumerate(names):
    gas = ct.Solution(n)                            # Import gas phases with mixture transport model

    p = 2e5                                         # pressure
    tin = 400.0                                     # unburned gas temperature
    phi = 1.1                                       # equivalence ratio

    fuel = {fuels[i]: 1}                            # Fuel composition
    oxidizer = {'O2': 1, 'N2': 3.76}                # Oxygen composition
    
    gas.TP = tin, p
    gas.set_equivalence_ratio(phi, fuel, oxidizer)
    
    f = ct.FreeFlame(gas, width=0.02)   # Create the free laminar premixed flame specifying the width of the grid
    f.inlet.X = gas.X                   # Inlet condition on mass fraction
    f.inlet.T = gas.T                   # Inlet condition on temperature
    
    f.energy_enabled = False                       # No energy for starters

    tol_ss = [1.0e-5, 1.0e-9]  # tolerance [rtol atol] for steady-state problem
    tol_ts = [1.0e-5, 1.0e-9]  # tolerance [rtol atol] for time stepping

    f.flame.set_steady_tolerances(default=tol_ss)
    f.flame.set_transient_tolerances(default=tol_ts)

    # Set calculation parameters
    f.set_max_jac_age(50, 50)                           # Maximum number of times the Jacobian will be used before it 
                                                        # must be re-evaluated
    f.set_time_step(1.0e-5, [2, 5, 10, 20, 80])         # Set time steps (in seconds) whenever Newton convergence fails 
    f.set_refine_criteria(ratio=10.0, slope=1, curve=1) # Refinement criteria

    # Calculation
    loglevel = 1                                        # amount of diagnostic output (0 to 5)
    refine_grid = 'disabled'                            # 'refine' or 'remesh' to enable refinement
                                                    # 'disabled' to disable

    f.solve(loglevel, refine_grid)                                      # solve the flame on the grid

    f.energy_enabled = True                                 # Energy equation enabled
    refine_grid = 'refine'                                 # Calculation and save of the results

    f.set_refine_criteria(ratio=5.0, slope=0.5, curve=0.5)  # Refinement criteria when energy equation is enabled

    f.solve(loglevel, refine_grid)

    f.set_refine_criteria(ratio=2.0, slope=0.05, curve=0.05) # Refinement criteria should be changed ...
    f.solve(loglevel, refine_grid)                           


    rcParams['figure.figsize'] = (14, 10)

    # Get the different arrays
    z = f.flame.grid
    u = f.u
    ifuel = gas.species_index(fuels[i])

    # create second subplot - velocity
    b=fig.add_subplot(111)
    b.plot(z,u)
    title(r'Velocity vs. Position', fontsize=25)
    xlabel(r'Position [m]', fontsize=15)
    ylabel("velocity [m/s]", fontsize=15)
    b.xaxis.set_major_locator(MaxNLocator(10)) 
plt.show()
fig.savefig('Comparison.png')
