#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright (c) 2022-2023 Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Default simulation configuration parameters.
'''

# ....................{ IMPORTS                            }....................
from beartype import beartype
from beartype.typing import Optional
from dataclasses import dataclass
import numpy as np

@beartype
@dataclass
class CalculionParams(object):
    '''
    Configuration parameters used in Calculion simulations.

    '''
    # Constants (never change):
    R: float = 8.314  # Ideal gas constant [J/(K*mol)]
    F: float = 96485.0  # Faraday's Constant [C/mol]
    e_o: float = 8.3145e-12  # Electric permittivity free space [F/m]
    z_Na: int = 1
    z_K: int = 1
    z_Cl: int = -1
    #---------
    # Solver properties and options:
    update_env: bool = False # Update the env. concs. (True) or assume cell surrounded by large media bath (False)?
    quasi_static_vmem: bool = True # Use the quasi-static Vmem approximation?
    iterative_solver: bool = False # Use the iterative solver (True) or solve an optimization problem (False)?
    steady_state_tol: float = 1.0e-7 # Convergence tolerance to use with iterative solver (change below taken to be ss).
    delta_t: float = 10.0 # Time step [s] for iterative solver
    N_iter: int = 50000 # Maximum number of iterative time-steps to run to seek steady-state
    target_Vmem = -10.0e-3 # Target Vmem in V

    # Cell and environmental properties and options:
    T_C: float = 37.0 # Temperature [degrees C]
    r_cell_um: float = 5.0 # cell radius [um]
    d_ecm_um: float = 1.0 # Excellular space width [um]
    d_mem_nm: float = 5.0 # Plasma membrane thickness [nm]
    e_r: float = 10.0 # cell membrane relative electric permittivity

    # FIXME: later make these selectable user profiles (e.g. "mammalian neuron", "xenopus", "squid axon" etc)
    P_Na_nm = 0.1 # Membrane permeability to sodium ions [nm/s]
    P_K_nm = 1.5 # Membrane permeability to potassium ions [nm/s]
    P_Cl_nm = 0.5 # Membrane permeability to chloride ions [nm/s]

    Na_o: float = 145.0 # Extracellular sodium ion concentration [mol/m^3 or mmol/L]
    Na_i: float = 10.0 # Intracellular sodium ion concentration [mol/m^3 or mmol/L]
    K_o: float = 5.0 # Extracellular potassium ion concentration [mol/m^3 or mmol/L]
    K_i: float = 139.0 # Intracellular potassium ion concentration [mol/m^3 or mmol/L]
    Cl_o: float = 115.0 # Extracellular chloride ion concentration [mol/m^3 or mmol/L]
    Cl_i: float = 5.0 # Intracellular chloride ion concentration [mol/m^3 or mmol/L]

    # Metabolic parameters:
    delGo_ATP: float = -32e3  # Standard Gibbs Free energy ATP hydrolysis [J/mol]
    ATP: float = 4.0  # ATP concentration range 1.0 to 5.0 [mol/m^3 or mmol/L]
    ADP: float = 0.01  # ADP concentration range 0.004 to 0.1 [mol/m^3 or mmol/L]
    P: float = 0.01  # ADP concentration range 0.004 to 0.1  [mol/m^3 or mmol/L]

    # Pump and transporter parameters:
    omega_NaK: float = 5.0e-13 # Rate constant for the Na-K-ATPase ion pump
    omega_NaKCl: float=1.0e-15 # Rate constant for the Na-K-2Cl cotransporter
    omega_KCl: float=1.0e-15 # Rate constant for the K-Cl symporter

    # Calculated parameters (used internally in calculations):
    r_cell = r_cell_um*1e-6 # radius of the cell [m]
    d_ecm = d_ecm_um*1e-6 # thickness of the extracellular space [m]
    r_env = r_cell + d_ecm # radius of the surrounding environment [m]
    d_mem = d_mem_nm*1e-9 # cell membrane thickness [m]
    c_mem = (e_o*e_r)/d_mem # Cell membrane patch capacitance [F/m^2]
    T = 273.15 + T_C  # Temperature [degrees K]
    alpha = F / (R * T)  # Alpha Constant [C/J]
    P_Na = P_Na_nm*1.0e-9 # Membrane permeability to sodium ions [m/s]
    P_K = P_K_nm*1.0e-9 # Membrane permeability to potassium ions [m/s]
    P_Cl = P_Cl_nm*1.0e-9 # Membrane permeability to chloride ions [m/s]
    vol_cell = np.pi*r_cell**3
    vol_env = np.pi*(r_env**3 - r_cell**3)

    Keqm_NaK = np.exp(delGo_ATP / (R * T)) # Equilibrium constant for ATP hydrolysis reaction

    # Display parameters:
    decimals_in_rounding: int = 2







