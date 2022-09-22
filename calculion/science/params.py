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

R = 8.314 # Ideal gas constant J/(K*mol)
F = 96485 # Faraday's Constant
e_o = 8.3145e-12 # Electric permittivity free space [F/m]

@beartype
@dataclass
class CalculionParams(object):
    '''
    Configuration parameters used in Calculion simulations.

    '''
    T_C: float = 37.0 # Temperature [degrees C]
    r_cell: float = 5.0e-6 # radius of the cell [m]
    r_env: Optional[float] = None # radius of the surrounding environment; if None no env updates [m]
    d_mem: float = 5.0e-9 # cell membrane thickness [m]
    e_r: float = 10.0 # cell membrane relative electric permittivity

    T = 273.15 + T_C  # Temperature [degrees K]
    alpha = F / (R * T)  # Alpha Constant [C/J]

    # FIXME: later make these selectable user profiles (e.g. "mammalian neuron", "xenopus", "squid axon" etc)
    P_Na = 1.0e-10 # Membrane permeability to sodium ions [m/s]
    P_K = 10.0e-10 # Membrane permeability to potassium ions [m/s]
    P_Cl = 5.0e-10 # Membrane permeability to chloride ions [m/s]

    Na_o = 145.0 # Extracellular sodium ion concentration [mol/m^3 or mmol/L]
    Na_i = 10.0 # Intracellular sodium ion concentration [mol/m^3 or mmol/L]
    K_o = 5.0 # Extracellular potassium ion concentration [mol/m^3 or mmol/L]
    K_i = 139.0 # Intracellular potassium ion concentration [mol/m^3 or mmol/L]
    Cl_o = 115.0 # Extracellular chloride ion concentration [mol/m^3 or mmol/L]
    Cl_i = 5.0 # Intracellular chloride ion concentration [mol/m^3 or mmol/L]

    # Metabolic parameters:
    delGo_ATP = -32e3  # Standard Gibbs Free energy ATP hydrolysis [J/mol]
    ATP = 4.0  # ATP concentration range 1.0 to 5.0 [mol/m^3 or mmol/L]
    ADP = 0.01  # ADP concentration range 0.004 to 0.1 [mol/m^3 or mmol/L]
    P = 0.01  # ADP concentration range 0.004 to 0.1  [mol/m^3 or mmol/L]
    Keqm_NaK = np.exp(delGo_ATP / (R * T)) # Equillibrium constant for ATP hydrolysis reaction

    omega_NaK = 1.0e-12 # Rate constant for the Na-K-ATPase ion pump
    # omega_nkcc = 0.0 # Rate constant for the NKCC ion symporter
    # omega_kcc = 0.0 # Rate constant for the KCC ion symporter



