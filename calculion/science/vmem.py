#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright (c) 2022-2023 Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Methods to calculate transmembrane voltage V_mem.
'''

import numpy as np
from beartype import beartype

@beartype
def vmem_ghk(P_Na: float, Na_o: float, Na_i: float,
                P_K: float, K_o: float, K_i: float,
                P_Cl: float, Cl_o: float, Cl_i: float,
                alpha: float):
    '''
    Estimate V_mem by solving the system for zero current (J=0), which signifies the
    point at which V_mem is not changing in time and coincides with the GHK equation.

    This variation includes only ion concentrations and membrane permeabilities and excludes
    contributions from electrogenic pumps and transporters.

    This is the well-known Goldman-Hodgkin-Katz (GHK) equation.
    '''

    num = Cl_i* P_Cl + K_o*P_K + Na_o*P_Na
    den = Cl_o* P_Cl + K_i*P_K + Na_i*P_Na

    vmem = (1/alpha) * np.log(num / den)

    return vmem


@beartype
def vmem_ghk_pump(P_Na: float, Na_o: float, Na_i: float,
                  P_K: float, K_o: float, K_i: float,
                  P_Cl: float, Cl_o: float, Cl_i: float,
                  Keqm_NaK: float, omega_pump: float,
                  ATP: float, ADP: float, P: float,
                  alpha: float):
    '''
    Estimate V_mem by solving the system for zero current (J=0), which signifies the
    point at which V_mem is not changing in time and coincides with the GHK equation.

    This variation includes computation of the impact of the electrogenic NaK-ATPase ion
    pump.
    '''

    term_num = (2*ADP*Keqm_NaK*(K_i**2)*(Na_o**3)*omega_pump*P -
                3*ADP*Keqm_NaK*(K_i**2)*(Na_o**3)*omega_pump*P +
                Cl_i*P_Cl + K_o*P_K + Na_o*P_Na)

    term_den = (2*ATP*(K_o**2)*(Na_i**3)*omega_pump -
               3*ATP*(K_o**2)*(Na_i**3)*omega_pump +
               Cl_o*P_Cl + K_i*P_K + Na_i*P_Na)

    return (1/alpha)*np.log(term_num/term_den)

@beartype
def vmem_bioe(Vm: float, J: float, c_mem: float):
    '''
    Update V_mem given a transmembrane current J (cell-inward and electrically positive is positive)
    with a cell membrane patch capacitance c_mem.
    '''
    Vm += (1/c_mem)*J

    return Vm
