#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright (c) 2022-2023 Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Classes and methods to calculate properties for the complete bioelectric system,
including optimization methods to solve for equilibrium concentrations and
target values.
'''

# ....................{ IMPORTS                            }....................
# FIXME: Add in the itterative temporal simulator
# FIXME: Add in the symbolic math

import numpy as np
from numpy import ndarray
from beartype import beartype
from beartype.typing import Optional, Union
from calculion.science.flux import compute_flux_vector, make_reaction_matrix
from calculion.science.params import CalculionParams
from calculion.science.vmem import vmem_ghk_pump
from scipy.optimize import minimize


def calc_param_change(
        paramv: Union[list, ndarray],
        P_Na: float, P_K: float, P_Cl: float,
        Na_o: float, K_o: float, Cl_o: float,
        ATP: float, ADP: float, P: float,
        Keqm_NaK: float, omega_pump: float,
        alpha: float, M_react_sys: ndarray
):
    '''
    Function that is used with scipy optimization to solve for the
    equillibrium cell concentrations.
    '''

    Na_i = paramv[0]
    K_i = paramv[1]
    Cl_i = paramv[2]

    # Calculate Vmem using the J=0 quasi-steady state assumption:
    Vm = vmem_ghk_pump(P_Na, Na_o, Na_i,
                       P_K, K_o, K_i,
                       P_Cl, Cl_o, Cl_i,
                       Keqm_NaK, omega_pump,
                       ATP, ADP, P,
                       alpha)

    # Calculate the fluxes for each ion:
    flux_v = compute_flux_vector(P_Na, Na_o, Na_i,
                                 P_K, K_o, K_i,
                                 P_Cl, Cl_o, Cl_i,
                                 ATP, ADP, P,
                                 Keqm_NaK, omega_pump,
                                 Vm, alpha)

    # Calculate the change vector:
    change_v = M_react_sys.dot(flux_v)

    # function to minimize will be the square of the parameter change vector:
    loss_f = np.sum(change_v ** 2)

    return loss_f