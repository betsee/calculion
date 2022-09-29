#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright (c) 2022-2023 Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Classes and methods to calculate properties for the complete single-cell bioelectric system,
including optimization methods to solve for equilibrium concentrations and
target values.
'''

# ....................{ IMPORTS                            }....................

import numpy as np
from numpy import ndarray
from beartype import beartype
from beartype.typing import Optional, Union
from calculion.science.flux import compute_flux_vector, make_reaction_matrix
from calculion.science.params import CalculionParams
from calculion.science.vmem import vmem_ghk_pump
from scipy.optimize import minimize

# ....................{ CLASSES                            }....................
@beartype
class Optimizer(object):
    '''
    **CalculIon single cell bioelectricity simulator** (i.e., high-level object
    performing core simulation computations).
    '''

    def __init__(self, p: Optional[CalculionParams] = None):
        '''
        Initialize this simulator
        '''

        # Create an instance of CalculionParams:
        if p is None:
            self.p = CalculionParams()
        else:
            self.p = p

        # Create a reaction matrix
        self.M_react_sys = make_reaction_matrix(self.p, quasi_static_vmem=True, update_env=False)


    def calc_param_change(
        self,
        paramv: Union[list, ndarray],
        P_Na: float, P_K: float, P_Cl: float,
        Na_o: float, K_o: float, Cl_o: float,
        ATP: float, ADP: float, P: float,
        Keqm_NaK: float, omega_pump: float,
        alpha: float
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
        change_v = self.M_react_sys.dot(flux_v)

        # function to minimize will be the square of the parameter change vector:
        loss_f = np.sum(change_v**2)

        return loss_f


    def calc_bioe_targets(
        self,
        paramv: Union[list, ndarray],
        Vm: float,
        Na_i: float, K_i: float, Cl_i: float,
        Na_o: float, K_o: float, Cl_o: float,
        ATP: float, ADP: float, P: float,
        Keqm_NaK: float, alpha: float
    ):
        '''
        Loss function used with scipy optimization to solve for the bioelectric
        sys params. Knowing target V_mem, ion concentrations, ATP
        concentrations, and temperature determine membrane permeabilities to
        best acheive the target state.
        '''

        P_Na = paramv[0]
        P_K = paramv[1]
        P_Cl = paramv[2]
        omega_pump = paramv[3]

        # Calculate the fluxes for each ion:
        flux_v = compute_flux_vector(P_Na, Na_o, Na_i,
                                     P_K, K_o, K_i,
                                     P_Cl, Cl_o, Cl_i,
                                     ATP, ADP, P,
                                     Keqm_NaK, omega_pump,
                                     Vm, alpha)

        # Calculate the change vector:
        change_v = self.M_react_sys.dot(flux_v)

        # function to minimize will be the square of the parameter changes:
        loss_f = np.sum(change_v**2)

        return loss_f


    def find_eqm_state(
        self,
        PNa: float, PK: float, PCl: float,
        Nai: float, Ki: float, Cli: float,
        Nao: float, Ko: float, Clo: float,
        ATP: float, ADP: float, P: float,
        Keqm_NaK: float, omega_pump: float,
        alpha: float,
    ):
        '''
        Given a set of membrane permeabilities and extracellular ion concentrations,
        compute the steady-state intracellular ion concentrations and V_mem for the
        system.

        '''
        self.paramv_o = np.asarray([Nai, Ki, Cli])  # Initial conditions

        sol0 = minimize(self.calc_param_change,
                        self.paramv_o,
                        args=(PNa, PK, PCl,
                              Nao, Ko, Clo,
                              ATP, ADP, P, Keqm_NaK, omega_pump,
                              alpha),
                        method='trust-constr',
                        jac=None,
                        hess=None,
                        hessp=None,
                        bounds=((0.0, 200.0), (0.0, 200.0), (0.0, 200.0)),
                        tol=None,
                        callback=None,
                        options=None)

        # Calculate Vmem at these conditions:
        self.V_mem = vmem_ghk_pump(PNa, Nao, sol0.x[0],
                              PK, Ko, sol0.x[1],
                              PCl, Clo, sol0.x[2],
                              Keqm_NaK, omega_pump,
                              ATP, ADP, P,
                              alpha)

        self.Na_i = sol0.x[0]
        self.K_i = sol0.x[1]
        self.Cl_i = sol0.x[2]

        self.sol0 = sol0

        return self.V_mem, self.Na_i, self.K_i, self.Cl_i


    def find_target_values(
        self,
        Vm: float,
        PNa: float, PK: float, PCl: float,
        Nai: float, Ki: float, Cli: float,
        Nao: float, Ko: float, Clo: float,
        ATP: float, ADP: float, P: float,
        Keqm_NaK: float, omega_pump: float,
        alpha: float,
    ):
        '''
        Given a set of membrane permeabilities and extracellular ion
        concentrations, compute the steady-state intracellular ion
        concentrations and V_mem for the system.
        '''

        self.paramv_o = np.asarray([PNa, PK, PCl, omega_pump])  # Initial conditions

        self.sol0 = minimize(self.calc_bioe_targets,
                        self.paramv_o,
                        args=(Vm, Nai, Ki, Cli,
                              Nao, Ko, Clo,
                              ATP, ADP, P,
                              Keqm_NaK, alpha),
                        method='trust-constr',
                        jac=None,
                        hess=None,
                        hessp=None,
                        bounds=((0.0, 1.0), (0.0, 1.0), (0.0, 1.0), (0.0, 1.0)),
                        tol=None,
                        callback=None,
                        options=None)

        self.P_Na = self.sol0.x[0]
        self.P_K = self.sol0.x[1]
        self.P_Cl = self.sol0.x[2]
        self.omega_pump = self.sol0.x[3]

        return self.P_Na, self.P_K, self.P_Cl, self.omega_pump
