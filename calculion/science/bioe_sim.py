#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright (c) 2022-2023 Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
This module creates system for computing steady-state and iterative bioelectric simulations.

'''

from beartype import beartype
import numpy as np
from scipy.optimize import minimize
from calculion.science.chem_opti import SteadyStateOpti
from calculion.science.reaction_system import ReactionSystem

@beartype
def solve_sys_steady_state(bes: ReactionSystem,
                            set_results: bool=False,
                            method: str = 'trust-constr'):
    '''
    Solve the steady-state optimization problem using Scipy's minimize.

    Parameters
    ----------
    bes : ReactionSystem
        An instance of ReactionSystem.

    set_results : bool
        If True, set the results of the optimization to the bes ReactionSystem.

    method : str
        The Scipy minimize optimization method to use.

    '''
    ss_params_list = []

    for chm in bes.chem_vect_reduced:
        ss_params_list.append(chm.symbol)

    ss_params_list.append(bes.V_mem_s)

    ss = SteadyStateOpti(bes, ss_params_list)

    bounds = [(0.0, 1.0e3) for i in ss_params_list]
    bounds[-1] = (-1.0, 1.0)

    sol0 = minimize(ss.divflux_opti_f,
                    ss.param_o,
                    args=(ss.args_vals),
                    jac=ss.jac_opti_f,
                    hess=ss.hess_opti_f,
                    method=method,
                    bounds=bounds
                    )

    if set_results:
        ss.set_results(bes, sol0.x)

    # Evaluate the Hessian at the solution:
    hess = ss.hess_opti_f(sol0.x, ss.args_vals)
    # Error on the parameters:
    xx_err = np.sqrt(np.diag(np.linalg.pinv(hess)))

    return ss, sol0.x, xx_err, sol0.fun

@beartype
def solve_sys_steady_state2(bes: ReactionSystem,
                           N_iter_max: int = 1000,
                           tol: float = 1.0e-20,
                           set_results: bool=False):
    '''
    Use the Gauss-Newton method to approximate the steady-state of the system.

    Parameters
    ----------
    bes : ReactionSystem
        An instance of ReactionSystem.

    N_iter_max : int
        Total maximum number of iterations the algorithm can run through.

    tol : float
        The tolerance, below which the solution is determined acheived.

    set_results : bool
        If True, set the results of the optimization to the bes ReactionSystem.

    '''
    ss_params_list = []

    for chm in bes.chem_vect_reduced:
        ss_params_list.append(chm.symbol)

    ss_params_list.append(bes.V_mem_s)

    ss = SteadyStateOpti(bes, ss_params_list)

    solss, sumsq = ss.gauss_newton(opti_N=N_iter_max, tol=tol)

    if set_results:
        ss.set_results(bes, solss)

    # Evaluate the Hessian at the solution:
    hess = ss.hess_opti_f(solss, ss.args_vals)
    # Error on the parameters:
    xx_err = np.sqrt(np.diag(np.linalg.pinv(hess)))

    return ss, solss, xx_err, sumsq

@beartype
def solve_pmem_steady_state(bes: ReactionSystem,
                            base_pmem: float=1.0,
                            set_results: bool=False):
    '''
    Get an estimate of the membrane permeabilities and pump/transporter rate constants
    needed to generate a particular profile of ion concentrations and V_mem.

    Parameters
    ----------
    bes : ReactionSystem
        An instance of ReactionSystem.

    set_results : bool
        If True, set the results of the optimization to the bes ReactionSystem.

    base_pmem : float
        The solution will be normalized to the first element of the optimized
        parameters list, which is typically the Na ion permeability. The base permeability
        is used to multiply the whole set of optimized membrane permeabilities to get them
        to useable units.
    '''

    ss_params_list = []

    for react in bes.react_vect:
        ss_params_list.append(react.P_s)

    ss = SteadyStateOpti(bes, ss_params_list)

    bounds = [(0.0, 1.0) for i in ss_params_list]

    sol0 = minimize(ss.divflux_opti_f,
                    ss.param_o,
                    args=(ss.args_vals),
                    jac=ss.jac_opti_f,
                    hess=ss.hess_opti_f,
                    method='trust-constr',
                    bounds=bounds
                    )

    xx = (sol0.x / sol0.x[0]) * base_pmem

    # Evaluate the Hessian at the solution:
    hess = ss.hess_opti_f(xx, ss.args_vals)
    # Error on the parameters:
    xx_err = np.sqrt(np.diag(np.linalg.pinv(hess)))

    if set_results:
        ss.set_results(bes, xx)

    return ss, xx, xx_err, sol0.fun