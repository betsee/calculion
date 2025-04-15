#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright (c) 2022-2025 Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
This module creates system for computing steady-state and iterative bioelectric simulations.

'''

import numpy as np
from scipy.optimize import minimize
from calculion.science.chem_opti import SteadyStateOpti
from calculion.science.reaction_system import ReactionSystem

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
