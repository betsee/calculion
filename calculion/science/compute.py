#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright (c) 2022-2023 Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Calculate steady-state equilibrium concentrations and Vmem for the complete bioelectric system from
user-specified properties. This represents a top-level function in Calculion.
'''
from beartype import beartype
from calculion.science.params import CalculionParams
from calculion.science.optimize import Optimizer
from calculion.science.time_solver import TimeSolver

@beartype
def get_steady_state(p: CalculionParams,
                     iterative_sol: bool=False,
                     update_env: bool=False,
                     quasi_static_vmem: bool=True,
                     ) -> dict:
    '''
    Function to compute the steady-state values of the bioelectrical system (i.e. the equilibrium ion
    concentrations and Vmem).

    Parameters
    -----------
    p  :  CalculionParams
        All parameters for the simulation, which may take alterations from the user.

    iterative_sol : bool
        Use the steady-state solver based on Scipy's Gauss-Newton solver method (False) or perform an iterative
        time simulation until it reaches steady-state (True). Iterative simulation is more accurate. If
        iterative_sol = False, update_env and quasi_static_vmem arguments are ignored as these are frozen as
        update_env=False and quasi_static_vmem = True in the non-iterative solver method.

    update_env : bool
        Simulate as if the cell is surrounded by an enormous bath of media (i.e. env concentrations will not change,
        update_env = False), or, use a true extracellular volume and update concentrations in that volume
        (update_env = True).This argument is ignored when iterative_sol = False.

    quasi_static_vmem : bool
        Use the quasi-static Vmem approximation, which assumes the electrical component reaches steady state many
        orders of magnitude faster than the chemical system, so that we can assume J=0 at the membrane at each
        iterative timestep. This approximation greatly improves stability, allowing for use of large time-steps.
        Alternatively, if quasi-static_vmem = False, then the membrane is assumed to be a capacitor and the non-zero
        transmembrane current is used to update Vmem at each time-step.

    '''
    if update_env or not(quasi_static_vmem):
        p.delta_t = 0.1 # alter the time constant to handle the more sensitive computation of a small extracell space

    if iterative_sol is False:
        sim = Optimizer(p) # Create an instance of the Calculion Optimizer

        # Using values from the params object, determine the steady state of the system, assuming the
        # cell is surrounded by a very large bath of extracellular media (i.e. no changes to extracellular
        # values):
        V_mem_eqo, Na_i_eq, K_i_eq, Cl_i_eq = sim.find_eqm_state(p.P_Na, p.P_K, p.P_Cl,
                                                                p.Na_i, p.K_i, p.Cl_i,
                                                                p.Na_o, p.K_o, p.Cl_o,
                                                                p.ATP, p.ADP, p.P,
                                                                p.Keqm_NaK, p.omega_NaK,
                                                                p.alpha)

        V_mem_eq = 1.0e3*V_mem_eqo

        # This steady-state solver method always assumes the extracellular concentrations are fixed:
        Na_o_eq = p.Na_o
        K_o_eq = p.K_o
        Cl_o_eq = p.Cl_o

    else:
        sim_time = TimeSolver(quasi_static_vmem=quasi_static_vmem,
                              update_env=update_env,
                              verbose=False)

        # Run the iterative solver until it reaches a steady-state or hits maximum number of iterations (p.N_iter):
        sim_time.time_loop(p)

        V_mem_eq = 1e3*sim_time.V_mem_time[-1]
        Na_i_eq = sim_time.Na_i_time[-1]
        K_i_eq = sim_time.K_i_time[-1]
        Cl_i_eq = sim_time.Cl_i_time[-1]

        if update_env:
            Na_o_eq = sim_time.Na_o_time[-1]
            K_o_eq = sim_time.K_o_time[-1]
            Cl_o_eq = sim_time.Cl_o_time[-1]

        else:
            Na_o_eq = p.Na_o
            K_o_eq = p.K_o
            Cl_o_eq = p.Cl_o

    steady_state_params_dict = {'Na_i': Na_i_eq,
                                'Na_o': Na_o_eq,
                                'K_i': K_i_eq,
                                'K_o': K_o_eq,
                                'Cl_i': Cl_i_eq,
                                'Cl_o': Cl_o_eq,
                                'V_mem': V_mem_eq}

    return steady_state_params_dict



