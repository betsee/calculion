#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright (c) 2021-2022 Ionovate.
# See "LICENSE" for further details.

'''
**Package metadata API** unit tests.

This submodule unit tests the public API of the
:mod:`calculion.meta` submodule.
'''

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To raise human-readable test errors, avoid importing from
# package-specific submodules at module scope.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# ....................{ TESTS                              }....................
def test_sim() -> None:
    '''
    Test the optimizer API of the :class:`calculion.science.optimize.Optimizer`
    class.
    '''

    # Defer test-specific imports.
    from calculion.science.optimize import Optimizer
    from calculion.science.params import CalculionParams

    p = CalculionParams()
    sim = Optimizer()

    # Test the optimizer's ability to find an equilibrium state:
    V_mem_eq, Na_i_eq, K_i_eq, Cl_i_eq = sim.find_eqm_state(p.P_Na, p.P_K, p.P_Cl,
                                                            p.Na_i, p.K_i, p.Cl_i,
                                                            p.Na_o, p.K_o, p.Cl_o,
                                                            p.ATP, p.ADP, p.P,
                                                            p.Keqm_NaK, p.omega_NaK,
                                                            p.alpha)


    # Test functionality that finds values of membrane permeability and pump rate that best match target vals:
    Vm_target = -10e-3
    P_Nat, P_Kt, P_Clt, omega_t = sim.find_target_values(Vm_target, p.P_Na, p.P_K, p.P_Cl, p.Na_i, p.K_i,
                                                         p.Cl_i, p.Na_o, p.K_o, p.Cl_o, p.ATP, p.ADP, p.P,
                                                         p.Keqm_NaK, p.omega_NaK,
                                                         p.alpha)

def test_tsim()->None:
    '''
    Test functionality associated with iterative simulations in time.
    '''

    from calculion.science.params import CalculionParams
    from calculion.science.time_solver import TimeSolver

    p = CalculionParams()
    p.delta_t = 1.0 # Adjust timestep and iterations to ensure that things run
    p.N_iter = 10

    sim_time = TimeSolver(quasi_static_vmem=True, update_env=False, verbose=True)

    sim_time.time_loop(p)
    sim_time.time_loop1(p)
    sim_time.time_loop2(p)

    sim_time = TimeSolver(quasi_static_vmem=False, update_env=False, verbose=True)
    sim_time.time_loop(p)
    sim_time.time_loop1(p)
    sim_time.time_loop2(p)

    sim_time = TimeSolver(quasi_static_vmem=True, update_env=True, verbose=True)
    sim_time.time_loop(p)
    sim_time.time_loop1(p)
    sim_time.time_loop2(p)

    sim_time = TimeSolver(quasi_static_vmem=False, update_env=True, verbose=True)
    sim_time.time_loop(p)
    sim_time.time_loop1(p)
    sim_time.time_loop2(p)

def test_calc() -> None:
    '''
    Test the top-level function that calculates bioelectrical steady-state parameters.

    '''
    from calculion.science.params import CalculionParams
    from calculion.science.compute import get_steady_state

    p = CalculionParams()  # Create a set of default parameters

    # Run through all the different ways to calculate the steady-state parameters:
    bioe_params = get_steady_state(p, iterative_sol=False, update_env=False, quasi_static_vmem=True)
    bioe_params = get_steady_state(p, iterative_sol=True, update_env=False, quasi_static_vmem=True)
    bioe_params = get_steady_state(p, iterative_sol=True, update_env=True, quasi_static_vmem=True)
    bioe_params = get_steady_state(p, iterative_sol=True, update_env=False, quasi_static_vmem=False)
    bioe_params = get_steady_state(p, iterative_sol=True, update_env=True, quasi_static_vmem=False)
