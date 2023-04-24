#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright (c) 2022-2023 Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Calculate steady-state equilibrium concentrations and Vmem for the complete bioelectric system from
user-specified properties. This represents a top-level function in Calculion.
'''
from beartype import beartype, BeartypeConf
from beartype.typing import Optional
import pandas as pd
from pandas import DataFrame
from calculion.scratch_science.params import CalculionParams
from calculion.scratch_science.optimize import Optimizer
from calculion.scratch_science.time_solver import TimeSolver
from calculion.scratch_science.vmem import vrev_na, vrev_k, vrev_cl, vmem_ghk_pump
# import streamlit as st

# @st.cache
#@beartype(conf=BeartypeConf(is_debug=True))
@beartype
def get_steady_state(p: CalculionParams,
                     iterative_sol: bool=False,
                     update_env: bool=False,
                     quasi_static_vmem: bool=True,
                     ) -> tuple[DataFrame, DataFrame, Optional[TimeSolver]]:
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
    # Calculate an initial V_mem estimate with the starting parameters (in units of mV):
    V_mem_o = 1e3*vmem_ghk_pump(p.P_Na, p.Na_o, p.Na_i,
                            p.P_K, p.K_o, p.K_i,
                            p.P_Cl, p.Cl_o, p.Cl_i,
                            p.Keqm_NaK, p.omega_NaK,
                            p.ATP, p.ADP, p.P, p.alpha)

    # Calculate initial reversal potentials for the ions (in units of mV):
    V_rev_Na_o = 1e3*vrev_na(p.Na_o, p.Na_i, p.alpha)
    V_rev_K_o = 1e3*vrev_k(p.K_o, p.K_i, p.alpha)
    V_rev_Cl_o = 1e3*vrev_cl(p.Cl_o, p.Cl_i, p.alpha)

    # if update_env or not(quasi_static_vmem):
    #     p.delta_t = 0.1 # alter the time constant to handle the more sensitive computation of a small extracell space

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

        sim_time = None

    else:
        sim_time = TimeSolver(quasi_static_vmem=quasi_static_vmem,
                              update_env=update_env,
                              verbose=False)

        # Run the iterative solver until it reaches a steady-state or hits maximum number of iterations (p.N_iter):
        sim_time.time_loop(p, convergence_tol=p.steady_state_tol)

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

        # Calculate final reversal potentials for the ions in a time-dependent manner :
        sim_time.V_rev_Na_time = vrev_na(sim_time.Na_o_time, sim_time.Na_i_time, p.alpha)
        sim_time.V_rev_K_time = vrev_k(sim_time.K_o_time, sim_time.K_i_time, p.alpha)
        sim_time.V_rev_Cl_time = vrev_cl(sim_time.Cl_o_time, sim_time.Cl_i_time, p.alpha)


    # Calculate final reversal potentials for the ions (in units of mV):
    V_rev_Na_eq = 1e3*vrev_na(Na_o_eq, Na_i_eq, p.alpha)
    V_rev_K_eq = 1e3*vrev_k(K_o_eq, K_i_eq, p.alpha)
    V_rev_Cl_eq = 1e3*vrev_cl(Cl_o_eq, Cl_i_eq, p.alpha)

    round_dec = p.decimals_in_rounding

    # See also this Wikipedia article on Unicode superscripts and subscripts:
    #     https://en.wikipedia.org/wiki/Unicode_subscripts_and_superscripts
    steady_state_ions_dict = {
                              'Na⁺ᵢₙ': (round(p.Na_i, round_dec), round(Na_i_eq, round_dec)),
                              'K⁺ᵢₙ': (round(p.K_i, round_dec), round(K_i_eq, round_dec)),
                              'Cl⁻ᵢₙ': (round(p.Cl_i, round_dec), round(Cl_i_eq, round_dec)),
                              'Na⁺ₒᵤₜ': (round(p.Na_o, round_dec), round(Na_o_eq, round_dec)),
                              'K⁺ₒᵤₜ': (round(p.K_o, round_dec), round(K_o_eq, round_dec)),
                              'Cl⁻ₒᵤₜ': (round(p.Cl_o, round_dec), round(Cl_o_eq, round_dec))
                            }

    ions_dataframe = pd.DataFrame.from_dict(steady_state_ions_dict,
                                    orient='index', columns=['Initial [mM]', 'Final [mM]'])
    # ions_dataframe.style.set_caption('Ion Concentrations')

    steady_state_elec_dict = {'Vₘₑₘ': (round(V_mem_o, round_dec), round(V_mem_eq, round_dec)),
                                'Vᵣₑᵥ Na⁺': (round(V_rev_Na_o, round_dec), round(V_rev_Na_eq, round_dec)),
                                'Vᵣₑᵥ K⁺': (round(V_rev_K_o, round_dec), round(V_rev_K_eq, round_dec)),
                                'Vᵣₑᵥ Cl⁻': (round(V_rev_Cl_o, round_dec), round(V_rev_Cl_eq, round_dec)),
                                }

    elec_dataframe = pd.DataFrame.from_dict(steady_state_elec_dict,
                                            orient='index', columns=['Initial [mV]', 'Final [mV]'])
    # elec_dataframe.style.set_caption('Bioelectrical Properties')

    return ions_dataframe, elec_dataframe, sim_time