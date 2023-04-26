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
    Test the optimizer API of the :class:`calculion.scratch_science.optimize.Optimizer`
    class.
    '''

    # Defer test-specific imports.
    from calculion.science.model_params import ModelParams
    from calculion.science.sim_params import SimParams
    from calculion.science.bioe_system import BioElectricSystem
    from calculion.science.bioe_sim import solve_sys_steady_state
    from calculion.science.channel_base import PulseFunctionChannel
    from calculion.science.chem_opti import IterSim

    p = ModelParams()  # Create a default parameters instance for model properties
    sim_p = SimParams() # Create default params for simulation properties

    # Turn on all pumps and transporters for maximal test coverage:
    p.use_NaK_ATPase = True
    p.use_NaKCl = True
    p.use_KCl = True

    p.update_parameters() # update all parameters

    sim = BioElectricSystem(p)  # Create the full bioelectrical study object
    besi = sim.bes  # alias to the ReactionSystem object

    # First calculate the steady-state of the system:
    besi.V_mem = 0.0

    # steady-state solver, solution, param err, total sum squares error:
    ss0, results_vect, x_err0, err0 = solve_sys_steady_state(besi,
                                                             method='trust-constr',
                                                             set_results=True)

    # Test the iterative solver:

    stepchan_Na = PulseFunctionChannel(['P_Na'],
                                       sim_p.perturb_PNa_multi * p.base_pmem,
                                       sim_p.perturb_PNa_start,
                                       sim_p.perturb_PNa_end)

    stepchan_K = PulseFunctionChannel(['P_K'],
                                       sim_p.perturb_PK_multi * p.base_pmem,
                                       sim_p.perturb_PK_start,
                                       sim_p.perturb_PK_end)

    stepchan_Cl = PulseFunctionChannel(['P_Cl'],
                                       sim_p.perturb_PCl_multi * p.base_pmem,
                                       sim_p.perturb_PCl_start,
                                       sim_p.perturb_PCl_end)

    chanlist = [stepchan_Na, stepchan_K, stepchan_Cl]

    isim = IterSim(besi, channels_list=chanlist)

    time, vm_time, chem_time = isim.run_iter_sim(sim_p.delta_t,
                                                 sim_p.N_iter,
                                                 use_quasi_static_approx=False,
                                                 use_hodgkin_huxley=False,
                                                 clamp_vmem_at=None,
                                                 sweep_vmem_vals=None)

    time, vm_time, chem_time = isim.run_iter_sim(sim_p.delta_t,
                                                 sim_p.N_iter,
                                                 use_quasi_static_approx=False,
                                                 use_hodgkin_huxley=True,
                                                 clamp_vmem_at=None,
                                                 sweep_vmem_vals=None)

    time, vm_time, chem_time = isim.run_iter_sim(sim_p.delta_t,
                                                 sim_p.N_iter,
                                                 use_quasi_static_approx=False,
                                                 use_hodgkin_huxley=False,
                                                 clamp_vmem_at=-0.05,
                                                 sweep_vmem_vals=None)

    time, vm_time, chem_time = isim.run_iter_sim(sim_p.delta_t,
                                                 sim_p.N_iter,
                                                 use_quasi_static_approx=False,
                                                 use_hodgkin_huxley=False,
                                                 clamp_vmem_at=None,
                                                 sweep_vmem_vals=(-0.08, 0.020))





