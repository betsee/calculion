#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright (c) 2022-2023 Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
**Core application entry-point** (i.e., submodule defining the root
:func:`main` function running this app, intended to be imported from elsewhere
within this codebase at both runtime and test-time).

Specifically, this submodule is imported by:

* The top-level :mod:`calculion.__main__` submodule, implicitly run by the
  active Python interpreter when passed the ``--m`` option on the command line
  (e.g., ``python3 -m calculion``).
* Integration tests programmatically exercising app functionality.
'''
import copy

import pandas as pd


# ....................{ KLUDGES ~ path                     }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# CAUTION: Kludge PYTHONPATH *BEFORE* importing from this package below.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Explicitly register all files and subdirectories of the parent directory
# containing this module to be importable modules and packages (respectively)
# for the remainder of this Python process if this directory has yet to be
# registered.
#
# Technically, this should *NOT* be required. Streamlit should implicitly
# guarantee this to be the case. Indeed, when Streamlit is run as a module
# (e.g., as "python3 -m streamlit run {app_name}/main.py"), Streamlit does just
# that. Unfortunately, when Streamlit is run as an external command (e.g., as
# "streamlit run {app_name}/main.py"), Streamlit does *NOT* guarantee this to be
# the case. Since Streamlit Cloud runs Streamlit as an external command rather
# than as a module, Streamlit Cloud effectively does *NOT* guarantee this to be
# the case as well.

# Isolate this kludge to a private function for safety.
def _register_dir() -> None:

    # Defer kludge-specific imports. Avert thy eyes, purist Pythonistas!
    from logging import info
    from pathlib import Path
    from sys import path as sys_path

    # Log this detection attempt.
    info('[APP] Detecting whether app package directory requires registration on "sys.path": %s', sys_path)
    # print('Registering app package directory for importation: %s')

    # Path object encapsulating the absolute filename of the file defining the
    # current module. Note that doing so may raise either:
    # * If this file inexplicably does *NOT* exist, "FileNotFoundError".
    # * If this file inexplicably resides under a directory subject to an
    #   infinite symbolic link loop, "RuntimeError".
    main_file = Path(__file__).resolve(strict=True)

    # Absolute dirname of the parent directory containing this app's top-level
    # package, which is guaranteed to be either:
    # * If this app is currently installed editably (e.g., "pip install -e ."),
    #   the repository directory containing the ".git/" directory for this app.
    # * If this app is currently installed non-editably (e.g., "pip install ."),
    #   the equivalent of the "site-packages/" directory for the active Python.
    package_dirname = str(main_file.parents[1])

    # If the current PYTHONPATH does *NOT* already contain this directory...
    if package_dirname not in sys_path:
        # Log this registration attempt.
        info('[APP] Registering app package directory for importation: %s', package_dirname)
        # print('Registering app package directory for importation: %s')

        # Append this directory to the current PYTHONPATH.
        sys_path.append(package_dirname)

# Kludge us up the bomb.
_register_dir()

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# CAUTION: Avoid importing anything at module scope *EXCEPT* from official
# Python modules in the standard library guaranteed to exist. Subsequent logic
# in the _main() function called below validates third-party runtime
# dependencies of this package to be safely importable, Before performing that
# validation, *NO* other modules are safely importable from.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# ....................{ MAIN                               }....................
def main() -> None:
    '''
    Core function running this Streamlit-based web app: **Calculion.**
    '''

    # ..................{ IMPORTS                            }..................
    import streamlit as st
    from PIL import Image
    import numpy as np
    from calculion.science.model_params import ModelParams
    from calculion.science.sim_params import SimParams
    from calculion.science.bioe_system import BioElectricSystem
    from calculion.science.bioe_sim import solve_sys_steady_state
    from calculion.science.channel_base import PulseFunctionChannel
    from calculion.science.chem_opti import IterSim
    from calculion.science.string_names import StringNames
    from calculion._util.path.utilpathself import (
        get_data_png_cell_network_schematic_0_file,
        get_data_png_cell_network_schematic_1_file,
        get_data_png_cell_network_schematic_2_file,
        get_data_png_cell_network_schematic_3_file,
        get_data_png_cell_network_schematic_4_file,
        get_data_png_cell_network_schematic_5_file,
        get_data_png_cell_network_schematic_6_file,
        get_data_png_cell_network_schematic_7_file,
        get_data_png_membrane_schematic_file,
        get_data_png_banner_file
    )
    # from pandas import DataFrame
    # from calculion.scratch_science.compute import get_steady_state
    from streamlit import (
        set_page_config
    )

    # ..................{ HEADERS                            }..................
    set_page_config(layout="wide") # set a wide page configuration?

    # Human-readable title of this web app.
    # title('CalculIon')
    banner_image_fn = str(get_data_png_banner_file())
    banner_image = Image.open(banner_image_fn)
    st.image(banner_image,
             use_column_width='always',
             output_format="PNG")

    # App subtitle, if we want it:
    # st.write('Calculating the *slow* changes of bioelectricity')

    # ..................{ LOCALS                             }..................
    p = ModelParams()  # Create a default parameters instance for model properties
    sim_p = SimParams() # Create default params for simulation properties
    l = StringNames() # string labels

    # ..................{ SIDEBAR                            }..................
    # The sidebar will contain all widgets to collect user-data for the
    # simulation. Create and name the sidebar.
    st.sidebar.header('Model Variables')

    # st.sidebar.write('#### Set simulation variables')
    with st.sidebar:

        st.write('**Initial Conditions**')

        # slider general settings:
        ion_min_val = 0.1  # min concentration of ions in mM
        ion_max_val = 250.0  # max concentration of ions in mM
        ion_slider_step = 0.1  # step-size for ion slider
        memp_min_val = 0.0 # min value for membrane permeability in nm/s
        memp_max_val = 200.0 # max value for membrane permeability in nm/s
        memp_slider_step = 0.01 # step-size for mem perm slider

        # Form of the widget to use for parameter value input from the user:
        # param_widget = st.slider # form of the widget to use for parameters
        param_widget = st.number_input # form of the widget to use for parameters

        # Initialize an expander block for each set of variables:
        default_expanded_state = False # default state of the expander is expanded?

        # Define another expander block collecting parameters for intracellular ion concentrations:
        ion_in_block = st.expander("Ion Concentrations Inside Cell",
                                   expanded=default_expanded_state,
                                   )

        with ion_in_block:

            p.cNa_i = param_widget('Na+ in [mM]',
                                  min_value=ion_min_val,
                                  max_value=ion_max_val,
                                  value=p.cNa_i,
                                  step=ion_slider_step,
                                  format='%f',
                                  key='slider_Na_i',
                                  label_visibility='visible',
                                  help='Set the starting value of Na⁺ in the cell.')

            p.cK_i = param_widget('K+ in [mM]',
                                 min_value=ion_min_val,
                                 max_value=ion_max_val,
                                 value=p.cK_i,
                                 step=ion_slider_step,
                                 format='%f',
                                 key='slider_K_i',
                                 label_visibility='visible',
                                 help='Set the starting value of K⁺ in the cell.')

            p.cCl_i = param_widget('Cl- in [mM]',
                                  min_value=ion_min_val,
                                  max_value=ion_max_val,
                                  value=p.cCl_i,
                                  step=ion_slider_step,
                                  format='%f',
                                  key='slider_Cl_i',
                                  label_visibility='visible',
                                  help='Set the starting value of Cl⁻ in the cell.')

        st.write('**Model Parameters**')

        # Define another expander block for specifying extracellular ion concentrations:
        ion_out_block = st.expander("Ion Concentrations Outside Cell", expanded=default_expanded_state)

        with ion_out_block:
            # Automatically reset parameter values in the parameters object p depending on user-selection:
            p.cNa_o = param_widget('Na+ out [mM]',
                             min_value=ion_min_val,
                             max_value=ion_max_val,
                             value=p.cNa_o,
                             step=ion_slider_step,
                             format='%f',
                             key='slider_Na_o',
                             label_visibility='visible',
                             help='Set the Na⁺ concentration outside the cell.')

            p.cK_o = param_widget('K+ out [mM]',
                             min_value=ion_min_val,
                             max_value=ion_max_val,
                             value=p.cK_o,
                             step=ion_slider_step,
                             format='%f',
                             key='slider_K_o',
                             label_visibility='visible',
                             help='Set the K⁺ concentration outside the cell.')

            p.cCl_o = param_widget('Cl- out [mM]',
                             min_value=ion_min_val,
                             max_value=ion_max_val,
                             value=p.cCl_o,
                             step=ion_slider_step,
                             format='%f',
                             key='slider_Cl_o',
                             label_visibility='visible',
                             help = 'Set the Cl⁻ concentration outside the cell.')

        memperm_block = st.expander("Cell Membrane Ion Permeabilities", expanded=default_expanded_state)

        with memperm_block:

            p.base_PNa = param_widget('Na+ Permeability [nm/s]:',
                               min_value=memp_min_val,
                               max_value=memp_max_val,
                               value=p.base_PNa,
                               step=memp_slider_step,
                               format='%f',
                               key='slider_P_Na',
                               label_visibility='visible',
                               help='Set the base membrane permeability to Na⁺.\n'
                                    '\nThis simulates cellular expression of Na⁺ leak channels.'
                               )

            p.base_PK = param_widget('K+ Permeability [nm/s]:',
                               min_value=memp_min_val,
                               max_value=memp_max_val,
                               value=p.base_PK,
                               step=memp_slider_step,
                               format='%f',
                               key='slider_P_K',
                               label_visibility='visible',
                               help='Set the base membrane permeability to K⁺.\n'
                                    '\nThis simulates cellular expression of K⁺ leak channels.'
                               )

            p.base_PCl = param_widget('Cl- Permeability [nm/s]:',
                               min_value=memp_min_val,
                               max_value=memp_max_val,
                               value=p.base_PCl,
                               step=memp_slider_step,
                               format='%f',
                               key='slider_P_Cl',
                               label_visibility='visible',
                               help='Set the base membrane permeability to Cl⁻.\n'
                                    '\nThis simulates cellular expression of Cl⁻ leak channels.'
                               )

            # # Update the main membrane permeabilities used in the simulation to m/s units:
            # p.PNa = p.base_pmem*p.base_PNa
            # p.PK = p.base_pmem*p.base_PK
            # p.PCl = p.base_pmem*p.base_PCl

        # Define another expander block for pump and transporter settings:
        pumps_block = st.expander("Ion Pump and Transporters", expanded=default_expanded_state)

        with pumps_block:

            # # Let the user specify whether they want the additional NaKCl cotransporter for secondary transport:
            NaKATP_on = st.checkbox(l.NaK_pump, value=True, help=f'Cell expresses {l.NaK_pump} ?')

            if NaKATP_on:

                p.use_NaK_ATPase = True

                p.base_NaKpump = param_widget('Na-K-ATPase Pump Rate [units]',
                                 min_value=0.0,
                                 max_value=1.0e6,
                                 value=p.base_NaKpump,
                                 step=0.01,
                                 format='%f',
                                 key='slider_omega_NaK',
                                 label_visibility='visible',
                                 help='Set the maximum rate of the Na-K-ATPase ion pump.')

                # update the na-k-atpase pump rate in units used in the simulation:
                # p.PNaK_ATPase = p.base_pmem*p.base_NaKpump*p.pump_unit_modifier
            else:
                p.use_NaK_ATPase = False
                p.base_NaKpump = 0.0 # otherwise set the rate to zero

            # # Let the user specify whether they want the additional NaKCl cotransporter for secondary transport:
            NaKCl_on = st.checkbox(l.NaKCl_cotrans_o,
                                   value=False,
                                   help=f'Cell expresses {l.NaKCl_cotrans} ?')

            if NaKCl_on:
                # Na-K-2Cl cotransporter properties:
                p.use_NaKCl = True

                p.base_NaKCl = param_widget('Na-K-2Cl Cotransporter Rate [units]',
                                 min_value=0.0,
                                 max_value=1.0e6,
                                 value=p.base_NaKCl,
                                 step=0.01,
                                 format='%f',
                                 key='slider_omega_NaKCl',
                                 label_visibility='visible',
                                 help='Set the maximum rate of the Na-K-2Cl cotransporter.')

                # update the na-k-2Cl cotransporter rate in units used in the simulation:
                # p.PNaKCl = p.base_pmem*p.base_NaKCl*p.pump_unit_modifier

            else:
                p.use_NaKCl = False
                p.base_NaKCl = 0.0 # otherwise set the rate to zero

            # # Let the user specify whether they want the additional KCl symporter for secondary transport:
            KCl_on = st.checkbox(l.KCl_symp_o,
                                 value=False,
                                 help=f'Cell expresses {l.KCl_symp} ?')

            if KCl_on:
                p.use_KCl = True
                # p.base_KCl = 0.0 # otherwise set the rate to zero
                # K-Cl symporter properties:
                p.base_KCl = param_widget('K-Cl Symporter Rate [units]',
                                 min_value=0.0,
                                 max_value=1.0e8,
                                 value=p.base_KCl,
                                 step=0.1,
                                 format='%f',
                                 key='slider_omega_KCl',
                                 label_visibility='visible',
                                 help='Set the maximum rate of the K-Cl symporter.'
                                           )

                # update the K-Cl symporter rate in units used in the simulation:
                # p.PKCl = p.base_pmem*p.base_KCl*p.pump_unit_modifier

            else:
                p.use_KCl = False
                p.base_KCl = 0.0 # otherwise set the rate to zero

        # Define another expander block for pump and transporter settings:
        metabolic_block = st.expander("Metabolic Settings", expanded=default_expanded_state)

        with metabolic_block:

            delGATP = param_widget('Gibbs Free Energy ATP Hydrolysis [kJ/mol]',
                             min_value=-34.0,
                             max_value=-28.0,
                             value=p.delG_ATP*1e-3,
                             step=0.1,
                             format='%f',
                             key='slider_delGATP',
                             label_visibility='visible',
                             help='Set the Gibbs standard free energy for ATP hydrolysis.')

            p.delG_ATP = delGATP*1e3 # convert to units J/mol from kJ/mol

            # Concentrations of metabolic items (these are not updated)
            p.cATP = param_widget('ATP Concentration [mM]',
                             min_value=0.001,
                             max_value=5.0,
                             value=p.cATP,
                             step=0.001,
                             format='%f',
                             key='slider_ATP',
                             label_visibility='visible',
                             help='Set the concentration of ATP in the cell.')

            p.cADP = param_widget('ADP Concentration [mM]',
                             min_value=0.001,
                             max_value=5.0,
                             value=p.cADP,
                             step=0.001,
                             format='%f',
                             key='slider_ADP',
                             label_visibility='visible',
                             help='Set the concentration of ADP in the cell.')

            p.cPi = param_widget('P Concentration [mM]',
                             min_value=0.001,
                             max_value=5.0,
                             value=p.cPi,
                             step=0.001,
                             format='%f',
                             key='slider_P',
                             label_visibility='visible',
                             help='Set the concentration of phosphorous (P) in the cell.')

        # Define an expander block for cell and env settings:
        cellenv_settings_block = st.expander("Cell and Environmental Settings",
                                             expanded=default_expanded_state)
        with cellenv_settings_block:
            p.T_C = param_widget('Temperature [degrees C]',
                             min_value=1.0,
                             max_value=60.0,
                             value=p.T_C,
                             step=0.1,
                             format='%f',
                             key='slider_T',
                             label_visibility='visible',
                             help='Set the temperature of the system.')

            # Update temperature in Kelvin used in simulations:
            p.T = p.T_C + 273.15 # temp in Kelvin

            p.r_cell_um = param_widget('Cell radius [um]',
                             min_value=1.0,
                             max_value=25.0,
                             value=p.r_cell_um,
                             step=0.1,
                             format='%f',
                             key='slider_rcell',
                             label_visibility='visible',
                             help='Set the radius of the cell.')

            # update r_cell to be in meters for simulations:
            p.r_cell = p.r_cell_um*1e-6


    # ..................{ RESULTS                            }..................
    # After collecting parameter values from the user, compute the steady-state
    # values for the bioelectrical system.
    # @st.cache_data
    def make_bioe_system(p):
        sim = BioElectricSystem(p)  # Create the full bioelectrical study object
        besi = sim.bes  # alias to the ReactionSystem object
        return besi

    # @st.cache_resource # Cache the results of this slower function
    def calculate_ss_results(_besi):

        # sim = BioElectricSystem(p)  # Create the full bioelectrical study object
        # bes = sim.bes  # alias to the ReactionSystem object

        # First calculate the steady-state of the system:
        _besi.V_mem = 0.0

        # steady-state solver, solution, param err, total sum squares error:
        ss0, results_vect, x_err0, err0 = solve_sys_steady_state(_besi,
                                                         method='trust-constr',
                                                         set_results=True)
        # Get the concentration ss dataframe:
        ion_ss = _besi.return_chem_props_dict()

        # Get the electrical properties dataframe:
        elec_ss = _besi.return_elec_props_dict()

        return ion_ss, elec_ss

    # @st.cache_data # Cache the results of this slower function
    def calculate_iter_results(_besi, simp, chanlist):

        isim = IterSim(_besi, channels_list=chanlist)

        time, vm_time, chem_time = isim.run_iter_sim(simp.delta_t,
                                                     simp.N_iter,
                                                     use_quasi_static_approx=False,
                                                     use_hodgkin_huxley=False,
                                                     clamp_vmem_at=None,
                                                     sweep_vmem_vals=None)


        return time, vm_time, chem_time, isim




    # ..................{ TABS                               }..................
    # Split the main area into these three tabs:
    # * The "Introduction" tab (tab1) will display a write-up of the theory
    #   behind Calculion.
    # * The "Simulation Results" tab will display the simulation results as
    #   tables and charts.
    # * The "Bioelectrical Network" tab will show a graphical depiction of the
    #   bioelectrical network.
    tab1, tab2, tab3 = st.tabs([
        'Introduction', 'Simulation', 'Bioelectrical Network'])

    with tab1:
        st.write('### Calculating the *slow* changes of bioelectricity.')
        # # This is how to serve up some html:
        # test_html_fn = get_data_html_test_file()
        # with open(test_html_fn,'r') as f:
        #     html_data = f.read()

        # st.components.v1.html(html_data, scrolling=True)
        mem_image_fn = str(get_data_png_membrane_schematic_file())
        mem_image = Image.open(mem_image_fn)
        st.image(mem_image,
                 # caption='Behold the Cellular Bioelectric Network!',
                 use_column_width='always',
                 output_format="PNG")

        st.write(f"Transmembrane electric potential, {l.Vmem}, plays a central "
                 f"role in numerous biological processes with a well-demonstrated "
                 f"capability to influence cell physiology, and potential "
                 f"applications in regenerative medicine and bioengineering. "
                 )

        st.write(f"Well-known models of {l.Vmem} and its dynamics "
                 f"(e.g. the Goldman-Hodgkin-Katz Voltage equation for {l.Vmem} and "
                 f"the Hodgkin-Huxley description of dynamic changes to {l.Vmem} "
                 f"with ion channel activity) typically neglect changes to ion "
                 f"concentrations inside and out of the cell, assuming they are "
                 f"constant with time."
                 )

        st.write(f"However, these slowly changing ion "
                 f"concentration gradients, which are primarily generated by "
                 f"the activity of ion pumps and transporters, are key specifiers "
                 f"of the steady-state value of {l.Vmem} (also known as the "
                 f"'resting potential') and to {l.Vmem} dynamics with ion channel "
                 f"activity.")

        st.write(f"We developed "
                 f"a comprehensive bioelectrical model for a single cell that "
                 f"examines changes to transmembrane ion concentrations and {l.Vmem}. ")

        st.write(f"This model is the engine behind this CalculIon web app.")

        st.write(f"The **sidebar** to the left allows you to specify a variety of "
                 f"system parameters to build a bioelectrical system with different "
                 f"leak channel, ion pump, and transporter expression levels."
                 )

        st.write(f"Use the **Simulation** tab to explore the outcomes of different "
                 f"system properties from both a steady-state and dynamic perspectives."
                 )

        st.write(f"The **Bioelectric Network** tab shows a schematic directed graph of "
                 f"the bioelectric system that is being modeled, with deeper "
                 f"explanations for the observed simulation results.")

        st.write(f"The details of our model emphasizing the *slow* changes "
                 f"that shape the electrochemical ion gradients fundamental "
                 f"to bioelectricity can be downloaded *here*.")

    with tab2:
        st.write("### Simulation")

        st.write("*Alter Model Variables in sidebar to explore outcomes...*")

        # Define a final expander block for simulator settings:
        sim_settings_block = st.expander("Simulation Settings",
                                         expanded=default_expanded_state)

        with sim_settings_block:

            # Iterative solver will not be used by default:
            itersol_checkbox = st.checkbox("Use iterative solver",
                                           value=sim_p.use_iterative_solver,
                                           key='checkbox_itersol',
                                           help='Use the iterative solver that integrates the '
                                                'system step-by-step in time?')

            if itersol_checkbox:
                sim_p.use_iterative_solver = True # Set the iterative solver parameter to True

                # Iterative solver time step:
                sim_p.delta_t = param_widget('Simulation time-step [s]',
                                          min_value=1.0e-4,
                                          max_value=1.0,
                                          value=sim_p.delta_t,
                                          step=1.0e-4,
                                          format='%f',
                                          key='slider_delta_t',
                                          label_visibility='visible',
                                          help='Set the time step for the iterative solver.')

                # Iterative solver max iterations:
                sim_p.start_time = param_widget('Start time for the iterative simulation.',
                                          min_value=0.0,
                                          max_value=1.0e5*sim_p.delta_t,
                                          value=sim_p.start_time,
                                          step=sim_p.delta_t,
                                          format='%f',
                                          key='slider_Niter',
                                          label_visibility='visible',
                                          help='Set the maximum number of timesteps that can be run.')

                # Iterative solver end time:
                sim_p.end_time = param_widget('End time',
                                          min_value=sim_p.start_time + sim_p.delta_t,
                                          max_value=sim_p.start_time + 1e6*sim_p.delta_t,
                                          value=sim_p.end_time,
                                          step=sim_p.delta_t,
                                          format='%f',
                                          key='slider_endt',
                                          label_visibility='visible',
                                          help='Set the simulation time at '
                                               'which the iterative solver stops.')


                # membrane permeability perturbations for sim:
                # FIXME: we need panels to harvest user pref on channel sim here:
                chan_list = []
                if sim_p.perturb_PNa:
                    stepchan_Na = PulseFunctionChannel(['P_Na'],
                                                       sim_p.perturb_PNa_multi*p.base_pmem,
                                                       sim_p.perturb_PNa_start,
                                                       sim_p.perturb_PNa_end)
                    chan_list.append(stepchan_Na)

                if sim_p.perturb_PK:
                    stepchan_K = PulseFunctionChannel(['P_K'],
                                                      sim_p.perturb_PK_multi*p.base_pmem,
                                                      sim_p.perturb_PK_start,
                                                      sim_p.perturb_PK_end)
                    chan_list.append(stepchan_K)

                if sim_p.perturb_PCl:
                    stepchan_Cl = PulseFunctionChannel(['P_Cl'],
                                                       sim_p.perturb_PCl_multi*p.base_pmem,
                                                       sim_p.perturb_PCl_start,
                                                       sim_p.perturb_PCl_end)
                    chan_list.append(stepchan_Cl)


                # Update all the calculated parameters with the new settings:
                # p.update_parameters()

                # isim = IterSim(bes, channels_list=[stepchan_Na, stepchan_K, stepchan_Cl])
                # isim2 = IterSim(bes, channels_list=[stepchan_Na, stepchan_K, stepchan_Cl])

        st.write("### Results")

        p.update_parameters() # Update parameters in case calculation forgotten
        sim_p.update_parameters() # Update parameters in case calculation forgotten

        # Generate the bioelectrical system model object:
        # sim = BioElectricSystem(p)  # Create the full bioelectrical study object
        bes = make_bioe_system(p)
        # bes = sim.bes  # alias to the ReactionSystem object
        beso = copy.deepcopy(bes)  # Make a copy to save initial state

        # Calculate the initial conditions of the system:
        # Get the concentration ss dataframe:
        ion_vals_init = beso.return_chem_props_dict()

        # Get the electrical properties dataframe:
        elec_vals_init = beso.return_elec_props_dict()

        # Calculate the steady-state results of the system:
        ion_vals_ss, elec_vals_ss = calculate_ss_results(bes)

        col1, col2 = st.columns(2)

        with col1:
            st.write('###### Initial Bioelectrical Potentials')
            st.dataframe(elec_vals_init.style.format("{:.1f}"))

            st.write('###### Steady-State Bioelectrical Potentials')
            st.dataframe(elec_vals_ss.style.format("{:.1f}"))

        with col2:
            st.write('###### Initial Ion Concentrations')
            st.dataframe(ion_vals_init.style.format("{:.1f}"))

            st.write('###### Steady-State Ion Concentrations')
            st.dataframe(ion_vals_ss.style.format("{:.1f}"))

        if itersol_checkbox:

            time, vm_time, chem_time, isim = calculate_iter_results(bes,
                                                                    sim_p,
                                                                    chan_list)

            pmem_dataframe = pd.DataFrame(np.column_stack((time[sim_p.starttime_plot_ind:],
                                                           (1/p.base_pmem)*isim.pmem_time[sim_p.starttime_plot_ind:])),
                                          columns=['Time', 'P_Na', 'P_K', 'P_Cl'])

            st.line_chart(pmem_dataframe, x='Time', y=['P_Na', 'P_K', 'P_Cl'])

            vmem_dataframe = pd.DataFrame(np.column_stack((time[sim_p.starttime_plot_ind:],
                                                           1e3*vm_time[sim_p.starttime_plot_ind:])),
                                          columns=['Time', l.Vmem_o])

            st.line_chart(vmem_dataframe, x='Time', y=l.Vmem_o)
            # st.dataframe(vmem_dataframe.style.format("{:.1f}"))

        # # Iterative solver results:
        # if itersol_checkbox:
        #
        #     time_props, volt_timedat, chem_timedat = calculate_iter_results(p)
        #
        #     st.write("#### Iterative Simulation Results")
        #
        #     # time = time_props['time_vect']
        #
        #     st.write('###### Bioelectrical Potentials')
        #
        #     all_V_series = [l.Vmem_o]
        #     all_chem_series = [l.Na_in_o, l.K_in_o, l.Cl_in_o]
        #
        #     shown_V_series = []
        #     shown_chem_series = []
        #
        #     show_Vmem = st.checkbox(l.Vmem, value=True,
        #                                 help=f'Show membrane potential, {l.Vmem}, on the graph?')
        #
        #     # Columns to arrange voltage series checkboxes horizontally:
        #     vc0, vc1, vc2, vc3, vc4, vc5 = st.columns(6)
        #
        #     with vc0:
        #         show_Ved_Na = st.checkbox(l.Ved_Na,
        #                                   value=False,
        #                                   help=f'Show Na+ electrochemical driving potential,'
        #                                        f' {l.Ved_Na}, on the graph?')
        #     with vc1:
        #         show_Ved_K = st.checkbox(l.Ved_K, value=False,
        #                                  help=f'Show K+ electrochemical driving potential, '
        #                                       f'{l.Ved_K}, on the graph?')
        #
        #     with vc2:
        #         show_Ved_Cl = st.checkbox(l.Ved_Cl, value=False,
        #                                   help=f'Show Cl- electrochemical driving potential, '
        #                                        f'{l.Ved_Cl}, on the graph?')
        #
        #     with vc3:
        #         show_Vrev_Na = st.checkbox(l.Vrev_Na, value=False,
        #                                    help=f'Show Na+ reversal potential, '
        #                                         f'{l.Vrev_Na}, on the graph?')
        #     with vc4:
        #         show_Vrev_K = st.checkbox(l.Vrev_K, value=False,
        #                                   help=f'Show K+ reversal potential, {l.Vrev_K}, on the graph?')
        #
        #     with vc5:
        #         show_Vrev_Cl = st.checkbox(l.Vrev_Cl, value=False,
        #                                    help=f'Show Cl- reversal potential, {l.Vrev_Cl}, on the graph?')
        #
        #     V_series_bools = [show_Vmem, show_Ved_Na, show_Ved_K, show_Ved_Cl, show_Vrev_Na, show_Vrev_K, show_Vrev_Cl]
        #
        #     for vbool, vname in zip(V_series_bools, all_V_series):
        #         if vbool:
        #             shown_V_series.append(vname)
        #
        #     # st.write(f'l.time: {l.time}')
        #     # st.write(f'volt_timedat: {volt_timedat}')
        #     st.line_chart(volt_timedat,
        #                   x=l.time,
        #                   y=shown_V_series,
        #                   use_container_width=True)
        #
        #     st.write('###### Intracellular Ion Concentrations')
        #
        #     # Columns to arrange concentration series checkboxes horizontally:
        #     cc0, cc1, cc2 = st.columns(3)
        #
        #     with cc0:
        #         show_Na_in = st.checkbox(l.Na_in, value=True,
        #                                  help=f'Show Na+ concentration in the cytoplasm, {l.Na_in}, on the graph?')
        #
        #     with cc1:
        #         show_K_in = st.checkbox(l.K_in, value=True, help=f'Show K+ concentration in the cytoplasm, '
        #                                                          f'{l.K_in}, on the graph?')
        #
        #     with cc2:
        #         show_Cl_in = st.checkbox(l.Cl_in, value=True, help=f'Show Cl- concentration in the cytoplasm,'
        #                                                            f' {l.Cl_in}, on the graph?')
        #
        #     chem_series_bools = [show_Na_in, show_K_in, show_Cl_in]
        #
        #     for cbool, cname in zip(chem_series_bools, all_chem_series):
        #         if cbool:
        #             shown_chem_series.append(cname)
        #
        #     st.line_chart(chem_timedat,
        #                   x=l.time,
        #                   y=shown_chem_series,
        #                   use_container_width=True)

    with tab3:
        col1, col2 = st.columns(2)

        with col1:

            if not NaKATP_on and not NaKCl_on and not KCl_on:
                st.write(f'Without active transport, all ions move by electrodiffusion '
                         f'to neutralize electrochemical gradients. Therefore, at steady state all '
                         f'electrochemical driving potentials are zero for all ions and {l.Vmem} '
                         f'depolarizes to near zero. '
                         f'At steady-state, changing membrane permeability to ions has no effect on {l.Vmem}.'
                         )

            elif not NaKATP_on and NaKCl_on and not KCl_on:
                st.write(f'Without primary active transport by the {l.NaK_pump}, '
                         f'there is no electrochemical driving potential on '
                         f'{l.Na}, meaning no energy to operate the {l.NaKCl_cotrans}. '
                         f'Therefore, at steady state, the '
                         f'electrochemical driving potentials for all ions remains zero, '
                         f'and {l.Vmem} remains depolarized near zero. '
                         f'At steady-state, changing membrane permeability to '
                         f'ions has no effect on {l.Vmem}.'
                         )

            elif not NaKATP_on and not NaKCl_on and KCl_on:
                st.write(f'Without primary active transport by the {l.NaK_pump}, '
                         f'there is no electrochemical driving potential on '
                         f'{l.K}, meaning no energy to operate the {l.KCl_symp}. '
                         f'Therefore at steady state, the '
                         f'electrochemical driving potentials for all ions remains zero, '
                         f'and {l.Vmem} remains depolarized near zero. '
                         f'At steady-state, changing membrane permeability to '
                         f'ions has no effect on {l.Vmem}.'
                         )

            elif not NaKATP_on and NaKCl_on and KCl_on:
                st.write(f'Without primary active transport by the {l.NaK_pump}, '
                         f'there is no electrochemical driving potential on '
                         f'{l.Na} or {l.K}, meaning no energy to operate '
                         f'the {l.NaKCl_cotrans} or the {l.KCl_symp}. '
                         f'Therefore at steady state, the '
                         f'electrochemical driving potentials for all ions remains zero, '
                         f'and {l.Vmem} remains depolarized near zero. '
                         f'At steady-state, changing membrane permeability to '
                         f'ions has no effect on {l.Vmem}.'
                         )

            elif NaKATP_on and not NaKCl_on and not KCl_on:
                st.write(f'The ubiquitous {l.NaK_pump} generates {l.Vmem} '
                         f'hyperpolarization by moving 3 {l.Na} out of the cell '
                         f'while bringing only 2 {l.K} in. The resulting electrochemical gradient '
                         f'for {l.Na} strongly favours movement of {l.Na} into the cell to neutralize both '
                         f'a transmembrane concentration and charge gradient. Therefore, {l.Na} has a '
                         f'strongly negative electrochemical driving force ({l.Ved_Na} = {elec_vals_ss.iloc[4, 0]} mV). '
                         f'When {l.Na} enters cells through open {l.Na} channels, {l.Vmem} depolarizes due '
                         f'to the entry of positive charge to the cell.'
                         )

                st.write(f'The electrochemical gradient '
                         f'for {l.K} favours movement of {l.K} out of the cell, against the electrical gradient, '
                         f'to neutralize the strong transmembrane {l.K} concentration gradient created by the {l.NaK_pump}. '
                         f'Therefore, {l.K} has a moderately '
                         f'positive electrochemical driving force ({l.Ved_K} = {elec_vals_ss.iloc[5, 0]} mV). '
                         f'When {l.K} enters cells through open {l.K} channels, {l.Vmem} hyperpolarizes '
                         f'due to the exit of positive charge.')

                st.write(f'When {l.Cl} is not subject to active transport, then when in steady-state, the '
                         f'electrochemical driving force on {l.Cl} will be zero '
                         f'({l.Ved_Cl} = {elec_vals_ss.iloc[6, 0]} mV) with {l.Cl} being in equilibrium with '
                         f'the transmembrane electrical and concentration gradient. When {l.Cl} channels are opened, '
                         f'there is no flux of {l.Cl} and no effect on cell {l.Vmem}.')

                st.write(f'Membrane ion permeabilities are taken to relate to leak channel '
                         f'expression levels. Therefore, increasing membrane {l.Na} and '
                         f'{l.K} permeability in the sidebar Model Variables should show '
                         f'{l.Vmem} depolarization and hyperpolarization, respectively.')

            elif NaKATP_on and NaKCl_on and not KCl_on:
                st.write(f'The ubiquitous {l.NaK_pump} generates {l.Vmem} '
                         f'hyperpolarization by moving 3 {l.Na} out of the cell '
                         f'while bringing only 2 {l.K} in. The resulting electrochemical gradient '
                         f'for {l.Na} strongly favours movement of {l.Na} into the cell to neutralize both '
                         f'a transmembrane concentration and charge gradient. Therefore, {l.Na} has a '
                         f'strongly negative electrochemical driving force ({l.Ved_Na} = {elec_vals_ss.iloc[4, 0]} mV). '
                         f'When {l.Na} enters cells through open {l.Na} channels, {l.Vmem} depolarizes due '
                         f'to the entry of positive charge to the cell.'
                         )

                st.write(f'The electrochemical gradient '
                         f'for {l.K} favours movement of {l.K} out of the cell, against the electrical gradient, '
                         f'to neutralize the strong transmembrane {l.K} concentration gradient created by the {l.NaK_pump}. '
                         f'Therefore, {l.K} has a moderately '
                         f'positive electrochemical driving force ({l.Ved_K} = {elec_vals_ss.iloc[5, 0]} mV). '
                         f'When {l.K} enters cells through open {l.K} channels, {l.Vmem} hyperpolarizes '
                         f'due to the exit of positive charge.')

                st.write(f'The {l.NaKCl_cotrans} uses the impetus for {l.Na} to enter the cell to bring '
                         f'{l.K} and {l.Cl} with it in secondary active transport. '
                         f'When {l.Cl} is subject to active transport by the {l.NaKCl_cotrans}, '
                         f'then when in steady-state the '
                         f'electrochemical driving force on {l.Cl} will be negative '
                         f'({l.Ved_Cl} = {elec_vals_ss.iloc[6, 0]} mV) with an impetus for {l.Cl} to leave the cell '
                         f'in favour with transmembrane electrical and concentration gradients. '
                         f'When {l.Cl} channels are opened, '
                         f'this outward flux of {l.Cl} leads to {l.Vmem} depolarization.')

                st.write(f'Membrane ion permeabilities are taken to relate to leak channel '
                         f'expression levels. Therefore, increasing membrane {l.Na} and '
                         f'{l.K} permeability in the sidebar Model Variables should show '
                         f'{l.Vmem} depolarization and hyperpolarization, respectively. '
                         f'With the {l.NaKCl_cotrans}, increasing {l.Cl} permeability should show '
                         f'{l.Vmem} depolarization.')

            elif NaKATP_on and not NaKCl_on and KCl_on:
                st.write(f'The ubiquitous {l.NaK_pump} generates {l.Vmem} '
                         f'hyperpolarization by moving 3 {l.Na} out of the cell '
                         f'while bringing only 2 {l.K} in. The resulting electrochemical gradient '
                         f'for {l.Na} strongly favours movement of {l.Na} into the cell to neutralize both '
                         f'a transmembrane concentration and charge gradient. Therefore, {l.Na} has a '
                         f'strongly negative electrochemical driving force ({l.Ved_Na} = {elec_vals_ss.iloc[4, 0]} mV). '
                         f'When {l.Na} enters cells through open {l.Na} channels, {l.Vmem} depolarizes due '
                         f'to the entry of positive charge to the cell.'
                         )

                st.write(f'The electrochemical gradient '
                         f'for {l.K} favours movement of {l.K} out of the cell, against the electrical gradient, '
                         f'to neutralize the strong transmembrane {l.K} concentration gradient created by the {l.NaK_pump}. '
                         f'Therefore, {l.K} has a moderately '
                         f'positive electrochemical driving force ({l.Ved_K} = {elec_vals_ss.iloc[5, 0]} mV). '
                         f'When {l.K} enters cells through open {l.K} channels, {l.Vmem} hyperpolarizes '
                         f'due to the exit of positive charge.')

                st.write(f'The {l.KCl_symp} uses the impetus for {l.K} to exit the cell to bring '
                         f'{l.Cl} with it in secondary active transport. '
                         f'When {l.Cl} is subject to active transport by the {l.KCl_symp}, '
                         f'then when in steady-state the '
                         f'electrochemical driving force on {l.Cl} will be positive '
                         f'({l.Ved_Cl} = {elec_vals_ss.iloc[6, 0]} mV) with an impetus for {l.Cl} to enter the cell '
                         f'in favour with transmembrane concentration gradients. '
                         f'When {l.Cl} channels are opened, '
                         f'this inward flux of {l.Cl} leads to {l.Vmem} hyperpolarization.')

                st.write(f'Membrane ion permeabilities are taken to relate to leak channel '
                         f'expression levels. Therefore, increasing membrane {l.Na} and '
                         f'{l.K} permeability in the sidebar Model Variables should show '
                         f'{l.Vmem} depolarization and hyperpolarization, respectively. '
                         f'With the {l.KCl_symp}, increasing {l.Cl} permeability should show '
                         f'{l.Vmem} hyperpolarization.')

            elif NaKATP_on and NaKCl_on and KCl_on:
                st.write(f'The ubiquitous {l.NaK_pump} generates {l.Vmem} '
                         f'hyperpolarization by moving 3 {l.Na} out of the cell '
                         f'while bringing only 2 {l.K} in. The resulting electrochemical gradient '
                         f'for {l.Na} strongly favours movement of {l.Na} into the cell to neutralize both '
                         f'a transmembrane concentration and charge gradient. Therefore, {l.Na} has a '
                         f'strongly negative electrochemical driving force ({l.Ved_Na} = {elec_vals_ss.iloc[4, 0]} mV). '
                         f'When {l.Na} enters cells through open {l.Na} channels, {l.Vmem} depolarizes due '
                         f'to the entry of positive charge to the cell.'
                         )

                st.write(f'The electrochemical gradient '
                         f'for {l.K} favours movement of {l.K} out of the cell, against the electrical gradient, '
                         f'to neutralize the strong transmembrane {l.K} concentration gradient created by the {l.NaK_pump}. '
                         f'Therefore, {l.K} has a moderately '
                         f'positive electrochemical driving force ({l.Ved_K} = {elec_vals_ss.iloc[5, 0]} mV). '
                         f'When {l.K} enters cells through open {l.K} channels, {l.Vmem} hyperpolarizes '
                         f'due to the exit of positive charge.')

                st.write(f'The {l.NaKCl_cotrans} uses the impetus for {l.Na} to enter the cell to bring '
                         f'{l.K} and {l.Cl} with it in secondary active transport. '
                         f'The {l.KCl_symp} uses the impetus for {l.K} to exit the cell to bring '
                         f'{l.Cl} with it in secondary active transport. '
                         f'When {l.Cl} is subject to active transport by both the '
                         f'{l.NaKCl_cotrans} and {l.KCl_symp}, '
                         f'then the steady-state '
                         f'electrochemical driving force on {l.Cl} ({l.Ved_Cl} = {elec_vals_ss.iloc[6, 0]} mV) '
                         f'will depend on which transporter flux is stronger.'
                         )

                st.write(f'Membrane ion permeabilities are taken to relate to leak channel '
                         f'expression levels. Therefore, increasing membrane {l.Na} and '
                         f'{l.K} permeability in the sidebar Model Variables should show '
                         f'{l.Vmem} depolarization and hyperpolarization, respectively. '
                         f'The effect of increasing {l.Cl} permeability will depend on the '
                         f'relative strength of the {l.NaKCl_cotrans} to the {l.KCl_symp}.')


        with col2:
            if not NaKATP_on and not NaKCl_on and not KCl_on:
                cell_net_image_fn = str(get_data_png_cell_network_schematic_4_file())
                cell_net_image = Image.open(cell_net_image_fn)
                st.image(cell_net_image,
                         caption=f'Cellular Bioelectric Network for this Simulation. '
                                 f'There is no active transport, and no electrochemical '
                                 f'driving forces on any ions. '
                                 f'{l.Vmem} = {elec_vals_ss.iloc[0, 0]} mV, '
                                 f'{l.Ved_Na} = {elec_vals_ss.iloc[4, 0]} mV, '
                                 f'{l.Ved_K} = {elec_vals_ss.iloc[5, 0]} mV, '
                                 f'{l.Ved_Cl} = {elec_vals_ss.iloc[6, 0]} mV.',
                         use_column_width='always',
                         output_format="PNG")

            elif not NaKATP_on and NaKCl_on and not KCl_on:
                cell_net_image_fn = str(get_data_png_cell_network_schematic_5_file())
                cell_net_image = Image.open(cell_net_image_fn)
                st.image(cell_net_image,
                         caption=f'Cellular Bioelectric Network for this Simulation. '
                                 f'There is no primary active transport, and no electrochemical '
                                 f'driving forces on any ions. '
                                 f'{l.Vmem} = {elec_vals_ss.iloc[0, 0]} mV, '
                                 f'{l.Ved_Na} = {elec_vals_ss.iloc[4, 0]} mV, '
                                 f'{l.Ved_K} = {elec_vals_ss.iloc[5, 0]} mV, '
                                 f'{l.Ved_Cl} = {elec_vals_ss.iloc[6, 0]} mV.',
                         use_column_width='always',
                         output_format="PNG")

            elif not NaKATP_on and not NaKCl_on and KCl_on:
                cell_net_image_fn = str(get_data_png_cell_network_schematic_6_file())
                cell_net_image = Image.open(cell_net_image_fn)
                st.image(cell_net_image,
                         caption=f'Cellular Bioelectric Network for this Simulation. '
                                 f'There is no primary active transport, and no electrochemical '
                                 f'driving forces on any ions. '
                                 f'{l.Vmem} = {elec_vals_ss.iloc[0, 0]} mV, '
                                 f'{l.Ved_Na} = {elec_vals_ss.iloc[4, 0]} mV, '
                                 f'{l.Ved_K} = {elec_vals_ss.iloc[5, 0]} mV, '
                                 f'{l.Ved_Cl} = {elec_vals_ss.iloc[6, 0]} mV.',
                         use_column_width='always',
                         output_format="PNG")

            elif not NaKATP_on and NaKCl_on and KCl_on:
                cell_net_image_fn = str(get_data_png_cell_network_schematic_7_file())
                cell_net_image = Image.open(cell_net_image_fn)
                st.image(cell_net_image,
                         caption=f'Cellular Bioelectric Network for this Simulation. '
                                 f'There is no primary active transport, and no electrochemical '
                                 f'driving forces on any ions. '
                                 f'{l.Vmem} = {elec_vals_ss.iloc[0, 0]} mV, '
                                 f'{l.Ved_Na} = {elec_vals_ss.iloc[4, 0]} mV, '
                                 f'{l.Ved_K} = {elec_vals_ss.iloc[5, 0]} mV, '
                                 f'{l.Ved_Cl} = {elec_vals_ss.iloc[6, 0]} mV.',
                         use_column_width='always',
                         output_format="PNG")

            elif NaKATP_on and NaKCl_on and KCl_on:
                cell_net_image_fn = str(get_data_png_cell_network_schematic_3_file())
                cell_net_image = Image.open(cell_net_image_fn)
                st.image(cell_net_image,
                         caption=f'Cellular Bioelectric Network for this Simulation. '
                                 f'Active transport is driven by the {l.NaK_pump}, '
                                 f'{l.NaKCl_cotrans}, and {l.KCl_symp}. '
                                 f'{l.Vmem} = {elec_vals_ss.iloc[0, 0]} mV, '
                                 f'{l.Ved_Na} = {elec_vals_ss.iloc[4, 0]} mV, '
                                 f'{l.Ved_K} = {elec_vals_ss.iloc[5, 0]} mV, '
                                 f'{l.Ved_Cl} = {elec_vals_ss.iloc[6, 0]} mV.',
                         use_column_width='always',
                         output_format="PNG")

            elif NaKATP_on and not NaKCl_on and not KCl_on:

                cell_net_image_fn = str(get_data_png_cell_network_schematic_0_file())
                cell_net_image = Image.open(cell_net_image_fn)
                st.image(cell_net_image,
                         caption=f'Cellular Bioelectric Network for this Simulation. '
                                 f'Active transport is driven by the {l.NaK_pump}. '
                                 f'{l.Vmem} = {elec_vals_ss.iloc[0, 0]} mV, '
                                 f'{l.Ved_Na} = {elec_vals_ss.iloc[4, 0]} mV, '
                                 f'{l.Ved_K} = {elec_vals_ss.iloc[5, 0]} mV, '
                                 f'{l.Ved_Cl} = {elec_vals_ss.iloc[6, 0]} mV.',
                         use_column_width='always',
                         output_format="PNG")

            elif NaKATP_on and NaKCl_on and not KCl_on:
                cell_net_image_fn = str(get_data_png_cell_network_schematic_1_file())
                cell_net_image = Image.open(cell_net_image_fn)
                st.image(cell_net_image,
                         caption=f'Cellular Bioelectric Network for this Simulation. '
                                 f'Active transport is driven by the {l.NaK_pump} and '
                                 f'the {l.NaKCl_cotrans}. '
                                 f'{l.Vmem} = {elec_vals_ss.iloc[0, 0]} mV, '
                                 f'{l.Ved_Na} = {elec_vals_ss.iloc[4, 0]} mV, '
                                 f'{l.Ved_K} = {elec_vals_ss.iloc[5, 0]} mV, '
                                 f'{l.Ved_Cl} = {elec_vals_ss.iloc[6, 0]} mV.',
                         use_column_width='always',
                         output_format="PNG")

            elif NaKATP_on and not NaKCl_on and KCl_on:
                cell_net_image_fn = str(get_data_png_cell_network_schematic_2_file())
                cell_net_image = Image.open(cell_net_image_fn)
                st.image(cell_net_image,
                         caption=f'Cellular Bioelectric Network for this Simulation. '
                                 f'Active transport is driven by the {l.NaK_pump} and '
                                 f'the {l.KCl_symp}. '
                                 f'{l.Vmem} = {elec_vals_ss.iloc[0, 0]} mV, '
                                 f'{l.Ved_Na} = {elec_vals_ss.iloc[4, 0]} mV, '
                                 f'{l.Ved_K} = {elec_vals_ss.iloc[5, 0]} mV, '
                                 f'{l.Ved_Cl} = {elec_vals_ss.iloc[6, 0]} mV.',
                         use_column_width='always',
                         output_format="PNG")

            else:
                raise Exception("Did not account for this scenario!")

# ....................{ MAIN ~ run                         }....................
# Run our Streamlit-based web app.
main()
