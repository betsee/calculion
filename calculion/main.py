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

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# CAUTION: Avoid importing anything at module scope *EXCEPT* from official
# Python modules in the standard library guaranteed to exist. Subsequent logic
# in the _main() function called below validates third-party runtime
# dependencies of this package to be safely importable, Before performing that
# validation, *NO* other modules are safely importable from.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# ....................{ MAIN                               }....................
#FIXME: Unit test us up, please.
def main() -> None:
    '''
    Core function running this Streamlit-based web app: **Calculion.**
    '''

    # ..................{ IMPORTS                            }..................
    import streamlit as st
    from PIL import Image
    from calculion.science.comp_sys import CompSys
    from calculion.science.params import CalculionParams
    from calculion.science.string_names import StringNames
    from calculion._util.path.utilpathself import (
        get_data_png_cell_network_schematic_0_file,
        get_data_png_cell_network_schematic_1_file,
        get_data_png_cell_network_schematic_2_file,
        get_data_png_cell_network_schematic_3_file,
        get_data_png_membrane_schematic_file,
    )
    from numpy import exp  #, column_stack
    # from pandas import DataFrame
    # from calculion.science.compute import get_steady_state
    from streamlit import (
        title,
        # set_page_config
    )

    # ..................{ HEADERS                            }..................
    # set_page_config(layout="wide") # set a wide page configuration?

    # Human-readable title of this web app.
    title('Calculion')

    # App subtitle, if we want it:
    # st.write('Calculating the *slow* changes of bioelectricity')

    # ..................{ LOCALS                             }..................
    p = CalculionParams()  # Create a default parameters instance
    l = StringNames() # string labels

    # ..................{ SIDEBAR                            }..................
    # The sidebar will contain all widgets to collect user-data for the
    # simulation. Create and name the sidebar.
    st.sidebar.header('Simulation Variables')

    # st.sidebar.write('#### Set simulation variables')
    with st.sidebar:

        st.write('**Initial Conditions**')

        # slider general settings:
        ion_min_val = 0.1  # min concentration of ions in mM
        ion_max_val = 250.0  # max concentration of ions in mM
        ion_slider_step = 0.1  # step-size for ion slider
        memp_min_val = 0.01 # min value for membrane permeability in nm/s
        memp_max_val = 5.0 # max value for membrane permeability in nm/s
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

            p.Na_i = param_widget('Na+ in [mM]',
                                  min_value=ion_min_val,
                                  max_value=ion_max_val,
                                  value=p.Na_i,
                                  step=ion_slider_step,
                                  format='%f',
                                  key='slider_Na_i',
                                  label_visibility='visible',
                                  help='Set the starting value of Na⁺ in the cell.')

            p.K_i = param_widget('K+ in [mM]',
                                 min_value=ion_min_val,
                                 max_value=ion_max_val,
                                 value=p.K_i,
                                 step=ion_slider_step,
                                 format='%f',
                                 key='slider_K_i',
                                 label_visibility='visible',
                                 help='Set the starting value of K⁺ in the cell.')

            p.Cl_i = param_widget('Cl- in [mM]',
                                  min_value=ion_min_val,
                                  max_value=ion_max_val,
                                  value=p.Cl_i,
                                  step=ion_slider_step,
                                  format='%f',
                                  key='slider_Cl_i',
                                  label_visibility='visible',
                                  help='Set the starting value of Cl⁻ in the cell.')

        st.write('**Simulation Parameters**')

        # Define another expander block for specifying extracellular ion concentrations:
        ion_out_block = st.expander("Ion Concentrations Outside Cell", expanded=default_expanded_state)

        with ion_out_block:
            # Automatically reset parameter values in the parameters object p depending on user-selection:
            p.Na_o = param_widget('Na+ out [mM]',
                             min_value=ion_min_val,
                             max_value=ion_max_val,
                             value=p.Na_o,
                             step=ion_slider_step,
                             format='%f',
                             key='slider_Na_o',
                             label_visibility='visible',
                             help='Set the Na⁺ concentration outside the cell.')

            p.K_o = param_widget('K+ out [mM]',
                             min_value=ion_min_val,
                             max_value=ion_max_val,
                             value=p.K_o,
                             step=ion_slider_step,
                             format='%f',
                             key='slider_K_o',
                             label_visibility='visible',
                             help='Set the K⁺ concentration outside the cell.')

            p.Cl_o = param_widget('Cl- out [mM]',
                             min_value=ion_min_val,
                             max_value=ion_max_val,
                             value=p.Cl_o,
                             step=ion_slider_step,
                             format='%f',
                             key='slider_Cl_o',
                             label_visibility='visible',
                             help = 'Set the Cl⁻ concentration outside the cell.')

        memperm_block = st.expander("Cell Membrane Ion Permeabilities", expanded=default_expanded_state)

        with memperm_block:

            p.P_Na_nm = param_widget('Na+ Permeability [nm/s]:',
                               min_value=memp_min_val,
                               max_value=memp_max_val,
                               value=p.P_Na_nm,
                               step=memp_slider_step,
                               format='%f',
                               key='slider_P_Na',
                               label_visibility='visible',
                               help='Set the base membrane permeability to Na⁺.\n'
                                    '\nThis simulates cellular expression of Na⁺ leak channels.'
                               )

            p.P_K_nm = param_widget('K+ Permeability [nm/s]:',
                               min_value=memp_min_val,
                               max_value=memp_max_val,
                               value=p.P_K_nm,
                               step=memp_slider_step,
                               format='%f',
                               key='slider_P_K',
                               label_visibility='visible',
                               help='Set the base membrane permeability to K⁺.\n'
                                    '\nThis simulates cellular expression of K⁺ leak channels.'
                               )

            p.P_Cl_nm = param_widget('Cl- Permeability [nm/s]:',
                               min_value=memp_min_val,
                               max_value=memp_max_val,
                               value=p.P_Cl_nm,
                               step=memp_slider_step,
                               format='%f',
                               key='slider_P_Cl',
                               label_visibility='visible',
                               help='Set the base membrane permeability to Cl⁻.\n'
                                    '\nThis simulates cellular expression of Cl⁻ leak channels.'
                               )

            # Update the main membrane permeabilities used in the simulation to m/s units:
            p.P_Na = 1.0e-9*p.P_Na_nm
            p.P_K = 1.0e-9*p.P_K_nm
            p.P_Cl = 1.0e-9*p.P_Cl_nm

        # Define another expander block for pump and transporter settings:
        pumps_block = st.expander("Ion Pump and Transporter Settings", expanded=default_expanded_state)

        with pumps_block:
            omega_NaK_o = param_widget('Na-K-ATPase Pump Rate [units]',
                             min_value=0.0,
                             max_value=1.0,
                             value=p.omega_NaK*1e12,
                             step=0.01,
                             format='%f',
                             key='slider_omega_NaK',
                             label_visibility='visible',
                             help='Set the maximum rate of the Na-K-ATPase ion pump.')

            # update the na-k-atpase pump rate in units used in the simulation:
            p.omega_NaK = omega_NaK_o*1e-12

            # # Let the user specify whether they want the additional NaKCl cotransporter for secondary transport:
            NaKCl_on = st.checkbox(l.NaKCl_cotrans_o, vale=False, help=f'Cell expresses {l.NaKCl_cotrans} ?')

            if NaKCl_on:
                # Na-K-2Cl cotransporter properties:
                omega_NaKCl_o = param_widget('Na-K-2Cl Cotransporter Rate [units]',
                                 min_value=0.0,
                                 max_value=1.0,
                                 value=p.omega_NaKCl*1e14,
                                 step=0.01,
                                 format='%f',
                                 key='slider_omega_NaKCl',
                                 label_visibility='visible',
                                 help='Set the maximum rate of the Na-K-2Cl cotransporter.')

                # update the na-k-2Cl cotransporter rate in units used in the simulation:
                p.omega_NaKCl = omega_NaKCl_o*1e-14

            # # Let the user specify whether they want the additional KCl symporter for secondary transport:
            KCl_on = st.checkbox(l.KCl_symp_o, vale=False, help=f'Cell expresses {l.KCl_symp} ?')

            if KCl_on:
                # K-Cl symporter properties:
                omega_KCl_o = param_widget('K-Cl Symporter Rate [units]',
                                 min_value=0.0,
                                 max_value=50.0,
                                 value=p.omega_KCl*1e12,
                                 step=0.1,
                                 format='%f',
                                 key='slider_omega_KCl',
                                 label_visibility='visible',
                                 help='Set the maximum rate of the K-Cl symporter.'
                                           )

                # update the K-Cl symporter rate in units used in the simulation:
                p.omega_KCl = omega_KCl_o*1e-12

        # Define another expander block for pump and transporter settings:
        metabolic_block = st.expander("Metabolic Settings", expanded=default_expanded_state)

        with metabolic_block:

            delGATP = param_widget('Gibbs Free Energy ATP Hydrolysis [kJ/mol]',
                             min_value=-34.0,
                             max_value=-28.0,
                             value=p.delGo_ATP*1e-3,
                             step=0.1,
                             format='%f',
                             key='slider_delGATP',
                             label_visibility='visible',
                             help='Set the Gibbs standard free energy for ATP hydrolysis.')

            p.delGo_ATP = delGATP*1e3 # convert to units J/mol from kJ/mol

            # Recalculate equilibrium constant for ATP hydrolysis reaction
            p.Keqm_NaK = exp(p.delGo_ATP / (p.R * p.T))

            p.ATP = param_widget('ATP Concentration [mM]',
                             min_value=0.1,
                             max_value=5.0,
                             value=p.ATP,
                             step=0.01,
                             format='%f',
                             key='slider_ATP',
                             label_visibility='visible',
                             help='Set the concentration of ATP in the cell.')

            p.ADP = param_widget('ADP Concentration [mM]',
                             min_value=0.001,
                             max_value=5.0,
                             value=p.ADP,
                             step=0.001,
                             format='%f',
                             key='slider_ADP',
                             label_visibility='visible',
                             help='Set the concentration of ADP in the cell.')

            p.P = param_widget('P Concentration [mM]',
                             min_value=0.001,
                             max_value=5.0,
                             value=p.P,
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

        # Define a final expander block for simulator settings:
        sim_settings_block = st.expander("Simulation Settings", expanded=default_expanded_state)

        with sim_settings_block:

            # Iterative solver will not be used by default:
            itersol_checkbox = st.checkbox("Use iterative solver",
                                           value=False,
                                           key='checkbox_itersol',
                                           help='Use the iterative solver that integrates the '
                                                'system step-by-step in time?')

            if itersol_checkbox:
                p.iterative_solver = True # Set the iterative solver parameter to True

                # Iterative solver time step:
                p.delta_t = param_widget('Simulation time-step [s]',
                                          min_value=1.0e-3,
                                          max_value=100.0,
                                          value=p.delta_t,
                                          step=0.01,
                                          format='%f',
                                          key='slider_delta_t',
                                          label_visibility='visible',
                                          help='Set the time step for the iterative solver.')

                # Iterative solver max iterations:
                p.N_iter = param_widget('Simulation max iterations',
                                          min_value=10,
                                          max_value=100000,
                                          value=p.N_iter,
                                          step=1,
                                          format='%d',
                                          key='slider_Niter',
                                          label_visibility='visible',
                                          help='Set the maximum number of timesteps that can be run.')

                # Iterative solver convergence tolerance:
                p.steady_state_tol = param_widget('Convergence tolerance',
                                          min_value=1e-20,
                                          max_value=1e-6,
                                          value=p.steady_state_tol,
                                          step=1e-15,
                                          format='%e',
                                          key='slider_tol',
                                          label_visibility='visible',
                                          help='Set the tolerance, below which the simulation will be'
                                               'assumed to be at steady-state.')



            else:
                p.iterative_solver = False # reset the iterative solver parameter to False

    # ..................{ RESULTS                            }..................
    # After collecting parameter values from the user, compute the steady-state
    # values for the bioelectrical system.

    @st.cache # Cache the results of this slower function
    def calculate_ss_results(p):

        sim = CompSys()  # Create an instance of the main computational simulator

        params_vect_o, consts_vect = sim.collect_params(p) # Get initial values of the computational simulator

        time_properties_dict = {}
        volt_timedat = {}
        chem_timedat = {}

        if p.iterative_solver is False:
            params_vect = sim.calc_steady_state(p)

        else:
            (params_vect_time,
             opti_funk_time,
             time_vect,
             print_message) = sim.calc_timestepped(p, N_iter=p.N_iter,
                                                    del_t=p.delta_t,
                                                    ti = 0.0,
                                                    tol=p.steady_state_tol)
            params_vect = params_vect_time[-1] # Get the last time-frame as the final parameters

            # Save time-dependent properties to the time_properties_dict:
            time_properties_dict['params_vect_time'] = params_vect_time
            time_properties_dict['opti_funk_time'] = opti_funk_time
            time_properties_dict['time_vect'] = time_vect

            # Use the time-dependent parameters to compute all electrical and chem properties as a function of time:
            volt_timedat = sim.calc_elec_param_set(params_vect_time, consts_vect, time_vect)
            chem_timedat = sim.calc_chem_param_set(params_vect_time, consts_vect, time_vect)

        # Get the concentration ss dataframe:
        ion_ss = sim.return_chem_props_dict(params_vect, consts_vect)

        # Get the electrical properties dataframe:
        elec_ss = sim.return_elec_props_dict(params_vect, consts_vect)


        return ion_ss, elec_ss, time_properties_dict, volt_timedat, chem_timedat

    # In the main area present the results of the simulation:
    ion_vals_ss, elec_vals_ss, time_props, volt_timedat, chem_timedat = calculate_ss_results(p)

    # ..................{ TABS                               }..................
    # Split the main area into these three tabs:
    # * The "Introduction" tab (tab1) will display a write-up of the theory
    #   behind Calculion.
    # * The "Simulation Results" tab will display the simulation results as
    #   tables and charts.
    # * The "Bioelectrical Network" tab will show a graphical depiction of the
    #   bioelectrical network.
    tab1, tab2, tab3 = st.tabs([
        'Introduction', 'Simulation Results', 'Bioelectrical Network'])

    with tab1:
        st.write('### Why Calculion?')
        # App subtitle, if we want it:
        st.write('#### Calculating the *slow* changes of bioelectricity')
        st.write('Here we will have a preamble describing the motivation and theory behind Calculion.')

        mem_image_fn = str(get_data_png_membrane_schematic_file())
        mem_image = Image.open(mem_image_fn)
        st.image(mem_image,
                 caption='Behold the Cellular Bioelectric Network!',
                 use_column_width='always',
                 output_format="PNG")

    with tab2:
        st.write("### Simulation Results")
        # st.write("*(Alter sidebar Simulation Variables to explore the possibilities...)*")
        col1, col2 = st.columns(2)

        with col1:
            st.write('###### Steady-State Bioelectrical Potentials')
            st.dataframe(elec_vals_ss.style.format("{:.1f}"))

        with col2:
            st.write('###### Steady-State Ion Concentrations')
            st.dataframe(ion_vals_ss.style.format("{:.1f}"))

        # Iterative solver results:
        if itersol_checkbox:
            st.write("#### Iterative Simulation Results")

            # time = time_props['time_vect']

            st.write('###### Bioelectrical Potentials')

            all_V_series = [l.Vmem_o, l.Ved_Na_o, l.Ved_K_o, l.Ved_Cl_o, l.Vrev_Na_o, l.Vrev_K_o, l.Vrev_Cl_o]
            all_chem_series = [l.Na_in_o, l.K_in_o, l.Cl_in_o]

            shown_V_series = []
            shown_chem_series = []

            show_Vmem = st.checkbox(l.Vmem, value=True, help=f'Show membrane potential, {l.Vmem}, on the graph?')
            show_Ved_Na = st.checkbox(l.Ved_Na, value=False, help=f'Show Na+ electrochemical driving potential, {l.Ved_Na}, on the graph?')
            show_Ved_K = st.checkbox(l.Ved_K, value=False, help=f'Show K+ electrochemical driving potential, {l.Ved_K}, on the graph?')
            show_Ved_Cl = st.checkbox(l.Ved_Cl, value=False, help=f'Show Cl- electrochemical driving potential, {l.Ved_Cl}, on the graph?')
            show_Vrev_Na = st.checkbox(l.Vrev_Na, value=False, help=f'Show Na+ reversal potential, {l.Vrev_Na}, on the graph?')
            show_Vrev_K = st.checkbox(l.Vrev_K, value=False, help=f'Show K+ reversal potential, {l.Vrev_K}, on the graph?')
            show_Vrev_Cl = st.checkbox(l.Vrev_Cl, value=False, help=f'Show Cl- reversal potential, {l.Vrev_Cl}, on the graph?')

            V_series_bools = [show_Vmem, show_Ved_Na, show_Ved_K, show_Ved_Cl, show_Vrev_Na, show_Vrev_K, show_Vrev_Cl]

            for vbool, vname in zip(V_series_bools, all_V_series):
                if vbool:
                    shown_V_series.append(vname)


            st.line_chart(volt_timedat,
                          x=l.time,
                          y=shown_V_series,
                          use_container_width=True)


            st.write('###### Intracellular Ion Concentrations')

            show_Na_in = st.checkbox(l.Na_in, value=True, help=f'Show Na+ concentration in the cytoplasm, {l.Na_in}, on the graph?')
            show_K_in = st.checkbox(l.K_in, value=True, help=f'Show K+ concentration in the cytoplasm, {l.K_in}, on the graph?')
            show_Cl_in = st.checkbox(l.Cl_in, value=True, help=f'Show Cl- concentration in the cytoplasm, {l.Cl_in}, on the graph?')

            chem_series_bools = [show_Na_in, show_K_in, show_Cl_in]

            for cbool, cname in zip(chem_series_bools, all_chem_series):
                if cbool:
                    shown_chem_series.append(cname)

            st.line_chart(chem_timedat,
                          x=l.time,
                          y=shown_chem_series,
                          use_container_width=True)

    with tab3:

        col1, col2 = st.columns(2)

        with col2:

            if NaKCl_on and KCl_on:
                cell_net_image_fn = str(get_data_png_cell_network_schematic_3_file())
                cell_net_image = Image.open(cell_net_image_fn)
                st.image(cell_net_image,
                         # caption='Behold the Cellular Bioelectric Network!',
                         use_column_width='always',
                         output_format="PNG")

            elif not NaKCl_on and not KCl_on:

                cell_net_image_fn = str(get_data_png_cell_network_schematic_0_file())
                cell_net_image = Image.open(cell_net_image_fn)
                st.image(cell_net_image,
                         caption='Behold the Cellular Bioelectric Network!',
                         use_column_width='always',
                         output_format="PNG")

            elif NaKCl_on and not KCl_on:
                cell_net_image_fn = str(get_data_png_cell_network_schematic_1_file())
                cell_net_image = Image.open(cell_net_image_fn)
                st.image(cell_net_image,
                         # caption='Behold the Cellular Bioelectric Network!',
                         use_column_width='always',
                         output_format="PNG")

            elif not NaKCl_on and KCl_on:
                cell_net_image_fn = str(get_data_png_cell_network_schematic_2_file())
                cell_net_image = Image.open(cell_net_image_fn)
                st.image(cell_net_image,
                         # caption='Behold the Cellular Bioelectric Network!',
                         use_column_width='always',
                         output_format="PNG")

            else:
                raise Exception("Did not account for this scenario!")
# ....................{ MAIN ~ run                         }....................
# Run our Streamlit-based web app.
main()
