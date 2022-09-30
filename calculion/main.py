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
    Script to run the Streamlit-based web app, Calculion.
    '''

    # ..................{ IMPORTS                            }..................
    from calculion.science.params import CalculionParams
    from numpy import exp, column_stack
    from pandas import DataFrame
    from calculion.science.compute import get_steady_state
    from streamlit import (
        title,
        # set_page_config
    )

    # ..................{ HEADERS                            }..................
    # set_page_config(layout="wide") # set a wide page configuration?

    # Human-readable title of this web app.
    title('Calculion')

    # ..................{ LOCALS                             }..................
    p = CalculionParams()  # Create a default parameters instance

    # ..................{ Calculion App             }..................
    import streamlit as st



    # App subtitle, if we want it:
    # st.write('Calculating the *slow* changes of bioelectricity')

    # def my_widget(key):
    #
    #     return st.button("Click me " + key)

    # The sidebar will contain all widgets to collect user-data for the simulation:
    # Create and name the sidebar:
    st.sidebar.header('Simulation Variables')
    # st.sidebar.write('#### Set simulation variables')
    #-----SIDEBAR AREA--------------------------------------------------------------------------------------------------
    with st.sidebar:

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

        memperm_block = st.expander("Cell Membrane Ion Permeabilities", expanded=default_expanded_state)

        with memperm_block:

            p.P_Na_nm = param_widget('Na+ Permeability [nm/s]:',
                               min_value=memp_min_val,
                               max_value=memp_max_val,
                               value=p.P_Na_nm,
                               step=memp_slider_step,
                               format='%f',
                               key='slider_P_Na',
                               label_visibility='visible'
                               )

            p.P_K_nm = param_widget('K+ Permeability [nm/s]:',
                               min_value=memp_min_val,
                               max_value=memp_max_val,
                               value=p.P_K_nm,
                               step=memp_slider_step,
                               format='%f',
                               key='slider_P_K',
                               label_visibility='visible'
                               )

            p.P_Cl_nm = param_widget('Cl- Permeability [nm/s]:',
                               min_value=memp_min_val,
                               max_value=memp_max_val,
                               value=p.P_Cl_nm,
                               step=memp_slider_step,
                               format='%f',
                               key='slider_P_Cl',
                               label_visibility='visible'
                               )

            # Update the main membrane permeabilities used in the simulation to m/s units:
            p.P_Na = 1.0e-9*p.P_Na_nm
            p.P_K = 1.0e-9*p.P_K_nm
            p.P_Cl = 1.0e-9*p.P_Cl_nm

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
                             label_visibility='visible')

            p.K_o = param_widget('K+ out [mM]',
                             min_value=ion_min_val,
                             max_value=ion_max_val,
                             value=p.K_o,
                             step=ion_slider_step,
                             format='%f',
                             key='slider_K_o',
                             label_visibility='visible')

            p.Cl_o = param_widget('Cl- out [mM]',
                             min_value=ion_min_val,
                             max_value=ion_max_val,
                             value=p.Cl_o,
                             step=ion_slider_step,
                             format='%f',
                             key='slider_Cl_o',
                             label_visibility='visible')

        # Define another expander block collecting parameters for intracellular ion concentrations:
        ion_in_block = st.expander("Ion Concentrations Inside Cell", expanded=default_expanded_state)

        with ion_in_block:

            p.Na_i = param_widget('Na+ in [mM]',
                             min_value=ion_min_val,
                             max_value=ion_max_val,
                             value=p.Na_i,
                             step=ion_slider_step,
                             format='%f',
                             key='slider_Na_i',
                             label_visibility='visible')

            p.K_i = param_widget('K+ in [mM]',
                             min_value=ion_min_val,
                             max_value=ion_max_val,
                             value=p.K_i,
                             step=ion_slider_step,
                             format='%f',
                             key='slider_K_i',
                             label_visibility='visible')

            p.Cl_i = param_widget('Cl- in [mM]',
                             min_value=ion_min_val,
                             max_value=ion_max_val,
                             value=p.Cl_i,
                             step=ion_slider_step,
                             format='%f',
                             key='slider_Cl_i',
                             label_visibility='visible')

        # Define another expander block for pump and transporter settings:
        pumps_block = st.expander("Ion Pump and Transporter Settings", expanded=default_expanded_state)

        with pumps_block:
            omega_NaK_o = param_widget('Na,K-ATPase Pump Rate [units]',
                             min_value=0.0,
                             max_value=1.0,
                             value=p.omega_NaK*1e12,
                             step=0.01,
                             format='%f',
                             key='slider_omega_NaK',
                             label_visibility='visible')

            # update the na-k-atpase pump rate in units used in the simulation:
            p.omega_NaK = omega_NaK_o*1e-12

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
                             label_visibility='visible')

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
                             label_visibility='visible')

            p.ADP = param_widget('ADP Concentration [mM]',
                             min_value=0.001,
                             max_value=5.0,
                             value=p.ADP,
                             step=0.001,
                             format='%f',
                             key='slider_ADP',
                             label_visibility='visible')

            p.P = param_widget('P Concentration [mM]',
                             min_value=0.001,
                             max_value=5.0,
                             value=p.P,
                             step=0.001,
                             format='%f',
                             key='slider_P',
                             label_visibility='visible')

        # Define an expander block for cell and env settings:
        cellenv_settings_block = st.expander("Cell and Environmental Settings", expanded=default_expanded_state)
        with cellenv_settings_block:
            p.T_C = param_widget('Temperature [degrees C]',
                             min_value=1.0,
                             max_value=60.0,
                             value=p.T_C,
                             step=1.0,
                             format='%f',
                             key='slider_T',
                             label_visibility='visible')

            # Update temperature in Kelvin used in simulations:
            p.T = p.T_C + 273.15 # temp in Kelvin

            p.r_cell_um = param_widget('Cell radius [um]',
                             min_value=1.0,
                             max_value=25.0,
                             value=p.r_cell_um,
                             step=1.0,
                             format='%f',
                             key='slider_rcell',
                             label_visibility='visible')

            # update r_cell to be in meters for simulations:
            p.r_cell = p.r_cell_um*1e-6

        # Define a final expander block for simulator settings:
        sim_settings_block = st.expander("Simulation Settings", expanded=default_expanded_state)

        with sim_settings_block:

            # Iterative solver will be used by default
            itersol_checkbox = st.checkbox("Use iterative solver", value=True, key='checkbox_itersol')

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
                                          label_visibility='visible')

                # Iterative solver max iterations:
                p.N_iter = param_widget('Simulation max iterations',
                                          min_value=10,
                                          max_value=100000,
                                          value=p.N_iter,
                                          step=1,
                                          format='%d',
                                          key='slider_Niter',
                                          label_visibility='visible')

                if st.checkbox('Update environment concentrations', value=False, key='checkbox_env_con'):
                    p.update_env = True

                    # Extracellular space properties
                    p.d_ecm_um = param_widget('Extracellular space thickness [um]',
                                              min_value=0.1,
                                              max_value=1000.0,
                                              value=p.d_ecm_um,
                                              step=0.1,
                                              format='%f',
                                              key='slider_decm',
                                              label_visibility='visible')

                    # Update extracellular space width to meters for simulations:
                    p.d_ecm = p.d_ecm_um * 1e-6

                else:
                    p.update_env = False

                if st.checkbox('Fully dynamic Vmem updates', value=False, key='checkbox_dyn_vmem'):
                    p.quasi_static_vmem = False

                    # Membrane relative electrical permittivity space properties
                    p.e_r = param_widget('Membrane relative permittivity',
                                         min_value=1.0,
                                         max_value=100.0,
                                         value=p.e_r,
                                         step=1.0,
                                         format='%f',
                                         key='slider_mem_er',
                                         label_visibility='visible')

                    # Update membrane thickness:
                    p.d_mem_nm = param_widget('Membrane thickness [nm]',
                                              min_value=1.0,
                                              max_value=10.0,
                                              value=p.d_mem_nm,
                                              step=1.0,
                                              format='%f',
                                              key='slider_dmem',
                                              label_visibility='visible')

                    # Update membrane thickness to meters for simulations:
                    p.d_mem = p.d_mem_nm * 1e-9

                else:
                    p.quasi_static_vmem = True

            else:
                p.iterative_solver = False # reset the iterative solver parameter to False


    #-----MAIN RESULTS AREA---------------------------------------------------------------------------------------------
    # After collecting parameter values from the user, compute the steady-state values for the bioelectrical system:

    @st.cache # Cache the results of this slower function
    def calculate_results(p):
        ion_ss, elec_ss, time_sols = get_steady_state(p,
                                       iterative_sol=p.iterative_solver,
                                       update_env=p.update_env,
                                       quasi_static_vmem=p.quasi_static_vmem)

        return ion_ss, elec_ss, time_sols

    # In the main area present the results of the simulation:
    ion_vals_ss, elec_vals_ss, time_sols = calculate_results(p)

    # Split the main area into three tabs:
    # Introduction tab (tab1) will display a write-up of the theory behind Calculion
    # Simulation Results tab will display the simulation results as tables and charts
    # Bioelectrical Network tab will show a graphical depiction of the bioelectrical network.
    tab1, tab2, tab3 = st.tabs(["Introduction", "Simulation Results", "Bioelectrical Network"])

    with tab1:
        st.write("### Why Calculion?")
        # App subtitle, if we want it:
        st.write('#### Calculating the *slow* changes of bioelectricity')
        st.write('Here we will have a preamble describing the motivation and theory behind Calculion.')

    with tab2:
        st.write("### Simulation Results")
        # st.write("*(Alter sidebar Simulation Variables to explore the possibilities...)*")
        col1, col2 = st.columns(2)

        with col1:
            st.write('###### Steady-State Bioelectrical Potentials')
            st.dataframe(elec_vals_ss)

        with col2:
            st.write('###### Steady-State Ion Concentrations')
            st.dataframe(ion_vals_ss)

        # Iterative solver results:
        if itersol_checkbox:
            st.write("#### Iterative Simulation Results")

            # Stack the results into a data stack and display as line chart:
            volt_dat = column_stack((time_sols.time/3600,
                                     1e3*time_sols.V_mem_time,
                                     1e3*time_sols.V_rev_Na_time,
                                     1e3*time_sols.V_rev_K_time,
                                     1e3*time_sols.V_rev_Cl_time
                                     ))

            volt_datframe =DataFrame(volt_dat, columns=['Time (hours)',
                                                        'V_mem',
                                                        'V_rev Na+',
                                                        'V_rev K+',
                                                        'V_rev Cl-'])

            st.write('###### Bioelectrical Potentials')
            st.line_chart(volt_datframe,
                          x='Time (hours)',
                          y=('V_mem', 'V_rev Na+', 'V_rev K+', 'V_rev Cl-'),
                          use_container_width=True)

            # Stack the results of ion concentration changes as a function of time:
            ion_dat = column_stack((time_sols.time / 3600,
                             time_sols.Na_i_time,
                             time_sols.K_i_time,
                             time_sols.Cl_i_time))

            ion_datframe = DataFrame(ion_dat, columns = ['Time (hours)',
                                                         'Na+ in',
                                                         'K+ in',
                                                         'Cl- in'])

            st.write('###### Intracellular Ion Concentrations')
            st.line_chart(ion_datframe,
                          x='Time (hours)',
                          y=('Na+ in', 'K+ in', 'Cl- in'),
                          use_container_width=True)

            if p.update_env:
                # Create a final line chart of Extracellular ion concentrations
                # Stack the results of ion concentration changes as a function of time:
                ion_dat = column_stack((time_sols.time / 3600,
                                        time_sols.Na_o_time,
                                        time_sols.K_o_time,
                                        time_sols.Cl_o_time))

                ion_datframe = DataFrame(ion_dat, columns=['Time (hours)',
                                                           'Na+ out',
                                                           'K+ out',
                                                           'Cl- out'])

                st.write('###### Extracellular Ion Concentrations')
                st.line_chart(ion_datframe,
                              x='Time (hours)',
                              y=('Na+ out', 'K+ out', 'Cl- out'),
                              use_container_width=True)




# Run our Streamlit-based web app.
main()