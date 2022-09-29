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
    Run our Streamlit-based web app.
    '''

    # ..................{ IMPORTS                            }..................
    from calculion.science.optimize import Optimizer
    from calculion.science.params import CalculionParams
    from numpy import exp
    from calculion.science.compute import get_steady_state
    from streamlit import (
        title,
    )

    # ..................{ HEADERS                            }..................
    # Human-readable title of this web app.
    title('Calculion')

    # ..................{ LOCALS                             }..................
    #FIXME: Enable after we fix CalculionSim(), please.
    #sim = CalculionSim()
    p = CalculionParams()  # Create a parameters instance for the web app

    # Populate Calculion tables and cells with default parameters from p:



    # ..................{ Calculion App             }..................
    #FIXME: Remove this facsimile content after actually implementing something.
    import streamlit as st


    def my_widget(key):
        st.subheader('Hello there!')
        return st.button("Click me " + key)

    # The sidebar will contain all widgets to collect user-data for the simulation:
    # Create and name the sidebar:
    st.sidebar.header('Simulation Variables')
    # st.sidebar.write('#### Set simulation variables')

    with st.sidebar:

        # slider general settings:
        ion_min_val = 0.1  # min concentration of ions in mM
        ion_max_val = 250.0  # max concentration of ions in mM
        ion_slider_step = 0.1  # step-size for ion slider
        memp_min_val = 0.1 # min value for membrane permeability in nm/s
        memp_max_val = 5.0 # max value for membrane permeability in nm/s
        memp_slider_step = 0.01 # step-size for mem perm slider

        # FIXME: in the code below, each named sidebar section could be an 'expander' panel that is
        # initially hidden (displaying only the title, e.g. 'Cell Membrane Ion Permeabilities") and
        # that can be expanded to reveal the sliders with a click.

        st.sidebar.write('#### Cell Membrane Ion Permeabilities:')

        p.P_Na_nm = st.slider('Na+ Permeability [nm/s]:',
                           min_value=memp_min_val,
                           max_value=memp_max_val,
                           value=p.P_Na_nm,
                           step=memp_slider_step,
                           format='%f',
                           key='slider_P_Na',
                           label_visibility='visible'
                           )

        p.P_K_nm = st.slider('K+ Permeability [nm/s]:',
                           min_value=memp_min_val,
                           max_value=memp_max_val,
                           value=p.P_K_nm,
                           step=memp_slider_step,
                           format='%f',
                           key='slider_P_K',
                           label_visibility='visible'
                           )

        p.P_Cl_nm = st.slider('Cl- Permeability [nm/s]:',
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


        st.sidebar.write('#### Ion Concentrations Outside Cell:')
        # Automatically reset parameter values in the parameters object p depending on user-selection:
        p.Na_o = st.slider('Na+ out [mM]',
                         min_value=ion_min_val,
                         max_value=ion_max_val,
                         value=p.Na_o,
                         step=ion_slider_step,
                         format='%f',
                         key='slider_Na_o',
                         label_visibility='visible')

        p.K_o = st.slider('K+ out [mM]',
                         min_value=ion_min_val,
                         max_value=ion_max_val,
                         value=p.K_o,
                         step=ion_slider_step,
                         format='%f',
                         key='slider_K_o',
                         label_visibility='visible')

        p.Cl_o = st.slider('Cl- out [mM]',
                         min_value=ion_min_val,
                         max_value=ion_max_val,
                         value=p.K_o,
                         step=ion_slider_step,
                         format='%f',
                         key='slider_Cl_o',
                         label_visibility='visible')

        st.sidebar.write('#### Ion Concentrations Inside Cell:')

        p.Na_i = st.slider('Na+ in [mM]',
                         min_value=ion_min_val,
                         max_value=ion_max_val,
                         value=p.Na_i,
                         step=ion_slider_step,
                         format='%f',
                         key='slider_Na_i',
                         label_visibility='visible')

        p.K_i = st.slider('K+ in [mM]',
                         min_value=ion_min_val,
                         max_value=ion_max_val,
                         value=p.K_i,
                         step=ion_slider_step,
                         format='%f',
                         key='slider_K_i',
                         label_visibility='visible')

        p.Cl_i = st.slider('Cl- in [mM]',
                         min_value=ion_min_val,
                         max_value=ion_max_val,
                         value=p.K_i,
                         step=ion_slider_step,
                         format='%f',
                         key='slider_Cl_i',
                         label_visibility='visible')

        st.sidebar.write('#### Ion Pump and Transporter Settings')
        omega_NaK_o = st.slider('Na,K-ATPase Pump Rate [units]',
                         min_value=0.0,
                         max_value=1.0,
                         value=p.omega_NaK*1e12,
                         step=0.01,
                         format='%f',
                         key='slider_omega_NaK',
                         label_visibility='visible')

        # update the na-k-atpase pump rate in units used in the simulation:
        p.omega_NaK = omega_NaK_o*1e-12

        st.sidebar.write('#### Metabolic Settings')

        delGATP = st.slider('Gibbs Free Energy ATP Hydrolysis [kJ/mol]',
                         min_value=-34.0,
                         max_value=-28.0,
                         value=p.delGo_ATP*1e-3,
                         step=0.1,
                         format='%f',
                         key='slider_delGATP',
                         label_visibility='visible')

        p.delGo_ATP = delGATP*1e3
        # Recalculate equilibrium constant for ATP hydrolysis reaction
        p.Keqm_NaK = exp(p.delGo_ATP / (p.R * p.T))

        p.ATP = st.slider('ATP Concentration [mM]',
                         min_value=0.1,
                         max_value=5.0,
                         value=p.ATP,
                         step=0.01,
                         format='%f',
                         key='slider_ATP',
                         label_visibility='visible')

        p.ADP = st.slider('ADP Concentration [mM]',
                         min_value=0.001,
                         max_value=5.0,
                         value=p.ADP,
                         step=0.001,
                         format='%f',
                         key='slider_ADP',
                         label_visibility='visible')

        p.P = st.slider('P Concentration [mM]',
                         min_value=0.001,
                         max_value=5.0,
                         value=p.P,
                         step=0.001,
                         format='%f',
                         key='slider_P',
                         label_visibility='visible')

    # In the main area present the results of the simulation:
    st.write('Na out:', p.Na_o, 'Na in:', p.Na_i)
    # This works in the main area
    clicked = my_widget("first")

    # And within an expander
    my_expander = st.expander("Expand", expanded=True)
    with my_expander:
        clicked = my_widget("second")


# Run our Streamlit-based web app.
main()
