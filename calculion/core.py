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

Caveats
----------
**This submodule is intentionally not named** ``calculion.main``. Why? Because
doing so would invite ambiguity with the existing :mod:`calculion.__main__`
dunder submodule and top-level ``main`` script residing in the root of the
repository declaring this package.
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
    from streamlit import (
        title,
    )

    # ..................{ HEADERS                            }..................
    # Human-readable title of this web app.
    title('Calculion')

    # ..................{ LOCALS                             }..................
    #FIXME: Enable after we fix CalculionSim(), please.
    #sim = CalculionSim()

    # ..................{ PLACEHOLDER                        }..................
    #FIXME: Remove this facsimile content after actually implementing something.
    import streamlit as st
    def my_widget(key):
        st.subheader('Hello there!')
        return st.button("Click me " + key)

    # This works in the main area
    clicked = my_widget("first")

    # And within an expander
    my_expander = st.expander("Expand", expanded=True)
    with my_expander:
        clicked = my_widget("second")

    # AND in st.sidebar!
    with st.sidebar:
        clicked = my_widget("third")


# Run our Streamlit-based web app.
main()
