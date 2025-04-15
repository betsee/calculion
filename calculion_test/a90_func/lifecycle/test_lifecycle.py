#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright (c) 2022-2025 Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
**App lifecycle** (i.e., full process of spinning up, running, and shutting down
this Streamlit-based web app) integration tests.

This submodule integration tests the public API of the
:mod:`calculion` package.
'''

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To raise human-readable test errors, avoid importing from
# package-specific submodules at module scope.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# ....................{ TESTS                              }....................
def test_app_lifecycle(capsys) -> None:
    '''
    Integration test exercising the **maximally trivial app lifecycle** (i.e.,
    workflow spinning up, running, and then immediately shutting down this app).

    Parameters
    ----------
    capsys
        Builtin fixture object permitting pytest-specific capturing (i.e.,
        hiding) of both standard output and error to be temporarily disabled.
    '''

    # Temporarily disable pytest-specific capturing (i.e., hiding) of both
    # standard output and error for the duration of this integration test. Why?
    # Because this test often takes an excessive amount of time and can,
    # moreover, fail to halt in various edge cases outside our control.
    # Capturing output unhelpfully obscures these concerns during local testing.
    with capsys.disabled():
        # Defer test-specific imports.
        #
        # Note that:
        # * Importing "calculion.main" has the beneficial side effect of running
        #   this Streamlit-based web app to completion.
        # * Attempting to run this app via the "streamlit run" subcommand does
        #   *NOT* halt as expected, as that subcommand understandably runs the
        #   passed app indefinitely.
        from calculion import main

        # Immediately return after doing so.
        return
