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
    Test the public API of the :class:`calculion.science.compute.CalculionSim`
    class.
    '''

    # Defer test-specific imports.
    from calculion.science.compute import CalculionSim

    #FIXME: Actually call methods of this simulator, please. \o/
    # App simulator.
    app_sim = CalculionSim()
