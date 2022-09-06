#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright (c) 2022-2023 Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.


'''
:mod:`pytest` **global test configuration** (i.e., early-time configuration
guaranteed to be run by :mod:`pytest` *after* passed command-line arguments are
parsed).

:mod:`pytest` implicitly imports *all* functionality defined by this module
into *all* submodules of this subpackage.
'''

# ....................{ FIXTURES                           }....................
# Provide app-specific pytest fixtures required by lower-level tests.
