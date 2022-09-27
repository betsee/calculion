#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright (c) 2022-2023 Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Enumerated types for the Calculion Science Core.
'''

from enum import Enum

class TimeSolverType(Enum):
    '''
    The type of temporal solver to use.
    '''
    explicit = 'explicit'
    implicit = 'implicit'
