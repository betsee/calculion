#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright (c) 2022-2023 Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Enumerations for the ionotronics module.

'''

from enum import Enum

class IonType(Enum):
    Na = 'Na+'
    K = 'K+'
    Cl = 'Cl-'
    Ca = 'Ca++'
    H = 'H+'
    HCO3 = 'HCO3-'
    Mg = 'Mg++'
    Zn = 'Zn++'
    Fe2 = 'Fe++'
    Fe3 = 'Fe+++'

class ReactionClass(Enum):
    '''
    Enumeration specifying the type of reaction.

    '''
    classic = 'classic'
    transport = 'transport'
    redox = 'redox'

class SystemDimension(Enum):
    '''
    Enumeration specifying the dimension of the system.
    '''
    d2 = '2D'
    d3 = '3D'

class PumpType(Enum):
    '''
    Enumeration for ion pump type.
    '''
    NaKATPase = 'Na,K-ATPase'
    HATPase = 'H-ATPase'
    CaATPase ='Ca-ATPase'
    HKATPase = 'HK-ATPase'

class TransporterType(Enum):
    NaK2Cl = 'NaKCC'
    KCl = 'KCl'

class ChannelType(Enum):
    NaPulse = 'Na Pulse'
    KPulse = 'K Pulse'
    ClPulse = 'Cl Pulse'

class ReactionType(Enum):
    Bicarb_buffer = 'Bicarbonate Buffer'
    ATP_hydrolysis = 'ATP Hydrolysis'
    Metabolism = 'Metabolism'