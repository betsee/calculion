#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright (c) 2022-2023 Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Creates a set of commonly-used string names for parameters and other values.
'''

from beartype import beartype
from dataclasses import dataclass

@beartype
@dataclass(frozen=True)
class StringNames(object):
    '''
    Immutable object defining commonly used string parameter names.
    '''

    # Fancy labels with subscripts and superscripts:
    Vmem: str = 'Vₘ'
    Vrev_Na: str = 'Vᵣ Na⁺'
    Vrev_K: str = 'Vᵣ K⁺'
    Vrev_Cl: str = 'Vᵣ Cl⁻'
    Ved_Na: str = 'Vₑ Na⁺'
    Ved_K: str = 'Vₑ K⁺'
    Ved_Cl: str = 'Vₑ Cl⁻'
    Na_in: str = 'Na⁺ᵢₙ'
    K_in: str = 'K⁺ᵢₙ'
    Cl_in: str = 'Cl⁻ᵢₙ'
    Na_out: str = 'Na⁺ₒᵤₜ'
    K_out: str = 'K⁺ₒᵤₜ'
    Cl_out: str = 'Cl⁻ₒᵤₜ'
    NaKCl_cotrans: str = 'Na⁺,K⁺,2Cl⁻ co-transporter'
    KCl_symp: str = 'K⁺,Cl⁻ symporter'

    time: str = 'Time (hours)'

    # Plain labels:
    Vmem_o: str = 'Vmem'
    Vrev_Na_o: str = 'Vrev Na+'
    Vrev_K_o: str = 'Vrev K+'
    Vrev_Cl_o: str = 'Vrev Cl-'
    Ved_Na_o: str = 'Ved Na+'
    Ved_K_o: str = 'Ved K+'
    Ved_Cl_o: str = 'Ved Cl-'
    Na_in_o: str = 'Na+ in'
    K_in_o: str = 'K+ in'
    Cl_in_o: str = 'Cl- in'
    Na_out_o: str = 'Na+ out'
    K_out_o: str = 'K+ out'
    Cl_out_o: str = 'Cl- out'
    NaKCl_cotrans_o: str = 'Na-K-Cl co-transporter'
    KCl_symp_o: str = 'K-Cl symporter'