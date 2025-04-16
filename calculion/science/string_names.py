#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright (c) 2022-2025 Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Creates a set of commonly-used string names for parameters and other values.
'''

from dataclasses import dataclass

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
    Na_in: str = 'Na⁺ in'
    K_in: str = 'K⁺ in'
    Cl_in: str = 'Cl⁻ in'
    Na_out: str = 'Na⁺ out'
    K_out: str = 'K⁺ out'
    Cl_out: str = 'Cl⁻ out'
    NaK_pump: str = 'Na⁺,K⁺ ATPase ion pump'
    NaKCl_cotrans: str = 'Na⁺,K⁺, 2Cl⁻ co-transporter'
    KCl_symp: str = 'K⁺, Cl⁻ symporter'
    Na: str = 'Na⁺'
    K: str = 'K⁺'
    Cl: str = 'Cl⁻'
    PNa: str = 'P_Na'
    PK: str = 'P_K'
    PCl: str = 'P_Cl'

    time: str = 'Time (seconds)'

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
    NaK_pump_o: str = 'Na-K-ATPase ion pump'
    NaKCl_cotrans_o: str = 'Na-K-Cl co-transporter'
    KCl_symp_o: str = 'K-Cl symporter'
