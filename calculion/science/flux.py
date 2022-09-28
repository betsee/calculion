#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright (c) 2022-2023 Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Methods to calculate ion flux.
'''

# ....................{ IMPORTS                            }....................
import numpy as np
from beartype import beartype
from calculion.science.params import CalculionParams
#-------------------------------------------------------------------------------

@beartype
def f_cation_flux(P_cat: float, cat_o: float, cat_i: float, V_mem: float, alpha: float, z_cat: float=1.0):
    '''
    Transmembrane flux for cations (+ ions).
    '''
    flux = P_cat*(cat_o - cat_i*np.exp(z_cat*V_mem*alpha))

    return flux

@beartype
def f_anion_flux(P_ani: float, ani_o: float, ani_i: float, V_mem: float, alpha: float, z_ani: float = -1.0):
    '''
    Transmembrane flux for anions (- ions).
    '''
    flux = -P_ani * (ani_i - ani_o*np.exp(-z_ani*V_mem*alpha))

    return flux

@beartype
def f_Na_flux(P_Na: float, Na_o: float, Na_i: float, V_mem: float, alpha: float):
    '''
    Transmembrane flux for the Na+ ion.
    '''
    flux = P_Na * (Na_o - Na_i * np.exp(V_mem * alpha))
    return flux

@beartype
def f_K_flux(P_K: float, K_o: float, K_i: float, V_mem: float, alpha: float):
    '''
    Transmembrane flux for the K+ ion.
    '''
    flux = P_K * (K_o - K_i * np.exp(V_mem * alpha))
    return flux

@beartype
def f_Cl_flux(P_Cl: float, Cl_o: float, Cl_i: float, V_mem: float, alpha: float):
    '''
    Transmembrane flux for the Cl- ion.
    '''
    flux = -P_Cl * (Cl_i - Cl_o * np.exp(V_mem * alpha))
    return flux

@beartype
def f_NaKpump_flux(ATP: float, ADP: float, P: float,
                   Na_i: float, Na_o: float, K_i: float, K_o: float,
                   V_mem: float, alpha: float,
                   K_eqm_NaK: float, omega_pump: float):
    '''
    Calculate flux for the NaK-ATPase pump.
    '''
    flux = omega_pump * (ADP * P * K_eqm_NaK * (K_i ** 2) * Na_o ** 3 -
                           ATP * (K_o ** 2) * (Na_i ** 3) * np.exp(V_mem * alpha))

    return flux

@beartype
def f_NKCC_flux(Na_i: float, Na_o: float,
                K_i: float, K_o: float,
                Cl_i: float, Cl_o: float,
                omega_nkcc: float):
    '''
    Calculate flux for the Na-K-Cl Cotransporter (NKCC) symporter.
    '''
    flux = omega_nkcc * (K_o * Na_o * Cl_o ** 2 - K_i * Na_i * Cl_i ** 2)

    return flux


def compute_flux_vector(P_Na: float, Na_o: float, Na_i: float,
                        P_K: float, K_o: float, K_i: float,
                        P_Cl: float, Cl_o: float, Cl_i: float,
                        ATP: float, ADP: float, P: float,
                        Keqm_NaK: float, omega_NaK: float,
                        V_mem: float, alpha: float):
    '''
    Given the input variables, computed the flux vector,
    flux_v = [flx_Na, flx_K, flx_Cl, f_X, flx_NaKpump], that governs time-dependent
    changes in the system.
    '''

    flx_Na = f_Na_flux(P_Na, Na_o, Na_i, V_mem, alpha)

    flx_K = f_K_flux(P_K, K_o, K_i, V_mem, alpha)

    flx_Cl = f_Cl_flux(P_Cl, Cl_o, Cl_i, V_mem, alpha)

    flx_NaKpump = f_NaKpump_flux(ATP, ADP, P,
                                 Na_i, Na_o,
                                 K_i, K_o,
                                 V_mem, alpha,
                                 Keqm_NaK, omega_NaK
                                 )

    return np.asarray([flx_Na, flx_K, flx_Cl, flx_NaKpump])

def make_reaction_matrix(p: CalculionParams, quasi_static_vmem: bool=True, update_env: bool=False):
    '''
    Return a matrix for computing the changes to parameters
    for the bioelectric system.

    When this reaction matrix is taken in a dot product with the flux vector,
    flux_v = [f_Na, f_K, f_Cl, f_NaKpump], the time change vector,
    dparams = [dNa_i, dNa_o, dK_i, dK_o, dCl_i, dCl_o] (for include_vmem = False) and
    dparams = [dNa_i, dNa_o, dK_i, dK_o, dCl_i, dCl_o, dVmem] (for include_vmem = True)
    is returned.
    '''

    div_i = 2 / p.r_cell # divergence term cell compartment
    div_o = -(2 * p.r_cell) / (p.r_env**2 - p.r_cell**2) # divergence term env compartment
    cm = p.c_mem # membrane patch capacitance [F/m^2]

    if quasi_static_vmem:
        if update_env:
            Msys = np.asarray([[div_i, 0, 0, 3*div_i],
                               [div_o, 0, 0, 3*div_o],
                               [0, div_i, 0, -2*div_i],
                               [0, div_o, 0, -2*div_o],
                               [0, 0, div_i, 0],
                               [0, 0, div_o, 0],
                               ])

        else:
            Msys = np.asarray([[div_i, 0, 0, 3*div_i],
                               [0, div_i, 0, -2*div_i],
                               [0, 0, div_i, 0],
                               ])

    else:
        if update_env:
            Msys = np.asarray([[div_i, 0, 0, 3*div_i],
                               [div_o, 0, 0, 3*div_o],
                               [0, div_i, 0, -2*div_i],
                               [0, div_o, 0, -2*div_o],
                               [0, 0, div_i, 0],
                               [0, 0, div_o, 0],
                               [p.F/cm, p.F/cm, -p.F/cm, p.F/cm]
                               ])
        else:
            Msys = np.asarray([[div_i, 0, 0, 3*div_i],
                               [0, div_i, 0, -2*div_i],
                               [0, 0, div_i, 0],
                               [p.F/cm, p.F/cm, -p.F/cm, p.F/cm]
                               ])

    return Msys