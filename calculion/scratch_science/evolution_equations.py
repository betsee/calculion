#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright (c) 2022-2023 Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Calculate updates to key parameters from implicit solutions of the system.
'''

# ....................{ IMPORTS                            }....................
from beartype import beartype
from beartype.typing import Union
import numpy as np
from numpy import ndarray
from calculion.scratch_science.params import CalculionParams
from calculion.scratch_science.vmem import vmem_ghk_pump

@beartype
def update_all_var(variables_list: Union[ndarray, list],
                   p: CalculionParams,
                   quasi_static_vmem: bool,
                   update_env: bool) -> tuple[list[float], list[float]]:
    '''
    Calculate updated values to V_mem and ion concentrations inside (and optionally outside) the cell.

    '''

    Na_i = variables_list[0]
    Na_o = variables_list[1]
    K_i = variables_list[2]
    K_o = variables_list[3]
    Cl_i = variables_list[4]
    Cl_o = variables_list[5]
    V_mem = variables_list[6]


    Na_o0 = Na_o
    Na_i0 = Na_i
    K_o0 = K_o
    K_i0 = K_i
    Cl_o0 = Cl_o
    Cl_i0 = Cl_i

    Delta_t = p.delta_t
    K_eqm_NaK = p.Keqm_NaK
    Omega_pump = p.omega_NaK
    ADP = p.ADP
    ATP = p.ATP
    P = p.P
    P_Na = p.P_Na
    P_K = p.P_K
    P_Cl = p.P_Cl
    z_Na = p.z_Na
    z_K = p.z_K
    z_Cl = p.z_Cl
    r_cell = p.r_cell
    d_ecm = p.d_ecm
    alpha = p.alpha
    F = p.F
    c_mem = p.c_mem

    V_mem0 = V_mem

    if quasi_static_vmem:
        numNai = (-6 * ADP * Cl_o0 * Delta_t * K_eqm_NaK * K_i0 ** 2 * Na_o0 ** 3 * Omega_pump * P * P_Cl * z_Cl +
                  6 * ADP * Delta_t * K_eqm_NaK * K_i0 ** 3 * Na_o0 ** 3 * Omega_pump * P * P_K * z_K +
                  4 * ADP * Delta_t * K_eqm_NaK * K_i0 ** 2 * Na_i0 * Na_o0 ** 3 * Omega_pump * P * P_Na * z_K +
                  6 * ATP * Cl_i0 * Delta_t * K_o0 ** 2 * Na_i0 ** 3 * Omega_pump * P_Cl * z_Cl -
                  6 * ATP * Delta_t * K_o0 ** 3 * Na_i0 ** 3 * Omega_pump * P_K * z_K -
                  4 * ATP * Delta_t * K_o0 ** 2 * Na_i0 ** 3 * Na_o0 * Omega_pump * P_Na * z_K -
                  2 * ATP * K_o0 ** 2 * Na_i0 ** 4 * Omega_pump * r_cell * z_K +
                  3 * ATP * K_o0 ** 2 * Na_i0 ** 4 * Omega_pump * r_cell * z_Na +
                  2 * Cl_i0 * Delta_t * Na_i0 * P_Cl * P_Na * z_Cl -
                  2 * Cl_o0 * Delta_t * Na_o0 * P_Cl * P_Na * z_Cl -
                  Cl_o0 * Na_i0 * P_Cl * r_cell * z_Cl +
                  2 * Delta_t * K_i0 * Na_o0 * P_K * P_Na * z_K -
                  2 * Delta_t * K_o0 * Na_i0 * P_K * P_Na * z_K +
                  K_i0 * Na_i0 * P_K * r_cell * z_K +
                  Na_i0 ** 2 * P_Na * r_cell * z_Na
                  )

        denNai = (r_cell * (-2 * ATP * K_o0 ** 2 * Na_i0 ** 3 * Omega_pump * z_K +
                            3 * ATP * K_o0 ** 2 * Na_i0 ** 3 * Omega_pump * z_Na -
                            Cl_o0 * P_Cl * z_Cl + K_i0 * P_K * z_K + Na_i0 * P_Na * z_Na)
                  )

        Nai1 = numNai / denNai

        # ----K in-------
        numKi = (4*ADP*Cl_o0*Delta_t*K_eqm_NaK*K_i0**2*Na_o0**3*Omega_pump*P*P_Cl*z_Cl -
                 6*ADP*Delta_t*K_eqm_NaK*K_i0**3*Na_o0**3*Omega_pump*P*P_K*z_Na -
                 4*ADP*Delta_t*K_eqm_NaK*K_i0**2*Na_i0*Na_o0**3*Omega_pump*P*P_Na*z_Na -
                 4*ATP*Cl_i0*Delta_t*K_o0**2*Na_i0**3*Omega_pump*P_Cl*z_Cl +
                 6*ATP*Delta_t*K_o0**3*Na_i0**3*Omega_pump*P_K*z_Na +
                 4*ATP*Delta_t*K_o0**2*Na_i0**3*Na_o0*Omega_pump*P_Na*z_Na -
                 2*ATP*K_i0*K_o0**2*Na_i0**3*Omega_pump*r_cell*z_K +
                 3*ATP*K_i0*K_o0**2*Na_i0**3*Omega_pump*r_cell*z_Na +
                 2*Cl_i0*Delta_t*K_i0*P_Cl*P_K*z_Cl -
                 2*Cl_o0*Delta_t*K_o0*P_Cl*P_K*z_Cl -
                 Cl_o0*K_i0*P_Cl*r_cell*z_Cl -
                 2*Delta_t*K_i0*Na_o0*P_K*P_Na*z_Na +
                 2*Delta_t*K_o0*Na_i0*P_K*P_Na*z_Na +
                 K_i0**2*P_K*r_cell*z_K +
                 K_i0*Na_i0*P_Na*r_cell*z_Na
)

        denKi = (r_cell*(-2*ATP*K_o0**2*Na_i0**3*Omega_pump*z_K +
                         3*ATP*K_o0**2*Na_i0**3*Omega_pump*z_Na -
                         Cl_o0*P_Cl*z_Cl + K_i0*P_K*z_K + Na_i0*P_Na*z_Na)
                 )

        Ki1 = numKi / denKi

        #---Cl in___
        numCli = (-4*ADP*Cl_o0*Delta_t*K_eqm_NaK*K_i0**2*Na_o0**3*Omega_pump*P*P_Cl*z_K +
                  6*ADP*Cl_o0*Delta_t*K_eqm_NaK*K_i0**2*Na_o0**3*Omega_pump*P*P_Cl*z_Na +
                  4*ATP*Cl_i0*Delta_t*K_o0**2*Na_i0**3*Omega_pump*P_Cl*z_K -
                  6*ATP*Cl_i0*Delta_t*K_o0**2*Na_i0**3*Omega_pump*P_Cl*z_Na -
                  2*ATP*Cl_i0*K_o0**2*Na_i0**3*Omega_pump*r_cell*z_K +
                  3*ATP*Cl_i0*K_o0**2*Na_i0**3*Omega_pump*r_cell*z_Na -
                  Cl_i0*Cl_o0*P_Cl*r_cell*z_Cl -
                  2*Cl_i0*Delta_t*K_i0*P_Cl*P_K*z_K -
                  2*Cl_i0*Delta_t*Na_i0*P_Cl*P_Na*z_Na +
                  Cl_i0*K_i0*P_K*r_cell*z_K +
                  Cl_i0*Na_i0*P_Na*r_cell*z_Na +
                  2*Cl_o0*Delta_t*K_o0*P_Cl*P_K*z_K +
                  2*Cl_o0*Delta_t*Na_o0*P_Cl*P_Na*z_Na)

        denCli = (r_cell*(-2*ATP*K_o0**2*Na_i0**3*Omega_pump*z_K +
                          3*ATP*K_o0**2*Na_i0**3*Omega_pump*z_Na -
                          Cl_o0*P_Cl*z_Cl + K_i0*P_K*z_K + Na_i0*P_Na*z_Na)
                    )

        Cli1 = numCli/denCli

        # ---Vmem -----
        numVm = np.log((-2*ADP*K_eqm_NaK*K_i0**2*Na_o0**3*Omega_pump*P*z_K +
                      3*ADP*K_eqm_NaK*K_i0**2*Na_o0**3*Omega_pump*P*z_Na -
                      Cl_i0*P_Cl*z_Cl +
                      K_o0*P_K*z_K +
                      Na_o0*P_Na*z_Na)/(-2*ATP*K_o0**2*Na_i0**3*Omega_pump*z_K +
                                        3*ATP*K_o0**2*Na_i0**3*Omega_pump*z_Na -
                                        Cl_o0*P_Cl*z_Cl +
                                        K_i0*P_K*z_K +
                                        Na_i0*P_Na*z_Na))


        denVm = (alpha)

        Vm1 = numVm/denVm

        if update_env:
            numNao = (6*ADP*Cl_o0*Delta_t*K_eqm_NaK*K_i0**2*Na_o0**3*Omega_pump*P*P_Cl*r_cell*z_Cl -
                      6*ADP*Delta_t*K_eqm_NaK*K_i0**3*Na_o0**3*Omega_pump*P*P_K*r_cell*z_K -
                      4*ADP*Delta_t*K_eqm_NaK*K_i0**2*Na_i0*Na_o0**3*Omega_pump*P*P_Na*r_cell*z_K -
                      6*ATP*Cl_i0*Delta_t*K_o0**2*Na_i0**3*Omega_pump*P_Cl*r_cell*z_Cl +
                      6*ATP*Delta_t*K_o0**3*Na_i0**3*Omega_pump*P_K*r_cell*z_K +
                      4*ATP*Delta_t*K_o0**2*Na_i0**3*Na_o0*Omega_pump*P_Na*r_cell*z_K -
                      2*ATP*K_o0**2*Na_i0**3*Na_o0*Omega_pump*d_ecm**2*z_K +
                      3*ATP*K_o0**2*Na_i0**3*Na_o0*Omega_pump*d_ecm**2*z_Na -
                      4*ATP*K_o0**2*Na_i0**3*Na_o0*Omega_pump*d_ecm*r_cell*z_K +
                      6*ATP*K_o0**2*Na_i0**3*Na_o0*Omega_pump*d_ecm*r_cell*z_Na -
                      2*Cl_i0*Delta_t*Na_i0*P_Cl*P_Na*r_cell*z_Cl +
                      2*Cl_o0*Delta_t*Na_o0*P_Cl*P_Na*r_cell*z_Cl -
                      Cl_o0*Na_o0*P_Cl*d_ecm**2*z_Cl -
                      2*Cl_o0*Na_o0*P_Cl*d_ecm*r_cell*z_Cl -
                      2*Delta_t*K_i0*Na_o0*P_K*P_Na*r_cell*z_K +
                      2*Delta_t*K_o0*Na_i0*P_K*P_Na*r_cell*z_K +
                      K_i0*Na_o0*P_K*d_ecm**2*z_K +
                      2*K_i0*Na_o0*P_K*d_ecm*r_cell*z_K +
                      Na_i0*Na_o0*P_Na*d_ecm**2*z_Na + 2*Na_i0*Na_o0*P_Na*d_ecm*r_cell*z_Na)

            denNao = (d_ecm*(-2*ATP*K_o0**2*Na_i0**3*Omega_pump*d_ecm*z_K +
                             3*ATP*K_o0**2*Na_i0**3*Omega_pump*d_ecm*z_Na -
                             4*ATP*K_o0**2*Na_i0**3*Omega_pump*r_cell*z_K +
                             6*ATP*K_o0**2*Na_i0**3*Omega_pump*r_cell*z_Na -
                             Cl_o0*P_Cl*d_ecm*z_Cl -
                             2*Cl_o0*P_Cl*r_cell*z_Cl +
                             K_i0*P_K*d_ecm*z_K +
                             2*K_i0*P_K*r_cell*z_K +
                             Na_i0*P_Na*d_ecm*z_Na +
                             2*Na_i0*P_Na*r_cell*z_Na))

            Nao1 = numNao/denNao

            # ----K out--------
            numKo = (-4*ADP*Cl_o0*Delta_t*K_eqm_NaK*K_i0**2*Na_o0**3*Omega_pump*P*P_Cl*r_cell*z_Cl +
                     6*ADP*Delta_t*K_eqm_NaK*K_i0**3*Na_o0**3*Omega_pump*P*P_K*r_cell*z_Na +
                     4*ADP*Delta_t*K_eqm_NaK*K_i0**2*Na_i0*Na_o0**3*Omega_pump*P*P_Na*r_cell*z_Na +
                     4*ATP*Cl_i0*Delta_t*K_o0**2*Na_i0**3*Omega_pump*P_Cl*r_cell*z_Cl -
                     6*ATP*Delta_t*K_o0**3*Na_i0**3*Omega_pump*P_K*r_cell*z_Na -
                     4*ATP*Delta_t*K_o0**2*Na_i0**3*Na_o0*Omega_pump*P_Na*r_cell*z_Na -
                     2*ATP*K_o0**3*Na_i0**3*Omega_pump*d_ecm**2*z_K +
                     3*ATP*K_o0**3*Na_i0**3*Omega_pump*d_ecm**2*z_Na -
                     4*ATP*K_o0**3*Na_i0**3*Omega_pump*d_ecm*r_cell*z_K +
                     6*ATP*K_o0**3*Na_i0**3*Omega_pump*d_ecm*r_cell*z_Na -
                     2*Cl_i0*Delta_t*K_i0*P_Cl*P_K*r_cell*z_Cl +
                     2*Cl_o0*Delta_t*K_o0*P_Cl*P_K*r_cell*z_Cl -
                     Cl_o0*K_o0*P_Cl*d_ecm**2*z_Cl -
                     2*Cl_o0*K_o0*P_Cl*d_ecm*r_cell*z_Cl +
                     2*Delta_t*K_i0*Na_o0*P_K*P_Na*r_cell*z_Na -
                     2*Delta_t*K_o0*Na_i0*P_K*P_Na*r_cell*z_Na +
                     K_i0*K_o0*P_K*d_ecm**2*z_K +
                     2*K_i0*K_o0*P_K*d_ecm*r_cell*z_K +
                     K_o0*Na_i0*P_Na*d_ecm**2*z_Na +
                     2*K_o0*Na_i0*P_Na*d_ecm*r_cell*z_Na)

            denKo = (d_ecm*(-2*ATP*K_o0**2*Na_i0**3*Omega_pump*d_ecm*z_K +
                            3*ATP*K_o0**2*Na_i0**3*Omega_pump*d_ecm*z_Na -
                            4*ATP*K_o0**2*Na_i0**3*Omega_pump*r_cell*z_K +
                            6*ATP*K_o0**2*Na_i0**3*Omega_pump*r_cell*z_Na -
                            Cl_o0*P_Cl*d_ecm*z_Cl -
                            2*Cl_o0*P_Cl*r_cell*z_Cl +
                            K_i0*P_K*d_ecm*z_K +
                            2*K_i0*P_K*r_cell*z_K +
                            Na_i0*P_Na*d_ecm*z_Na +
                            2*Na_i0*P_Na*r_cell*z_Na))

            Ko1 = numKo/denKo

            # ---Cl out___
            numClo = (4*ADP*Cl_o0*Delta_t*K_eqm_NaK*K_i0**2*Na_o0**3*Omega_pump*P*P_Cl*r_cell*z_K -
                      6*ADP*Cl_o0*Delta_t*K_eqm_NaK*K_i0**2*Na_o0**3*Omega_pump*P*P_Cl*r_cell*z_Na -
                      4*ATP*Cl_i0*Delta_t*K_o0**2*Na_i0**3*Omega_pump*P_Cl*r_cell*z_K +
                      6*ATP*Cl_i0*Delta_t*K_o0**2*Na_i0**3*Omega_pump*P_Cl*r_cell*z_Na -
                      2*ATP*Cl_o0*K_o0**2*Na_i0**3*Omega_pump*d_ecm**2*z_K +
                      3*ATP*Cl_o0*K_o0**2*Na_i0**3*Omega_pump*d_ecm**2*z_Na -
                      4*ATP*Cl_o0*K_o0**2*Na_i0**3*Omega_pump*d_ecm*r_cell*z_K +
                      6*ATP*Cl_o0*K_o0**2*Na_i0**3*Omega_pump*d_ecm*r_cell*z_Na +
                      2*Cl_i0*Delta_t*K_i0*P_Cl*P_K*r_cell*z_K +
                      2*Cl_i0*Delta_t*Na_i0*P_Cl*P_Na*r_cell*z_Na -
                      Cl_o0**2*P_Cl*d_ecm**2*z_Cl -
                      2*Cl_o0**2*P_Cl*d_ecm*r_cell*z_Cl -
                      2*Cl_o0*Delta_t*K_o0*P_Cl*P_K*r_cell*z_K -
                      2*Cl_o0*Delta_t*Na_o0*P_Cl*P_Na*r_cell*z_Na +
                      Cl_o0*K_i0*P_K*d_ecm**2*z_K +
                      2*Cl_o0*K_i0*P_K*d_ecm*r_cell*z_K +
                      Cl_o0*Na_i0*P_Na*d_ecm**2*z_Na +
                      2*Cl_o0*Na_i0*P_Na*d_ecm*r_cell*z_Na)

            denClo = (d_ecm*(-2*ATP*K_o0**2*Na_i0**3*Omega_pump*d_ecm*z_K +
                             3*ATP*K_o0**2*Na_i0**3*Omega_pump*d_ecm*z_Na -
                             4*ATP*K_o0**2*Na_i0**3*Omega_pump*r_cell*z_K +
                             6*ATP*K_o0**2*Na_i0**3*Omega_pump*r_cell*z_Na -
                             Cl_o0*P_Cl*d_ecm*z_Cl -
                             2*Cl_o0*P_Cl*r_cell*z_Cl +
                             K_i0*P_K*d_ecm*z_K +
                             2*K_i0*P_K*r_cell*z_K +
                             Na_i0*P_Na*d_ecm*z_Na +
                             2*Na_i0*P_Na*r_cell*z_Na)
                        )

            Clo1 = numClo / denClo

        else: # They're equal to the same as the inputs (no change)
            Nao1 = Na_o0
            Ko1 = K_o0
            Clo1 = Cl_o0

    else:
        # ---Na in___
        numNai = (6*ADP*Delta_t*K_eqm_NaK*K_i0**2*Na_o0**3*Omega_pump*P -
                  6*ATP*Delta_t*K_o0**2*Na_i0**3*Omega_pump*np.exp(V_mem0*alpha) -
                  2*Delta_t*Na_i0*P_Na*np.exp(V_mem0*alpha) +
                  2*Delta_t*Na_o0*P_Na + Na_i0*r_cell)

        denNai = (r_cell)

        Nai1 = numNai / denNai

        # ---K in___
        numKi = (-4*ADP*Delta_t*K_eqm_NaK*K_i0**2*Na_o0**3*Omega_pump*P +
                 4*ATP*Delta_t*K_o0**2*Na_i0**3*Omega_pump*np.exp(V_mem0*alpha) -
                 2*Delta_t*K_i0*P_K*np.exp(V_mem0*alpha) +
                 2*Delta_t*K_o0*P_K + K_i0*r_cell)

        denKi = (r_cell)

        Ki1 = numKi / denKi

        # ---Cl in___
        numCli = (-2*Cl_i0*Delta_t*P_Cl + Cl_i0*r_cell + 2*Cl_o0*Delta_t*P_Cl*np.exp(V_mem0*alpha))

        denCli = (r_cell)

        Cli1 = numCli / denCli

        # ---Vmem -----
        Vm1 = (-2*ADP*Delta_t*F*K_eqm_NaK*K_i0**2*Na_o0**3*Omega_pump*P*z_K/c_mem +
               3*ADP*Delta_t*F*K_eqm_NaK*K_i0**2*Na_o0**3*Omega_pump*P*z_Na/c_mem +
               2*ATP*Delta_t*F*K_o0**2*Na_i0**3*Omega_pump*z_K*np.exp(V_mem0*alpha)/c_mem -
               3*ATP*Delta_t*F*K_o0**2*Na_i0**3*Omega_pump*z_Na*np.exp(V_mem0*alpha)/c_mem -
               Cl_i0*Delta_t*F*P_Cl*z_Cl/c_mem + Cl_o0*Delta_t*F*P_Cl*z_Cl*np.exp(V_mem0*alpha)/c_mem -
               Delta_t*F*K_i0*P_K*z_K*np.exp(V_mem0*alpha)/c_mem + Delta_t*F*K_o0*P_K*z_K/c_mem -
               Delta_t*F*Na_i0*P_Na*z_Na*np.exp(V_mem0*alpha)/c_mem + Delta_t*F*Na_o0*P_Na*z_Na/c_mem + V_mem0)

        if update_env:
            # ---Na out___
            numNao = (-6*ADP*Delta_t*K_eqm_NaK*K_i0**2*Na_o0**3*Omega_pump*P*r_cell +
                      6*ATP*Delta_t*K_o0**2*Na_i0**3*Omega_pump*r_cell*np.exp(V_mem0*alpha) +
                      2*Delta_t*Na_i0*P_Na*r_cell*np.exp(V_mem0*alpha) -
                      2*Delta_t*Na_o0*P_Na*r_cell + Na_o0*d_ecm**2 + 2*Na_o0*d_ecm*r_cell)

            denNao = (d_ecm**2 + 2*d_ecm*r_cell)

            Nao1 = numNao / denNao

            # ---K out___
            numKo = (4*ADP*Delta_t*K_eqm_NaK*K_i0**2*Na_o0**3*Omega_pump*P*r_cell -
                     4*ATP*Delta_t*K_o0**2*Na_i0**3*Omega_pump*r_cell*np.exp(V_mem0*alpha) +
                     2*Delta_t*K_i0*P_K*r_cell*np.exp(V_mem0*alpha) -
                     2*Delta_t*K_o0*P_K*r_cell + K_o0*d_ecm**2 + 2*K_o0*d_ecm*r_cell)

            denKo = (d_ecm**2 + 2*d_ecm*r_cell)

            Ko1 = numKo / denKo

            # ---Cl out___
            numClo = (2*Cl_i0*Delta_t*P_Cl*r_cell -
                      2*Cl_o0*Delta_t*P_Cl*r_cell*np.exp(V_mem0*alpha) +
                      Cl_o0*d_ecm**2 + 2*Cl_o0*d_ecm*r_cell
)

            denClo = (d_ecm**2 + 2*d_ecm*r_cell)

            Clo1 = numClo / denClo
        else:
            Nao1 = Na_o0
            Ko1 = K_o0
            Clo1 = Cl_o0

    delNai = Nai1 - Na_i
    delKi = Ki1 - K_i
    delCli = Cli1 - Cl_i
    delVm = Vm1 - V_mem

    if update_env:
        delNao = Nao1 - Na_o
        delKo = Ko1 - K_o
        delClo = Clo1 - Cl_o

        change_list = [delNai, delNao, delKi, delKo, delCli, delClo, delVm]

    else:
        change_list = [delNai, delKi, delCli, delVm]

    return [Nai1, Nao1, Ki1, Ko1, Cli1, Clo1, Vm1], change_list

def unit_loop(var_list: list[float],
              p: CalculionParams,
              quasi_static_vmem: bool,
              update_env: bool):
    '''
    Unit loop based on a var_list of [Na_i, Na_o, K_i, K_o, Cl_i, Cl_o, V_mem].
    Returns the sum of squared parameter changes to each variable in the list.

    '''
    # compute change to the variables:
    _, change_list = update_all_var(var_list,
                                           p,
                                           quasi_static_vmem=quasi_static_vmem,
                                           update_env=update_env)
    # Return the sum of squares:
    return np.sum(np.asarray(change_list)**2)

# FIXME: the following don't really work, so commenting them out for now
# def unit_loop_inverse_all(var_list: list[float],
#               p: CalculionParams,
#               quasi_static_vmem: bool,
#               update_env: bool):
#     '''
#     Unit loop based on a var_list of var_list = [P_Na, P_K, P_Cl, omega_NaK, Keqm_NaK, ATP, ADP, P], to
#     assist in solving the inverse-optimization problem.
#     Returns the sum of squared parameter changes to system steady state.
#
#     '''
#
#     # Change values of parameters according to the input variables list:
#     p.P_Na = var_list[0]
#     p.P_K = var_list[1]
#     p.P_Cl = var_list[2]
#     p.omega_NaK = var_list[3]
#     p.Keqm_NaK = var_list[4]
#     p.ATP = var_list[5]
#     p.ADP = var_list[6]
#     p.P = var_list[7]
#
#     param_list = [p.Na_i, p.Na_o, p.K_i, p.K_o, p.Cl_i, p.Cl_o, p.target_Vmem]
#     # compute change to the variables:
#     _, change_list = update_all_var(param_list,
#                                            p,
#                                            quasi_static_vmem=quasi_static_vmem,
#                                            update_env=update_env)
#     # Return the sum of squares:
#     return np.sum(np.asarray(change_list)**2)
#
# def unit_loop_inverse_mem(var_list: list[float],
#               p: CalculionParams,
#               quasi_static_vmem: bool,
#               update_env: bool):
#     '''
#     Unit loop based on a var_list of var_list = [P_Na, P_K, P_Cl], to
#     assist in solving the inverse-optimization problem focusing on membrane permeabilities.
#     Returns the sum of squared parameter changes to system steady state.
#
#     '''
#
#     # Change values of parameters according to the input variables list:
#     p.P_Na = var_list[0]
#     p.P_K = var_list[1]
#     p.P_Cl = var_list[2]
#     # p.omega_NaK = var_list[3]
#     # p.Keqm_NaK = var_list[4]
#     # p.ATP = var_list[5]
#     # p.ADP = var_list[6]
#     # p.P = var_list[7]
#
#     param_list = [p.Na_i, p.Na_o, p.K_i, p.K_o, p.Cl_i, p.Cl_o, p.target_Vmem]
#     # compute change to the variables:
#     _, change_list = update_all_var(param_list,
#                                            p,
#                                            quasi_static_vmem=quasi_static_vmem,
#                                            update_env=update_env)
#     # Return the sum of squares:
#     return np.sum(np.asarray(change_list)**2)
#
# def unit_loop_inverse_pump(var_list: list[float],
#               p: CalculionParams,
#               quasi_static_vmem: bool,
#               update_env: bool):
#     '''
#     Unit loop based on a var_list of var_list = [omega_NaK, Keqm_NaK], to
#     assist in solving the inverse-optimization problem focusing on NaK-ATPase pump properties.
#     Returns the sum of squared parameter changes to system steady state.
#
#     '''
#
#     # Change values of parameters according to the input variables list:
#     p.omega_NaK = var_list[0]
#     p.Keqm_NaK = var_list[1]
#
#     param_list = [p.Na_i, p.Na_o, p.K_i, p.K_o, p.Cl_i, p.Cl_o, p.target_Vmem]
#     # compute change to the variables:
#     _, change_list = update_all_var(param_list,
#                                            p,
#                                            quasi_static_vmem=quasi_static_vmem,
#                                            update_env=update_env)
#     # Return the sum of squares:
#     return np.sum(np.asarray(change_list)**2)
#
# def unit_loop_inverse_met(var_list: list[float],
#               p: CalculionParams,
#               quasi_static_vmem: bool,
#               update_env: bool):
#     '''
#     Unit loop based on a var_list of var_list = [omega_NaK, Keqm_NaK, ATP, ADP, P], to
#     assist in solving the inverse-optimization problem focusing on NaK-ATPase and metabolic parameters.
#     Returns the sum of squared parameter changes to system steady state.
#
#     '''
#
#     # Change values of parameters according to the input variables list:
#     p.omega_NaK = var_list[0]
#     p.Keqm_NaK = var_list[1]
#     p.ATP = var_list[2]
#     p.ADP = var_list[3]
#     p.P = var_list[4]
#
#     param_list = [p.Na_i, p.Na_o, p.K_i, p.K_o, p.Cl_i, p.Cl_o, p.target_Vmem]
#     # compute change to the variables:
#     _, change_list = update_all_var(param_list,
#                                            p,
#                                            quasi_static_vmem=quasi_static_vmem,
#                                            update_env=update_env)
#     # Return the sum of squares:
#     return np.sum(np.asarray(change_list)**2)