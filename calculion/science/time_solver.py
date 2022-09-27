#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright (c) 2022-2023 Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Class to generate temporal updates to transmembrane potential Vmem and ion concentrations inside and out of the cell.
'''

# ....................{ IMPORTS                            }....................
from beartype import beartype
import numpy as np
import sympy as sp
from calculion.science.sci_enum import TimeSolverType


class TimeSolver(object):
    '''

    '''

    def __int__(self):
        '''

        '''
        pass


    def _compute_update_eq(self, solver_type: TimeSolverType=TimeSolverType.implicit):
        '''

        '''

        # Define symbols:
        V_mem, z_Na, z_K, z_Cl = sp.symbols('V_mem, z_Na, z_K, z_Cl', real=True)

        d_ecm, r_cell, c_mem, delta_t, F, R, T, P_Na, P_K, P_Cl = sp.symbols('d_ecm, r_cell, c_mem, delta_t, '
                                                                             'F, R, T, '
                                                                             'P_Na, P_K, P_Cl',
                                                                             real=True, positive=True)

        alpha = F/(R*T)

        Na_o0, Na_i0, Na_o1, Na_i1 = sp.symbols('Na_o0, Na_i0, Na_o1, Na_i1', real=True, positive=True)
        K_o0, K_i0, K_o1, K_i1 = sp.symbols('K_o0, K_i0, K_o1, K_i1', real=True, positive=True)
        Cl_o0, Cl_i0, Cl_o1, Cl_i1 = sp.symbols('Cl_o0, Cl_i0, Cl_o1, Cl_i1', real=True, positive=True)
        V_mem0, V_mem1 = sp.symbols('V_mem0, V_mem1', real=True)

        # Divergence operators (assuming constant flux across the membrane):
        div_i = (2/r_cell) # divergence from membrane flux into cell
        div_o = ((2*r_cell)/(2*r_cell*d_ecm + d_ecm**2)) # divergence from membrane flux into extracellular space


        if solver_type is TimeSolverType.explicit:
            f_Na0 = P_Na*(Na_o0 - Na_i0 * sp.exp(alpha * V_mem0))  # Na flux in terms of past Vmem, past Na in and out
            f_K0 = P_Na*(K_o0 - K_i0 * sp.exp(alpha * V_mem0))  # K flux in terms of past Vmem, past K in and out
            f_Cl0 = -P_Cl*(Cl_i0 - Cl_o0 * sp.exp(alpha * V_mem0)) # Cl flux in terms of past Vmem, past Cl in and out

            # Conduction current:
            J_cond0 = z_Na*F*f_Na0 + z_K*F*f_K0 + z_Cl*F*f_Cl0  # current in terms of past concs. and past Vmem

            # Concentration updates in time, inside and out of the cell:
            eq_Nai_t0 = sp.Eq((Na_i1 - Na_i0) / delta_t, div_i * f_Na0)
            eq_Nao_t0 = sp.Eq((Na_o1 - Na_o0) / delta_t, -div_o * f_Na0)

            eq_Ki_t0 = sp.Eq((K_i1 - K_i0) / delta_t, div_i * f_K0)
            eq_Ko_t0 = sp.Eq((K_o1 - K_o0) / delta_t, -div_o * f_K0)

            eq_Cli_t0 = sp.Eq((Cl_i1 - Cl_i0) / delta_t, div_i * f_Cl0)
            eq_Clo_t0 = sp.Eq((Cl_o1 - Cl_o0) / delta_t, -div_o * f_Cl0)

            # Voltage update in time:
            eq_Vm_t0 = sp.Eq((V_mem1 - V_mem0) / delta_t, (1 / c_mem) * J_cond0)

            # Fully explicit time stepping: solve concentrations in terms of past Vmem and past concentrations;
            # solve Vmem in terms of past concentrations:
            sol0 = sp.solve((eq_Nai_t0, eq_Nao_t0, eq_Ki_t0, eq_Ko_t0, eq_Cli_t0, eq_Clo_t0, eq_Vm_t0),
                       (Na_i1, Na_o1, K_i1, K_o1, Cl_i1, Cl_o1, V_mem1)
                       )

            sol0_list = []

            sol0sub = []
            for k, v in sol0.items():
                sol0sub.append(sol0[k].simplify())

            sol0_list.append(sol0sub)

        elif solver_type is TimeSolverType.implicit:
            # Flux equations:
            f_Na01 = P_Na*(Na_o0 - Na_i0*sp.exp(alpha*V_mem1))  # Na flux in terms of future Vmem, past Na in and out
            f_K01 = P_Na*(K_o0 - K_i0 * sp.exp(alpha*V_mem1))  # K flux in terms of future Vmem, past K in and out
            f_Cl01 = -P_Cl*(Cl_i0 - Cl_o0 * sp.exp(alpha*V_mem1))  # Cl flux in terms of future Vmem, past Cl in and out

            # Conduction current:
            J_cond01 = z_Na*F*f_Na01 + z_K*F*f_K01 + z_Cl*F*f_Cl01  # current in terms of past concs. and future Vmem

            # Concentration updates in time, inside and out of the cell:
            eq_Nai_t01 = sp.Eq((Na_i1 - Na_i0) / delta_t, div_i * f_Na01)
            eq_Nao_t01 = sp.Eq((Na_o1 - Na_o0) / delta_t, -div_o * f_Na01)

            eq_Ki_t01 = sp.Eq((K_i1 - K_i0) / delta_t, div_i * f_K01)
            eq_Ko_t01 = sp.Eq((K_o1 - K_o0) / delta_t, -div_o * f_K01)

            eq_Cli_t01 = sp.Eq((Cl_i1 - Cl_i0) / delta_t, div_i * f_Cl01)
            eq_Clo_t01 = sp.Eq((Cl_o1 - Cl_o0) / delta_t, -div_o * f_Cl01)

            # Voltage update in time:
            eq_Vm_t01 = sp.Eq((V_mem1 - V_mem0) / delta_t, (1 / c_mem) * J_cond01)

            # Partial implicit: solve concentrations in terms of past concentration and future Vmem;
            # solve Vmem in terms of past concentrations and future Vmem:
            sol0_list = sp.solve((eq_Nai_t01, eq_Nao_t01, eq_Ki_t01, eq_Ko_t01, eq_Cli_t01, eq_Clo_t01, eq_Vm_t01),
                       (Na_i1, Na_o1, K_i1, K_o1, Cl_i1, Cl_o1, V_mem1)
                       )

        else:
            raise Exception("Only implicit and explicit time solver types are supported.")


        # Create Lambdified functions for numerical computations:
        # Future Na+ inside cell:
        Na_i1_funk = sp.lambdify(((F, R, T, d_ecm, r_cell, c_mem, z_Na, z_K, z_Cl, delta_t),
                                  (P_Na, P_K, P_Cl), (Na_o0, Na_i0, K_o0, K_i0, Cl_o0, Cl_i0, V_mem0)),
                                 sol0_list[0][0].simplify())
        # Future Na+ outside cell:
        Na_o1_funk = sp.lambdify(((F, R, T, d_ecm, r_cell, c_mem, z_Na, z_K, z_Cl, delta_t),
                                  (P_Na, P_K, P_Cl), (Na_o0, Na_i0, K_o0, K_i0, Cl_o0, Cl_i0, V_mem0)),
                                 sol0_list[0][1].simplify())
        # Future K+ inside cell:
        K_i1_funk = sp.lambdify(((F, R, T, d_ecm, r_cell, c_mem, z_Na, z_K, z_Cl, delta_t),
                                  (P_Na, P_K, P_Cl), (Na_o0, Na_i0, K_o0, K_i0, Cl_o0, Cl_i0, V_mem0)),
                                 sol0_list[0][2].simplify())
        # Future K+ outside cell:
        K_o1_funk = sp.lambdify(((F, R, T, d_ecm, r_cell, c_mem, z_Na, z_K, z_Cl, delta_t),
                                  (P_Na, P_K, P_Cl), (Na_o0, Na_i0, K_o0, K_i0, Cl_o0, Cl_i0, V_mem0)),
                                 sol0_list[0][3].simplify())
        # Future Cl- inside cell:
        Cl_i1_funk = sp.lambdify(((F, R, T, d_ecm, r_cell, c_mem, z_Na, z_K, z_Cl, delta_t),
                                  (P_Na, P_K, P_Cl), (Na_o0, Na_i0, K_o0, K_i0, Cl_o0, Cl_i0, V_mem0)),
                                 sol0_list[0][4].simplify())
        # Future Cl- outside cell:
        Cl_o1_funk = sp.lambdify(((F, R, T, d_ecm, r_cell, c_mem, z_Na, z_K, z_Cl, delta_t),
                                  (P_Na, P_K, P_Cl), (Na_o0, Na_i0, K_o0, K_i0, Cl_o0, Cl_i0, V_mem0)),
                                 sol0_list[0][5].simplify())
        # Future Vmem:
        Vmem_funk = sp.lambdify(((F, R, T, d_ecm, r_cell, c_mem, z_Na, z_K, z_Cl, delta_t),
                                  (P_Na, P_K, P_Cl), (Na_o0, Na_i0, K_o0, K_i0, Cl_o0, Cl_i0, V_mem0)),
                                 sol0_list[0][6].simplify())

        param_update_funk_list = [Na_i1_funk, Na_o1_funk,
                                  K_i1_funk, K_o1_funk,
                                  Cl_i1_funk, Cl_o1_funk,
                                  Vmem_funk]

        return param_update_funk_list

    @beartype
    def time_loop(self, dt: float, N_iter: int=100):
        '''

        '''


# Voltage, measurements, constants:
# V_mem = psi_i - psi_o

# Concentrations inside and out of the cell, membrane permeabilities:
# Na_o, Na_i, z_Na, P_Na = sp.symbols('Na_o, Na_i, z_Na, P_Na')
# K_o, K_i, z_K, P_K = sp.symbols('K_o, K_i, z_K, P_K')
# Cl_o, Cl_i, z_Cl, P_Cl = sp.symbols('Cl_o, Cl_i, z_Cl, P_Cl')
# # X_o, X_i, z_X, P_X = sp.symbols('X_o, X_i, z_X, P_X') # Proteins

# Na_o, Na_i, K_o, K_i, Cl_o, Cl_i, ATP, ADP, P = sp.symbols('Na_o, Na_i, K_o, K_i, Cl_o, Cl_i, ATP, ADP, P',
#                                               real=True, positive=True)



# X_o, X_i, z_X, P_X = sp.symbols('X_o, X_i, z_X, P_X') # Proteins

# ATP concentrations, free energies of hydrolysis, rates of Na-K-ATPase and metabolism:
# K_eqm_NaK, omega_pump, omega_met = sp.symbols('K_eqm_NaK, Omega_pump, Omega_met', real=True)

# Individual ion fluxes:
# Cation fluxes:
# f_Na = P_Na*(Na_o - Na_i*sp.exp(alpha*V_mem)) # Na_o <--> Na_i  delGo_Na = 0
# f_K = P_K*(K_o - K_i*sp.exp(alpha*V_mem)) # K_o <--> K_i  delGo_K = 0
# Cl (anion) reaction is written backwards to ge the same exp(alpha*V_mem) term in the final expression
# f_Cl = -P_Cl*(Cl_i - Cl_o*sp.exp(alpha*V_mem)) # Cl_i <--> Cl_o  delGo_K = 0
# f_X = -P_X*(X_i - X_o*sp.exp(alpha*V_mem)) # X_i <--> X_o  delGo_K = 0

# Conduction current:
# J_cond = z_Na*F*f_Na + z_K*F*f_K + z_Cl*F*f_Cl
# J_cond = z_Na*F*f_Na + z_K*F*f_K + z_Cl*F*f_Cl + z_X*F*f_X

# # Solve for V_mem when J_cond = 0 yeilds the GHK equation:
# eq_GHK = sp.solve(J_cond, V_mem)[0].simplify()