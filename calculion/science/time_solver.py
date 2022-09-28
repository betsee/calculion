#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright (c) 2022-2023 Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Class to generate temporal updates to transmembrane potential Vmem and ion concentrations inside and out of the cell.
'''

# ....................{ IMPORTS                            }....................
from beartype import beartype
from beartype.typing import Optional
import numpy as np
import sympy as sp
from calculion.science.sci_enum import TimeSolverType
from calculion.science.params import CalculionParams
from calculion.science.vmem import vmem_ghk_pump
from calculion.science.flux import make_reaction_matrix, compute_flux_vector
from calculion.science.evolution_equations import update_all_var


class TimeSolver(object):
    '''

    '''

    def __init__(self, quasi_static_vmem: bool=True, update_env: bool=True):
        '''

        '''



        self._quasi_static_vmem = quasi_static_vmem
        self._update_env = update_env

        self.update_eqn_dict = None
        self.Msys = None

        # Create indices to access variables in parameter lists:
        self.i_Nai = 0
        self.i_Nao = 1
        self.i_Ki = 2
        self.i_Ko = 3
        self.i_Cli = 4
        self.i_Clo = 5
        self.i_Vmem = 6


    @beartype
    def _get_update_eqns(self, quasi_static_vmem: bool=True, update_env: bool=True):
        '''

        '''

        update_eqn_dict = {}

        V_mem, alpha, d_ecm, r_cell, c_mem, F, R, T = sp.symbols('V_mem, alpha, d_ecm, r_cell, c_mem, F, R, T',
                                                                 real=True)

        r_env = r_cell + d_ecm

        Na_o, Na_i, K_o, K_i, Cl_o, Cl_i, ATP, ADP, P = sp.symbols('Na_o, Na_i, K_o, K_i, Cl_o, Cl_i, ATP, ADP, P',
                                                                   real=True, positive=True)

        z_Na, z_K, z_Cl, P_Na, P_K, P_Cl = sp.symbols('z_Na, z_K, z_Cl, P_Na, P_K, P_Cl', real=True)

        Na_o0, Na_i0, Na_o1, Na_i1 = sp.symbols('Na_o0, Na_i0, Na_o1, Na_i1', real=True, positive=True)
        K_o0, K_i0, K_o1, K_i1 = sp.symbols('K_o0, K_i0, K_o1, K_i1', real=True, positive=True)
        Cl_o0, Cl_i0, Cl_o1, Cl_i1 = sp.symbols('Cl_o0, Cl_i0, Cl_o1, Cl_i1', real=True, positive=True)
        V_mem0, V_mem1 = sp.symbols('V_mem0, V_mem1', real=True)

        # ATP concentrations, free energies of hydrolysis, rates of Na-K-ATPase and metabolism:
        (K_eqm_NaK,
         omega_pump,
         omega_met) = sp.symbols('K_eqm_NaK, Omega_pump, Omega_met', real=True)

        # Divergence operator and other symbols:
        delta_t = sp.symbols('Delta_t', real=True, positive=True)

        # Divergence operators (assuming constant flux across the membrane):
        div_i = (2 / r_cell)
        div_o = -((2 * r_cell) / (r_env ** 2 - r_cell ** 2))

        # Flux equations:
        # Explicit (in terms of past values):
        f_Na_exp = P_Na * (Na_o0 - Na_i0 * sp.exp(alpha * V_mem))
        f_K_exp = P_K * (K_o0 - K_i0 * sp.exp(alpha * V_mem))
        f_Cl_exp = -P_Cl * (Cl_i0 - Cl_o0 * sp.exp(alpha * V_mem))  # Cl_i <--> Cl_o  delGo_K = 0

        # Define symbols:
        f_NaKpump_exp = (omega_pump * (K_eqm_NaK * ADP * P * (K_i0 ** 2) * (Na_o0 ** 3) -
                                        ATP * (K_o0 ** 2) * (Na_i0 ** 3) * sp.exp(alpha * V_mem))
                         )

        f_NaKpump_Na_exp = 3 * f_NaKpump_exp
        f_NaKpump_K_exp = -2 * f_NaKpump_exp

        # Add pump contribution to original ion fluxes:
        f_Na_pump_exp = f_Na_exp + f_NaKpump_Na_exp
        f_K_pump_exp = f_K_exp + f_NaKpump_K_exp
        f_Cl_pump_exp = f_Cl_exp

        J_cond_pump_exp = (z_Na * F * f_Na_pump_exp +
                           z_K * F * f_K_pump_exp +
                           z_Cl * F * f_Cl_pump_exp
                           )

        # Concentration updates in time, inside and out of the cell:
        eq_Nai_t_pump_exp = sp.Eq((Na_i1 - Na_i0) / delta_t, div_i * f_Na_pump_exp)
        eq_Nao_t_pump_exp = sp.Eq((Na_o1 - Na_o0) / delta_t, div_o * f_Na_pump_exp)

        eq_Ki_t_pump_exp = sp.Eq((K_i1 - K_i0) / delta_t, div_i * f_K_pump_exp)
        eq_Ko_t_pump_exp = sp.Eq((K_o1 - K_o0) / delta_t, div_o * f_K_pump_exp)

        eq_Cli_t_pump_exp = sp.Eq((Cl_i1 - Cl_i0) / delta_t, div_i * f_Cl_exp)
        eq_Clo_t_pump_exp = sp.Eq((Cl_o1 - Cl_o0) / delta_t, div_o * f_Cl_exp)

        # Here we use the quasi-static condition of zero membrane current instead of updating Vmem:
        eq_J_pump_exp = sp.Eq(J_cond_pump_exp, 0)

        sym_param_list = [Na_i0, Na_o0, K_i0, K_o0, Cl_i0, Cl_o0, V_mem0,
                          ATP, ADP, P, K_eqm_NaK, omega_pump,
                          P_Na, P_K, P_Cl, z_Na, z_K, z_Cl, alpha, F,
                          r_cell, d_ecm, c_mem, delta_t]

        if quasi_static_vmem:

            if update_env:

                # Solve the system of equations explicitly using the J=0 approximation:
                sol_pump_exp = sp.solve((eq_Nai_t_pump_exp,
                                         eq_Nao_t_pump_exp,
                                         eq_Ki_t_pump_exp,
                                         eq_Ko_t_pump_exp,
                                         eq_Cli_t_pump_exp,
                                         eq_Clo_t_pump_exp,
                                         eq_J_pump_exp),
                                        (Na_i1, Na_o1, K_i1, K_o1, Cl_i1, Cl_o1, V_mem)
                                        )

                # Lambdify all the resulting expressions:
                Nai_pump_funk = sp.lambdify((sym_param_list),
                                            sol_pump_exp[0][0])

                Nao_pump_funk = sp.lambdify((sym_param_list),
                                            sol_pump_exp[0][1])

                Ki_pump_funk = sp.lambdify((sym_param_list),
                                           sol_pump_exp[0][2])

                Ko_pump_funk = sp.lambdify((sym_param_list),
                                           sol_pump_exp[0][3])

                Cli_pump_funk = sp.lambdify((sym_param_list),
                                            sol_pump_exp[0][4])

                Clo_pump_funk = sp.lambdify((sym_param_list),
                                            sol_pump_exp[0][5])

                Vmem_pump_funk = sp.lambdify((sym_param_list),
                                             sol_pump_exp[0][6])

                update_eqn_dict['Nao'] = Nao_pump_funk
                update_eqn_dict['Ko'] = Ko_pump_funk
                update_eqn_dict['Clo'] = Clo_pump_funk

            else:
                # Solve the system of equations explicitly using the J=0 approximation but ignoring env updates:
                sol_pump_exp = sp.solve((eq_Nai_t_pump_exp,
                                         eq_Ki_t_pump_exp,
                                         eq_Cli_t_pump_exp,
                                         eq_J_pump_exp),
                                        (Na_i1, K_i1, Cl_i1, V_mem)
                                        )

                Nai_pump_funk = sp.lambdify((sym_param_list),
                                            sol_pump_exp[0][0])

                Ki_pump_funk = sp.lambdify((sym_param_list),
                                           sol_pump_exp[0][1])

                Cli_pump_funk = sp.lambdify((sym_param_list),
                                            sol_pump_exp[0][2])

                Vmem_pump_funk = sp.lambdify((sym_param_list),
                                             sol_pump_exp[0][3])

        else:

            # Fully explicit: solve concentrations in terms of past Vmem and past concentrations; solve Vmem
            # in terms of past concentrations; This does not use J=0 approximatino at mem:
            f_Na_exp2 = P_Na * (Na_o0 - Na_i0 * sp.exp(alpha * V_mem0))
            f_K_exp2 = P_K * (K_o0 - K_i0 * sp.exp(alpha * V_mem0))
            f_Cl_exp2 = -P_Cl * (Cl_i0 - Cl_o0 * sp.exp(alpha * V_mem0))  # Cl_i <--> Cl_o  delGo_K = 0

            f_NaKpump_exp2 = (omega_pump * (K_eqm_NaK * ADP * P * (K_i0 ** 2) * (Na_o0 ** 3) -
                                             ATP * (K_o0 ** 2) * (Na_i0 ** 3) * sp.exp(alpha * V_mem0))
                              )

            f_NaKpump_Na_exp2 = 3 * f_NaKpump_exp2
            f_NaKpump_K_exp2 = -2 * f_NaKpump_exp2

            f_Na_pump_exp2 = f_Na_exp2 + f_NaKpump_Na_exp2
            f_K_pump_exp2 = f_K_exp2 + f_NaKpump_K_exp2
            f_Cl_pump_exp2 = f_Cl_exp2  # Cl_i <--> Cl_o  delGo_K = 0

            eq_Nai_t_pump_exp2 = sp.Eq((Na_i1 - Na_i0) / delta_t, div_i * f_Na_pump_exp2)
            eq_Nao_t_pump_exp2 = sp.Eq((Na_o1 - Na_o0) / delta_t, div_o * f_Na_pump_exp2)

            eq_Ki_t_pump_exp2 = sp.Eq((K_i1 - K_i0) / delta_t, div_i * f_K_pump_exp2)
            eq_Ko_t_pump_exp2 = sp.Eq((K_o1 - K_o0) / delta_t, div_o * f_K_pump_exp2)

            eq_Cli_t_pump_exp2 = sp.Eq((Cl_i1 - Cl_i0) / delta_t, div_i * f_Cl_pump_exp2)
            eq_Clo_t_pump_exp2 = sp.Eq((Cl_o1 - Cl_o0) / delta_t, div_o * f_Cl_pump_exp2)

            J_cond_pump_exp2 = (z_Na * F * f_Na_pump_exp2 +
                                z_K * F * f_K_pump_exp2 +
                                z_Cl * F * f_Cl_pump_exp2
                                )

            eq_Vm_t_pump_exp2 = sp.Eq((V_mem1 - V_mem0) / delta_t, (1 / c_mem) * J_cond_pump_exp2)

            if update_env:
                sol_exp2 = sp.solve((eq_Nai_t_pump_exp2,
                                     eq_Nao_t_pump_exp2,
                                     eq_Ki_t_pump_exp2,
                                     eq_Ko_t_pump_exp2,
                                     eq_Cli_t_pump_exp2,
                                     eq_Clo_t_pump_exp2,
                                     eq_Vm_t_pump_exp2),
                                    (Na_i1, Na_o1, K_i1, K_o1, Cl_i1, Cl_o1, V_mem1)
                                    )
                Nao_pump_funk = sp.lambdify((sym_param_list),
                                            sol_exp2[Na_o1])

                Ko_pump_funk = sp.lambdify((sym_param_list),
                                           sol_exp2[K_o1])

                Clo_pump_funk = sp.lambdify((sym_param_list),
                                            sol_exp2[Cl_o1])

                update_eqn_dict['Nao'] = Nao_pump_funk
                update_eqn_dict['Ko'] = Ko_pump_funk
                update_eqn_dict['Clo'] = Clo_pump_funk


            else:
                sol_exp2 = sp.solve((eq_Nai_t_pump_exp2,
                                     eq_Ki_t_pump_exp2,
                                     eq_Cli_t_pump_exp2,
                                     eq_Vm_t_pump_exp2),
                                    (Na_i1, K_i1, Cl_i1, V_mem1)
                                    )

            Nai_pump_funk = sp.lambdify((sym_param_list),
                                          sol_exp2[Na_i1])

            Ki_pump_funk = sp.lambdify((sym_param_list),
                                          sol_exp2[K_i1])

            Cli_pump_funk = sp.lambdify((sym_param_list),
                                          sol_exp2[Cl_i1])

            Vmem_pump_funk = sp.lambdify((sym_param_list),
                                          sol_exp2[V_mem1])

        update_eqn_dict['Nai'] = Nai_pump_funk
        update_eqn_dict['Ki'] = Ki_pump_funk
        update_eqn_dict['Cli'] = Cli_pump_funk
        update_eqn_dict['Vmem'] = Vmem_pump_funk

        self.update_eqn_dict = update_eqn_dict

    @beartype
    def time_loop(self, p: CalculionParams):
        '''

        '''

        Na_i_time = []
        K_i_time = []
        Cl_i_time = []
        V_mem_time = []
        time = []

        param_change_time = []

        Na_o_time = []
        K_o_time = []
        Cl_o_time = []

        tt = 0.0  # First time value

        # Estimate initial Vmem from parameters assuming system is at electrical steady state (J=0):
        V_mem = vmem_ghk_pump(p.P_Na, p.Na_o, p.Na_i,
                              p.P_K, p.K_o, p.K_i,
                              p.P_Cl, p.Cl_o, p.Cl_i,
                              p.Keqm_NaK, p.omega_NaK,
                              p.ATP, p.ADP, p.P, p.alpha)

        # Initialize the list of variables:
        var_list = [p.Na_i, p.Na_o, p.K_i, p.K_o, p.Cl_i, p.Cl_o, V_mem]

        for ii in range(p.N_iter):

            # Append all values to the time-storage arrays:
            Na_o_time.append(var_list[self.i_Nao])
            K_o_time.append(var_list[self.i_Ko])
            Cl_o_time.append(var_list[self.i_Clo])
            Na_i_time.append(var_list[self.i_Nai])
            K_i_time.append(var_list[self.i_Ki])
            Cl_i_time.append(var_list[self.i_Cli])
            V_mem_time.append(var_list[self.i_Vmem])
            time.append(tt)

            # update the variables:
            var_list, change_list = update_all_var(var_list,
                                                       p,
                                                       self._quasi_static_vmem,
                                                       self._update_env)


            # update the time-step:
            tt += p.delta_t

            # Calculate the change in parameters from last time step:
            self.param_change = np.asarray(change_list)
            param_change_time.append(self.param_change)

        # Assign storage arrays to the object:
        self.Na_i_time = np.asarray(Na_i_time)
        self.Na_o_time = np.asarray(Na_o_time)
        self.K_i_time = np.asarray(K_i_time)
        self.K_o_time = np.asarray(K_o_time)
        self.Cl_i_time = np.asarray(Cl_i_time)
        self.Cl_o_time = np.asarray(Cl_o_time)
        self.V_mem_time = np.asarray(V_mem_time)
        self.time = np.asarray(time)

        self.param_list = var_list

        self.param_change_time = param_change_time


    @beartype
    def time_loop1(self, p: CalculionParams):
        '''

        '''

        if self.update_eqn_dict is None:
            self._get_update_eqns(quasi_static_vmem = self._quasi_static_vmem,
                              update_env = self._update_env)

        Na_i_time = []
        K_i_time = []
        Cl_i_time = []
        V_mem_time = []
        time = []

        param_change_time = []

        Na_o_time = []
        K_o_time = []
        Cl_o_time = []

        tt = 0.0  # First time value

        # Make the reaction (stoichiometry) matrix:
        # Msys = sim.make_reaction_matrix(p.r_cell, r_env, p.d_mem, p.e_r, include_vmem=t_update_vmem)

        # Estimate initial Vmem from parameters assuming system is at electrical steady state (J=0):
        V_mem = vmem_ghk_pump(p.P_Na, p.Na_o, p.Na_i,
                              p.P_K, p.K_o, p.K_i,
                              p.P_Cl, p.Cl_o, p.Cl_i,
                              p.Keqm_NaK, p.omega_NaK,
                              p.ATP, p.ADP, p.P, p.alpha)

        # Initialize the parameter list:
        param_list = [p.Na_i, p.Na_o, p.K_i, p.K_o, p.Cl_i, p.Cl_o, V_mem,
                          p.ATP, p.ADP, p.P, p.Keqm_NaK, p.omega_NaK,
                          p.P_Na, p.P_K, p.P_Cl, p.z_Na, p.z_K, p.z_Cl, p.alpha, p.F,
                          p.r_cell, p.d_ecm, p.c_mem, p.delta_t]

        for ii in range(p.N_iter):

            param_list_o = param_list*1

            # Append all values to the time-storage arrays:
            Na_o_time.append(param_list[self.i_Nao])
            K_o_time.append(param_list[self.i_Ko])
            Cl_o_time.append(param_list[self.i_Clo])
            Na_i_time.append(param_list[self.i_Nai])
            K_i_time.append(param_list[self.i_Ki])
            Cl_i_time.append(param_list[self.i_Cli])
            V_mem_time.append(param_list[self.i_Vmem])
            time.append(tt)

            # update values:
            param_list[self.i_Nai] = self.update_eqn_dict['Nai'](*param_list)
            param_list[self.i_Ki] = self.update_eqn_dict['Ki'](*param_list)
            param_list[self.i_Cli] = self.update_eqn_dict['Cli'](*param_list)
            param_list[self.i_Vmem] = self.update_eqn_dict['Vmem'](*param_list)

            if self._update_env:
                param_list[self.i_Nao] = self.update_eqn_dict['Nao'](*param_list)
                param_list[self.i_Ko] = self.update_eqn_dict['Ko'](*param_list)
                param_list[self.i_Clo] = self.update_eqn_dict['Clo'](*param_list)

            # update the time-step:
            tt += p.delta_t

            # Calculate the change in parameters from last time step:
            change_rate = np.asarray(param_list[0:7]) - np.asarray(param_list_o[0:7])
            self.param_change = change_rate
            param_change_time.append(change_rate)

        # Assign storage arrays to the object:
        self.Na_i_time = np.asarray(Na_i_time)
        self.Na_o_time = np.asarray(Na_o_time)
        self.K_i_time = np.asarray(K_i_time)
        self.K_o_time = np.asarray(K_o_time)
        self.Cl_i_time = np.asarray(Cl_i_time)
        self.Cl_o_time = np.asarray(Cl_o_time)
        self.V_mem_time = np.asarray(V_mem_time)
        self.time = np.asarray(time)

        self.param_list = param_list

        self.param_change_time = param_change_time

    def time_loop2(self, p: CalculionParams):
        '''

        '''

        if self.Msys is None:
            self.Msys = make_reaction_matrix(p, quasi_static_vmem=self._quasi_static_vmem)

        Na_i_time = []
        K_i_time = []
        Cl_i_time = []
        V_mem_time = []
        time = []

        param_change_time = []

        Na_o_time = []
        K_o_time = []
        Cl_o_time = []

        tt = 0.0  # First time value

        # Initial values:
        Na_o = p.Na_o
        Na_i = p.Na_i
        K_o = p.K_o
        K_i = p.K_i
        Cl_o = p.Cl_o
        Cl_i = p.Cl_i


        # Make the reaction (stoichiometry) matrix:
        # Msys = sim.make_reaction_matrix(p.r_cell, r_env, p.d_mem, p.e_r, include_vmem=t_update_vmem)

        # Estimate initial Vmem from parameters assuming system is at electrical steady state (J=0):
        V_mem = vmem_ghk_pump(p.P_Na, p.Na_o, p.Na_i,
                              p.P_K, p.K_o, p.K_i,
                              p.P_Cl, p.Cl_o, p.Cl_i,
                              p.Keqm_NaK, p.omega_NaK,
                              p.ATP, p.ADP, p.P, p.alpha)

        for ii in range(p.N_iter):

            Na_o_time.append(Na_o)
            Na_i_time.append(Na_i)
            K_o_time.append(K_o)
            K_i_time.append(K_i)
            Cl_o_time.append(Cl_o)
            Cl_i_time.append(Cl_i)
            V_mem_time.append(V_mem)
            time.append(tt)

            V_mem_o = V_mem*1

            # Compute the flux vector for the present value of parameters:
            [f_Na, f_K, f_Cl, f_NaKpump] = compute_flux_vector(p.P_Na, Na_o, Na_i,
                                                               p.P_K, K_o, K_i,
                                                               p.P_Cl, Cl_o, Cl_i,
                                                               p.ATP, p.ADP, p.P,
                                                               p.Keqm_NaK, p.omega_NaK, V_mem, p.alpha)

            if self._quasi_static_vmem is True:

                [dNa_i, dNa_o, dK_i, dK_o, dCl_i, dCl_o] = self.Msys.dot([f_Na, f_K, f_Cl, f_NaKpump])

                # Estimate Vmem from parameters assuming system comes to rapid electrical steady state (Jmem=0):
                V_mem = vmem_ghk_pump(p.P_Na, Na_o, Na_i,
                                      p.P_K, K_o, K_i,
                                      p.P_Cl, Cl_o, Cl_i,
                                      p.Keqm_NaK, p.omega_NaK,
                                      p.ATP, p.ADP, p.P, p.alpha)

                dVm = V_mem - V_mem_o


            else:
                [dNa_i, dNa_o, dK_i, dK_o, dCl_i, dCl_o, dVm] = self.Msys.dot([f_Na, f_K, f_Cl, f_NaKpump])

                V_mem += p.delta_t * dVm


            if self._update_env:
                Na_o += p.delta_t * dNa_o
                K_o += p.delta_t * dK_o
                Cl_o += p.delta_t * dCl_o

            Na_i += p.delta_t * dNa_i
            K_i += p.delta_t * dK_i
            Cl_i +=p.delta_t * dCl_i

            tt += p.delta_t

            if self._update_env:
                self.param_change = np.asarray([dNa_i, dNa_o, dK_i, dK_o, dCl_i, dCl_o, dVm])
            else:
                self.param_change = np.asarray([dNa_i, 0.0, dK_i, 0.0, dCl_i, 0.0, dVm])

            param_change_time.append(self.param_change)

        self.Na_o_time = np.asarray(Na_o_time)
        self.Na_i_time = np.asarray(Na_i_time)
        self.K_o_time = np.asarray(K_o_time)
        self.K_i_time = np.asarray(K_i_time)
        self.Cl_o_time = np.asarray(Cl_o_time)
        self.Cl_i_time = np.asarray(Cl_i_time)
        self.V_mem_time = np.asarray(V_mem_time)
        self.time = np.asarray(time)
        self.param_change_time = np.asarray(param_change_time)




