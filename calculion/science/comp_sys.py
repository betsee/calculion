#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright (c) 2022-2023 Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Creates a computational system to calculate steady-state equilibrium concentrations and Vmem for the complete
bioelectric system from user-specified properties, as well as time-dependent concentrations and Vmem. Calculations are
made from analytical equations.
'''
from beartype import beartype, BeartypeConf
from beartype.typing import Optional
import sympy as sp
import numpy as np
from scipy.optimize import minimize
import pandas as pd
from pandas import DataFrame
from calculion.science.params import CalculionParams

@beartype
class CompSys(object):
    '''
    A computational system to calculating steady-state concentrations and Vmem for the complete
    bioelectric system from user-specified properties, as well as time-dependent concentrations and Vmem.
    Calculations are made from analytical equations. Analytical equations can be viewed as LaTex expressions.

    '''

    def __int__(self, p: Optional[CalculionParams]=None):
        '''
        Initialize the computational system with default parameters.
        '''

        if p is None:
            self._p = CalculionParams()
        else:
            self._p = p

    def _analytical_sys(self):
        '''
        Set up the analytical (symbolic) system.
        '''

        # Define symbols:

        # Voltage, measurements, constants:
        V_mem, z_Na, z_K, z_Cl, = sp.symbols('V_mem, z_Na, z_K, z_Cl', real=True)

        alpha, d_ecm, r_cell, c_mem, F, R, T = sp.symbols('alpha, d_ecm, r_cell, c_mem, F, R, T',
                                                          real=True, positive=True)
        r_env = r_cell + d_ecm

        # Divergence operators (assuming constant flux across the membrane):
        div_i = (2 / r_cell)
        div_o = -((2 * r_cell) / (r_env ** 2 - r_cell ** 2))

        # Concentrations inside and out of the cell:
        Na_o, Na_i, K_o, K_i, Cl_o, Cl_i, ATP, ADP, P = sp.symbols('Na_o, Na_i, K_o, K_i, Cl_o, Cl_i, ATP, ADP, P',
                                                                   real=True, positive=True)
        # Membrane permeabilities
        P_Na, P_K, P_Cl = sp.symbols('P_Na, P_K, P_Cl', real=True, positive=True)

        # ATP concentrations, free energies of hydrolysis, rates of Na-K-ATPase and metabolism:
        (K_eqm_NaK,
         omega_pump,
         omega_met) = sp.symbols('K_eqm_NaK, Omega_pump, Omega_met', real=True)

        # Individual ion fluxes: # NOTE: into the cell direction is positive flux
        # Cation fluxes:
        # Note Cl (anion) reaction is written backwards to ge the same exp(alpha*V_mem) term in the final expression
        f_Na = P_Na * (Na_o - Na_i * sp.exp(alpha * V_mem))  # Na_o <--> Na_i  delGo_Na = 0
        f_K = P_K * (K_o - K_i * sp.exp(alpha * V_mem))  # K_o <--> K_i  delGo_K = 0
        f_Cl = -P_Cl * (Cl_i - Cl_o * sp.exp(alpha * V_mem))  # Cl_i <--> Cl_o  delGo_K = 0

        # Nernst potentials (reversal potentials) of the individual ions:
        Na_rev_Eq = sp.solve(f_Na, V_mem)[0]
        K_rev_Eq = sp.solve(f_K, V_mem)[0]
        Cl_rev_Eq = sp.solve(f_Cl, V_mem)[0]

        # Pump and transporter fluxes:

        # Na,K-ATPase pump flux:
        # pump reaction is written backwards to get the same exp(alpha*V_mem) term in the final expression
        # 3 Na_o + 2 K_i + ADP + P <--> 3 Na_i + 2 K_o + ATP,  -delGoATP = +28e3 to +34e3 J/mol
        f_NaKpump = (omega_pump * (K_eqm_NaK * ADP * P * (K_i ** 2) * (Na_o ** 3) -
                                   ATP * (K_o ** 2) * (Na_i ** 3) * sp.exp(alpha * V_mem))
                     )
        # Equation for the Na-K-Cl transporter:
        # 1 Na_o + 1 K_o + 2 Cl_o <--> 1 Na_i + 1 K_i + 2 Cl_i
        omega_NaKCl = sp.symbols('Omega_NaKCl', real=True)
        f_NaKCl = omega_NaKCl * (Na_o * K_o * Cl_o ** 2 - Na_i * K_i * Cl_i ** 2)

        # Equation for the K-Cl symporter:
        # 1 K_o + Cl_o <--> 1 K_i + Cl_i
        omega_KCl = sp.symbols('Omega_KCl', real=True)
        f_KCl = omega_KCl * (K_o * Cl_o - K_i * Cl_i)

        # Full system: Na, K-ATPase pumps, Na-K-2Cl cotransporter, K-Cl symporter; Solving for ions in cell only:

        flux_vect = sp.Matrix([f_Na, f_K, f_Cl, f_NaKpump, f_NaKCl, f_KCl])

        M_sys = sp.Matrix([[div_i, 0, 0, 3*div_i, div_i, 0],
                           [0, div_i, 0, -2*div_i, div_i, div_i],
                           [0, 0, div_i, 0, 2*div_i, div_i],
                           [F/c_mem, F/c_mem, -F/c_mem, F/c_mem, 0, 0]
                                   ])

        consts_vect = [Na_o, K_o, Cl_o, F, P_Na, P_K, P_Cl, r_cell, c_mem, alpha,
                       omega_pump, K_eqm_NaK, ATP, ADP, P, omega_NaKCl, omega_KCl]


        # Common to fixed extracellular, solving for all ions in and Vmem:

        zer_vect = sp.Matrix([0,0,0,0])

        params_vect = [Na_i, K_i, Cl_i]

        par_names = ['Na_i', 'K_i', 'Cl_i']

        self.bounds_vect = [(0.0, 1000.0), (0.0, 1000.0), (0.0, 1000.0)]

        # Evolution equations for the system:
        evo_Eq_o = (M_sys * flux_vect)

        dt_Vmem_Eq = evo_Eq_o[-1]

        # Vmem calculated using a quasi-steady-state approximation (solving J=0 for V_mem)
        # this is equivalent to the GHK equation for case of no ion pumps/transporters:
        Vmem_Eq = sp.solve(dt_Vmem_Eq, V_mem)[0]

        # Substitute in the quasi-static solution of Vmem to the system and truncate evolution equations
        # to only deal with changes to concentrations:
        evo_Eq = sp.Matrix(evo_Eq_o.subs(V_mem, Vmem_Eq)[0:-1])

        opti_Eq = ((evo_Eq.T * evo_Eq)[0])  # Objective function to minimize (symbolic version)
        all_params_vect = [params_vect, consts_vect]  # all model parameters (symbolic)
        jac_Eq = [sp.diff(opti_Eq, pi, 1) for pi in params_vect]  # symbolic system Jacobian
        hess_Eq = [[opti_Eq.diff(pj).diff(pi) for pj in params_vect] for pi in params_vect]  # symbolic system Hessian

        # Lambdify the main expressions:
        self.evo_funk = sp.lambdify(all_params_vect, evo_Eq)
        self.opti_funk = sp.lambdify(all_params_vect, opti_Eq)
        self.jac_funk = sp.lambdify(all_params_vect, jac_Eq)
        self.hess_funk = sp.lambdify(all_params_vect, hess_Eq)
        self.vmem_funk = sp.lambdify(all_params_vect, Vmem_Eq)


    def calc_steady_state(self, p):

        p_params_vect = np.asarray([p.Na_i, p.K_i, p.Cl_i])

        p_consts_vect = [p.Na_o, p.K_o, p.Cl_o, p.F, p.P_Na, p.P_K, p.P_Cl, p.r_cell, p.c_mem, p.alpha,
                        p.omega_NaK, p.Keqm_NaK, p.ATP, p.ADP, p.P, p.omega_NaKCl, p.omega_KCl]

        vmem_o = self.vmem_funk(p_params_vect, p_consts_vect)  # initialize Vmem

        sol0 = minimize(self.opti_funk,
                        p_params_vect,
                        args=(p_consts_vect),
                        method='trust-constr',
                        jac=self.jac_funk,
                        hess=self.hess_funk,
                        tol=1.0e-30,
                        bounds=self.bounds_vect
                        )



    def calc_timestepped(self, p, N_iter, del_t, ti: float = 0.0):

        parami_vect = np.asarray([p.Na_i, p.K_i, p.Cl_i])

        consti_vect = [p.Na_o, p.K_o, p.Cl_o, p.F, p.P_Na, p.P_K, p.P_Cl, p.r_cell, p.c_mem, p.alpha,
                        p.omega_NaK, p.Keqm_NaK, p.ATP, p.ADP, p.P, p.omega_NaKCl, p.omega_KCl]

        param_time = []
        opti_funk_time = []
        time_vect = []


        paramj_vect = np.copy(parami_vect)

        for ni in range(N_iter):
            ti += del_t

            paramj_dt = self.evo_funk(paramj_vect, consti_vect).flatten()

            paramj_vect += del_t * paramj_dt

            opti_funki = self.opti_funk(paramj_vect, consti_vect)

            param_time.append(paramj_vect * 1)
            opti_funk_time.append(opti_funki * 1)
            time_vect.append(ti)

        param_time = np.asarray(param_time)
        opti_funk_time = np.asarray(opti_funk_time)

        return param_time, opti_funk_time, time_vect