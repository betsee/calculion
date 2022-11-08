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

    def __int__(self, tol=1.0e-15):
        '''
        Initialize the computational system with default parameters.
        '''
        self._tol = tol

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

        # zer_vect = sp.Matrix([0,0,0,0])

        params_vect = [Na_i, K_i, Cl_i]

        self.param_names = ['Na_i', 'K_i', 'Cl_i']

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

        # Numerical versions of reversal potentials:
        self.Na_rev_funk = sp.lambdify(all_params_vect, Na_rev_Eq)
        self.K_rev_funk = sp.lambdify(all_params_vect, K_rev_Eq)
        self.Cl_rev_funk = sp.lambdify(all_params_vect, Cl_rev_Eq)

    def collect_params(self, p: CalculionParams):
        '''
        Collect and return a parameters vector and a constants vector from a Calculion Params object.
        '''

        params_vect = np.asarray([p.Na_i, p.K_i, p.Cl_i])

        consts_vect = [p.Na_o, p.K_o, p.Cl_o, p.F, p.P_Na, p.P_K, p.P_Cl, p.r_cell, p.c_mem, p.alpha,
                        p.omega_NaK, p.Keqm_NaK, p.ATP, p.ADP, p.P, p.omega_NaKCl, p.omega_KCl]

        return params_vect, consts_vect

    def calc_steady_state(self, p: CalculionParams):
        '''
        Use Scipy's Trust-Const minimization routine to find the solution of the system at steady-state.

        Parameters
        -----------
        '''

        # Initial electrical properties of the system:
        # self.elec_props_o = self.calc_electrical_properties(p_params_vect, p_consts_vect)
        # Create the parameters and constants vector from the Calculion Params instance:
        p_params_vect, p_consts_vect = self.collect_params(p)

        self.sol0 = minimize(self.opti_funk,
                        p_params_vect,
                        args=(p_consts_vect),
                        method='trust-constr',
                        jac=self.jac_funk,
                        hess=self.hess_funk,
                        tol=self._tol,
                        bounds=self.bounds_vect
                        )

        # Final electrical properties of the system:
        # self.elec_props = self.calc_electrical_properties(self.sol0.x, p_consts_vect)

        return self.sol0.x

    def calc_timestepped(self,
                         p: CalculionParams,
                         N_iter: int=1000,
                         del_t: float=1.0,
                         ti: float = 0.0,
                         tol: float=1.0e-9):
        '''
        Use explicit Euler method to solve the time-dependent system.

        Parameters
        -----------
        p: CalculionParams
            An instance of Calculion Parameters object.
        N_iter: int
            The maximum number of itterations to run the system.
        del_t: float
            The time-step over which to integrate each step of the iterative solution.
        ti: float
            The start-time of the simulation.
        tol: float
            The stopping tolerance -- when changes between each step are lower than tol the simulation is
            assumed to be at steady-state and is stopped.
        '''
        print_message = ''

        # Create the parameters and constants vector from the Calculion Params instance:
        parami_vect, consti_vect = self.collect_params(p)

        param_time = []
        opti_funk_time = []
        time_vect = []

        paramj_vect = np.copy(parami_vect)

        for ii, ni in enumerate(range(N_iter)):
            ti += del_t # Advance the time value by the time-step

            paramj_dt = self.evo_funk(paramj_vect, consti_vect).flatten() # Calculate the change to each parameter

            paramj_vect += del_t * paramj_dt # Integrate each parameter

            opti_funki = self.opti_funk(paramj_vect, consti_vect) # Compute the value of the optimization function

            param_time.append(paramj_vect * 1)
            opti_funk_time.append(opti_funki * 1)
            time_vect.append(ti)

            if opti_funki < tol:
                print_message = f"System is at steady state after {ii} itterations."
                break

        param_time = np.asarray(param_time)
        opti_funk_time = np.asarray(opti_funk_time)

        return param_time, opti_funk_time, time_vect, print_message

    def return_elec_props_dict(self, parami_v, consti_v):
        '''
        Returns the electrical properties of the system as a Pandas Dataframe
         given a parameters and constants vector.

        '''
        Vmem = self.vmem_funk(parami_v, consti_v)  # Vmem
        Na_rev = self.Na_rev_funk(parami_v, consti_v) # Na reversal
        K_rev = self.K_rev_funk(parami_v, consti_v) # K reversal
        Cl_rev = self.Cl_rev_funk(parami_v, consti_v) # Cl reversal
        Na_ed = Vmem - Na_rev # Electrochemical driving force Na
        K_ed = Vmem - K_rev # Electrochemical driving force K
        Cl_ed = Vmem - Cl_rev # Electrochemical driving force Cl

        unit_conv = 1.0e3 # unit conversion to mv

        # Electrical properties dictionary
        elec_props = {'Vₘ': np.round(Vmem*unit_conv, 1),
                      'Vᵣ Na⁺': np.round(Na_rev*unit_conv, 1),
                      'Vᵣ K⁺': np.round(K_rev*unit_conv, 1),
                      'Vᵣ Cl⁻': np.round(Cl_rev*unit_conv, 1),
                      'Vₑ Na⁺': np.round(Na_ed*unit_conv, 1),
                      'Vₑ K⁺': np.round(K_ed*unit_conv, 1),
                      'Vₑ Cl⁻': np.round(Cl_ed*unit_conv, 1)
        }

        elec_dataframe = pd.DataFrame.from_dict(elec_props,
                                                orient='index', columns=['Voltage [mV]'])

        return elec_dataframe

    def return_chem_props_dict(self, parami_v, consti_v, round_dec: int=2):
        '''
        Returns a Pandas dataframe of the ion concentrations inside and outside of the cell
        given a parameters and constants vector.
        '''

        steady_state_ions_dict = {
            'Na⁺ᵢₙ': round(parami_v[0], round_dec),
            'K⁺ᵢₙ': round(parami_v[1], round_dec),
            'Cl⁻ᵢₙ': round(parami_v[2], round_dec),
            'Na⁺ₒᵤₜ': round(consti_v[0], round_dec),
            'K⁺ₒᵤₜ': round(consti_v[1], round_dec),
            'Cl⁻ₒᵤₜ': round(consti_v[2], round_dec)
                 }

        ions_dataframe = pd.DataFrame.from_dict(steady_state_ions_dict,
                                                orient='index', columns=['Concentration [mM]'])

        return ions_dataframe

    def calc_elec_param_set(self, params_vect_set, consti_v, time):
        '''
        Calculate electrical properties from a set of parameters taken over different times or conditions.
        '''
        Vmem_vect = []
        Na_rev_vect = []
        K_rev_vect = []
        Cl_rev_vect = []
        Na_ed_vect = []
        K_ed_vect = []
        Cl_ed_vect = []

        for parami_v in params_vect_set:

            Vmem = self.vmem_funk(parami_v, consti_v)  # Vmem
            Na_rev = self.Na_rev_funk(parami_v, consti_v) # Na reversal
            K_rev = self.K_rev_funk(parami_v, consti_v) # K reversal
            Cl_rev = self.Cl_rev_funk(parami_v, consti_v) # Cl reversal
            Na_ed = Vmem - Na_rev # Electrochemical driving force Na
            K_ed = Vmem - K_rev # Electrochemical driving force K
            Cl_ed = Vmem - Cl_rev # Electrochemical driving force Cl

            Vmem_vect.append(Vmem)
            Na_rev_vect.append(Na_rev)
            K_rev_vect.append(K_rev)
            Cl_rev_vect.append(Cl_rev)
            Na_ed_vect.append(Na_ed)
            K_ed_vect.append(K_ed)
            Cl_ed_vect.append(Cl_ed)

        Vmem_vect = np.asarray(Vmem_vect)
        Na_rev_vect = np.asarray(Na_rev_vect)
        K_rev_vect = np.asarray(K_rev_vect)
        Cl_rev_vect = np.asarray(Cl_rev_vect)
        Na_ed_vect = np.asarray(Na_ed_vect)
        K_ed_vect = np.asarray(K_ed_vect)
        Cl_ed_vect = np.asarray(Cl_ed_vect)

        # Stack the results of time-dependent properties:
        volt_dat = np.column_stack((time / 3600,
                                 1e3 * Vmem_vect,
                                 1e3 * Na_rev_vect,
                                 1e3 * K_rev_vect,
                                 1e3 * Cl_rev_vect,
                                 1e3 * Na_ed_vect,
                                 1e3 * K_ed_vect,
                                 1e3 * Cl_ed_vect,
                                 ))

        volt_timedat = DataFrame(volt_dat, columns=['Time (hours)',
                                                    'Vₘ',
                                                    'Vᵣ Na⁺',
                                                    'Vᵣ K⁺',
                                                    'Vᵣ Cl⁻',
                                                    'Vₑ Na⁺',
                                                    'Vₑ K⁺',
                                                    'Vₑ Cl⁻'
                                                    ])

        return volt_timedat


    def calc_chem_param_set(self, params_vect_set, consti_v, time):
        '''
        Compile ion concentrations in cells from a set of parameters taken over different times or conditions.
        '''
        Na_vect = []
        K_vect = []
        Cl_vect = []

        for parami in params_vect_set:
            Na_vect.append(parami[0])
            K_vect.append(parami[1])
            Cl_vect.append(parami[2])

        Na_vect = np.asarray(Na_vect)
        K_vect = np.asarray(K_vect)
        Cl_vect = np.asarray(Cl_vect)

        # Stack the results of time-dependent properties:
        chem_dat = np.column_stack((time / 3600, Na_vect, K_vect, Cl_vect))

        chem_timedat = DataFrame(chem_dat, columns=['Time (hours)',
                                                    'Na⁺ᵢₙ (mM)',
                                                    'K⁺ᵢₙ (mM)',
                                                    'Cl⁻ᵢₙ (mM)'
                                                    ])

        return chem_timedat


