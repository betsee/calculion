#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright (c) 2022-2023 Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
This module creates systems for working with ions and ionic reactions in IonSpire.

'''
import copy
from beartype import beartype
from beartype.typing import Optional, Union
import sympy as sp
from sympy.core.symbol import Symbol
import numpy as np
from numpy import ndarray
from scipy.optimize import minimize, basinhopping
from calculion.science.channel_base import DynamicChannelABC
from calculion.science.chem_base import ReactionABC
from calculion.science.reaction_system import ReactionSystem
from calculion.science.model_params import ModelParams


@beartype
class SteadyStateOpti(object):
    '''
    Use the Gauss-Newton method to find the steady-state values of a bioelectric
    system's concentrations and transmembrane voltage.

    Public Attributes
    ------------------
    args_list : list[Symbol]
        List of constant parameters (as Sympy Symbols) that don't change in the optimization.

    args_vals : list[Union[float, int]]
        List of values for constant parameters.

    divflux_f : function
        Function that computes the divergence of the flux and transmembrane current, which is
        also the residuals vector of the optimization. This has been created from the analytic
        Sympy expression divflux_s.

    divflux_params_list : list
        Sympy Symbols of the ordered parameters taken by the divflux_f function.

    jac_f : function
        Function that computes the Jacobian, which is used in each step of the optimization.
        This has been created from the analytic Sympy expression jac_s.

    jac_params_o : list[Union[float, int]]
        List of parameter values for the Jacobian.

    keep_positive : list[bool]
        For each parameter listed in ss_params_list, keep_positive lists a boolean,
        which is True if the respective parameter should remain above zero.

    param_o : list
        List of parameter values that are optimized by the solver.

    ss_params_list : list[Symbol]
        List of Sympy Symbols that will be the free parameters of the optimization.
        Typicaly this is a list of symbols from the list bes.chem_vect_reduced,
        augmented by V_mem_s.

    Private Attributes
    ------------------

    _bes : ReactionSystem
        Instance of a ReactionSystem with parameters to be optimized.

    _ss_jac_s : MutableDenseMatrix
        Analytic representation of the jacobian matrix for the optimization system.

    '''

    def __init__(self,
                 bes: ReactionSystem,
                 ss_params_list: list[Symbol],
                 keep_positive: Optional[list[bool]]=None,
                 norm_params: float=1.0
                 ):
        '''
        Initialize the SteadyStateOpti object.

        Parameters
        ------------
        bes : ReactionSystem
            Instance of reaction system.

        ss_params_list : list[Symbol]
            List of Sympy Symbols that will be the free parameters of the optimization.
            Typicaly this is a list of symbols from the list bes.chem_vect_reduced,
            augmented by V_mem_s.

        keep_positive: Optional[list[bool]]
            For each parameter listed in ss_params_list, keep_positive lists a boolean,
            which is True if the respective parameter should remain above zero.

        '''
        # obtain a deep copy of the bioelectric system to avoid modifying values
        self._bes = copy.deepcopy(bes)

        # Parameters to be optimized in steady-state search:
        # Initial value of parameters:
        self.param_o = (1/norm_params)*np.asarray(bes._get_param_vals(ss_params_list))

        if keep_positive is None:
            self.keep_positive = [False for i in ss_params_list]
        else:
            self.keep_positive = keep_positive

        self.divflux_params_list = bes.divflux_params_list
        self.divflux_f = bes.divflux_f

        # Constant parameters that don't change
        self.args_list = np.setdiff1d(np.asarray(self.divflux_params_list), np.asarray(ss_params_list),
                                 assume_unique=True).tolist()

        # Value of constants:
        self.args_vals = bes._get_param_vals(self.args_list)

        # We need to calculate the Jacobian:
        ss_jac_s = bes.divflux_s.jacobian(sp.Matrix(ss_params_list))
        self._ss_jac_s = ss_jac_s

        # Lambdify (get numerical function) for the V_flux vector:
        jac_params_list = list(ss_jac_s.free_symbols)

        self.jac_f = sp.lambdify((ss_params_list, self.args_list), ss_jac_s)

        # Initialize the float values to use in function call to div_flux_f:
        self.jac_params_o = bes._get_param_vals(jac_params_list)

        self.ss_params_list = ss_params_list

        # Create a second version of divfluxf that can be used with scipy's optimization
        # functions:
        # Symbolic optimization function, 1 output variable:
        self.divflux_opti_s = (bes.divflux_s*bes.divflux_s.T)
        # Numerical optimization function, 1 output variable:
        self.divflux_opti_f = sp.lambdify((ss_params_list, self.args_list),
                                          self.divflux_opti_s[0])
        # Symbolic Jacobian for 1 output variable optimization function:
        self.jac_opti_s = self.divflux_opti_s.jacobian(ss_params_list)
        # Numerical Jacobian for 1 output variable optimization function:
        self.jac_opti_f = sp.lambdify((ss_params_list, self.args_list),
                                          sp.Array(self.jac_opti_s)[0, :])

        # Get the symbolic Hessian for the 1 output variable optimization function:
        self.hess_opti_s = self.jac_opti_s.jacobian(ss_params_list)
        # And the numerical hessian:
        self.hess_opti_f = sp.lambdify((ss_params_list, self.args_list), self.hess_opti_s)

    def print_results(self, ss_params: Union[list, ndarray]):
        '''
        Prints the results of the optimization.

        Parameters
        ----------
        ss_params: Union[list, ndarray]
            Values of the optimized parameters, ordered according to self.ss_params_list.

        '''
        for i, (pnme, pval) in enumerate(zip(self.ss_params_list, ss_params)):
            if i == len(ss_params) - 1:
                print(pnme, pval * 1e3)
            else:
                print(pnme, pval)

    def set_results(self, bes, ss_params):
        '''
        Write the results in ss_params to the ReactionSystem object.

        Parameters
        ----------
        bes: ReactionSystem
            Instance of ReactionSystem to which parameter values are updated.

        ss_params: Union[list, ndarray]
            Values of the optimized parameters, ordered according to self.ss_params_list.
        '''

        # Update all parameters in the model:
        for i, (pnme, pval) in enumerate(zip(self.ss_params_list, ss_params)):
            setattr(bes, pnme.name, pval)

        # update the full chem vector:
        bes.chem_vals = bes._get_chem_vals(bes.chem_vect)

        # And the reduced chem vector:
        bes.chem_vals_reduced = bes._get_param_vals(bes.chem_vect_reduced)

class IterSim(object):
    '''
    Run an iterative simulation in time using the explicit Euler technique.

    Public Attributes
    -----------------
    time : ndarray
        Stores the time value at each iteration of the simulation.

    vm_time : ndarray
        Stores the Vmem value at each iteration of the simulation.

    jc_time : ndarray
        Stores the conduction current value at each iteration of the simulation.

    chem_time : ndarray
        Stores the concentration values at each iteration of the simulation.

    pmem_time : ndarray
        Stores the membrane permeabilities/rate constants at each iteration of the simulation.

    chem_sim : list
        Value of concentrations at active iteration step.

    chem_sim_name : list
        Name of concentrations.

    t : float
        Value of time at active iteration step.

    Private Attributes
    -------------------
    _bes : ReactionSystem
        An instance of ReactionSystem, defines the system being simulated by explicit Euler.

    _channels_list : NoneType or list
        List of all channels being simulated in the simulation as DynamicChannelABC objects.

    '''

    # @beartype
    def __init__(self,
                 bes: ReactionSystem,
                 channels_list: Optional[list[DynamicChannelABC]]=None,
                 ):
        '''
        bes : ReactionSystem
            An instance of ReactionSystem defining all reactions of the system to be
        optimized.

        channels_list : list[DynamicChannelABC]
            A list of all the channels that are to be included in the dynamic simulation.

        '''

        # Make a copy of the bioelectrical system:
        self._bes = copy.deepcopy(bes)
        self._channels_list = channels_list

        self.t = 0.0

        self.chem_sim = []
        self.chem_sim_name = []

        # Get the string names of the simulated chemicals:
        for chm in self._bes.chem_vect_reduced:
            self.chem_sim.append(chm)
            self.chem_sim_name.append(chm.name)


    def run_iter_sim(self,
                     dt: float,
                     N_iter: int,
                     use_quasi_static_approx: bool=False,
                     use_hodgkin_huxley: bool = False,
                     clamp_vmem_at: Optional[float]=None,
                     sweep_vmem_vals: Optional[tuple[float, float]]=None):
        '''
        Run an iterative simulation using the system-calculated
        equations to update chemical concentrations and Vmem.

        Parameters
        ---------------
        dt : float
            Time step to use in the iterative simulation.

        N_iter : int
            Total number of iterations to calculate.

        use_quasi_static_approx : bool=False
            Use the Jc=0 "quasi-static" approximation for Vmem?

        use_hodgkin_huxley : bool = False
            Use the Hodgkin-Huxley equation to compute changes to Vmem?

        clamp_vmem_at : Optional[float]=None
            Clamp Vmem at this voltage? If None, Vmem changes freely.

        sweep_vmem_vals : Optional[tuple[float, float]]=None
            Sweep Vmem from the first voltage through to the second voltage?
            If None, Vmem changes freely.


        '''

        # If using the Hodgkin-Huxley model, don't update concentrations and calculate
        # conduction current using a different equation:
        self.use_hodgkin_huxley = use_hodgkin_huxley

        # Calculate the divergence of the flux for this system:
        vm_time = []
        jc_time = [] # Conduction current density as a function of time
        chem_time = []
        time = []

        self.resting_Pmem_dict = {}
        self.channel_pmem_names = set()

        if self._channels_list is not None:
            for chan in self._channels_list:
                for ion_pmem_name in chan.ion_perm_list:
                    self.channel_pmem_names.add(ion_pmem_name)

                # collect resting permittivity values
                self.resting_Pmem_dict = self._bes.get_pmem_vals(chan.ion_perm_list,
                                                                 self.resting_Pmem_dict)

                # Initialize the channel to the initial Vmem of the BES:
                chan.init_channel(self._bes.V_mem)

        pmem_time = []

        if clamp_vmem_at is not None:
            self._bes.V_mem = clamp_vmem_at # Initialize Vmem to the clamp voltage

        if sweep_vmem_vals is not None:
            v0 = sweep_vmem_vals[0]
            v1 = sweep_vmem_vals[1]
            dvm = (v1 - v0)/N_iter
            self._bes.V_mem = v0

        else:
            dvm = 0.0

        for ti in range(N_iter):
            # print(ti)
            self.t = ti*dt # update the time parameter

            # If ion channels are supplied compute their influence on BioElectrical System membrane permeabilities:
            if self._channels_list is not None:
                chan_pmem_mod_dict = {}
                for k, v in self.resting_Pmem_dict.items():
                    # FIXME: if vmem is a float, this works but not for an array
                    chan_pmem_mod_dict[k] = 0.0 # set entries for necessary ions to zero
                for chan in self._channels_list:
                    # get default values of membrane permeabilities acted on by this channel

                    perm_mod = chan.run_channel(self._bes.V_mem,
                                                dt,
                                                self.t) # channel updates the permittivity

                    for ion_pmem_name in chan.ion_perm_list:
                        chan_pmem_mod_dict[ion_pmem_name] += perm_mod

                for ion_pmem_name in self.channel_pmem_names: # for this ion's pmem...
                    pi_rest = self.resting_Pmem_dict[ion_pmem_name] # get resting pmem
                    pi_chan = chan_pmem_mod_dict[ion_pmem_name] # get pmem mod by channel
                    pi_total = pi_rest + pi_chan # add them together
                    setattr(self._bes, ion_pmem_name, pi_total) # set it to the BES

            if use_hodgkin_huxley is False:

                # Calculate the change array vector.
                change_array = self._bes.divflux_f(*self._bes.divflux_params)[0]

                self._bes.chem_vals_reduced += change_array[0:-1]*dt

                # Update chemical values on the main object:
                self._bes._set_chem_vals()

                if use_quasi_static_approx is True:

                    if clamp_vmem_at is None and sweep_vmem_vals is None:
                        self._bes.V_mem = self._bes.solve_ss_vmem()

                    elif sweep_vmem_vals is not None:
                        self._bes.V_mem += dvm

                    self.jc = 0.0

                else:
                    self.jc = self._bes.c_mem*change_array[-1]

                    if clamp_vmem_at is None and sweep_vmem_vals is None:
                        # Calculate V_mem using the steady-state approximation (Jc=0) derived formula:
                        self._bes.V_mem += change_array[-1]*dt

                    elif sweep_vmem_vals is not None:
                        self._bes.V_mem += dvm

                # Update the input parameters to divflux_f and vss_f functions:
                self._bes.divflux_params = self._bes._get_param_vals(self._bes.divflux_params_list)

            # If we're using the Hodgkin-Huxley approximation for conduction current, skip updating
            # ion concentrations entirely and calculate conduction current using a linearized equation:
            else:
                self.jc = -self._bes.get_jc_hh()

                if clamp_vmem_at is None and sweep_vmem_vals is None:
                    self._bes.V_mem += (1/self._bes.c_mem)*self.jc*dt

                elif sweep_vmem_vals is not None:
                    self._bes.V_mem += dvm


            time.append(self.t)
            vm_time.append(self._bes.V_mem * 1)
            jc_time.append(self.jc*1)
            chem_time.append(self._bes.chem_vals_reduced * 1)
            pmem_time.append(np.asarray([self._bes.P_Na, self._bes.P_K, self._bes.P_Cl]))
            # pmem_time.append(perm_vals)

        self.time = np.asarray(time)
        self.vm_time = np.asarray(vm_time)
        self.jc_time = np.asarray(jc_time)
        self.chem_time = np.asarray(chem_time)
        self.pmem_time = np.asarray(pmem_time)

        return self.time, self.vm_time, self.chem_time