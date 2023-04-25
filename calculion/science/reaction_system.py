#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright (c) 2022-2023 Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
This module creates systems for working with ions and ionic reactions in IonSpire.

'''

from beartype import beartype
from beartype.typing import Optional, Union
import sympy as sp
from sympy.core.symbol import Symbol
import numpy as np
from numpy import ndarray
from scipy.optimize import minimize
import pandas as pd
from calculion.science.model_params import ModelParams
from calculion.science.chem_base import (Chemical, TransportReaction, ReactionABC)
from calculion.science.string_names import StringNames
import pydot
from pydot import Dot

class ReactionSystem(object):
    '''
    Create a bioelectric system characterized by a system of reactions and ions
    that are traversing across a thin membrane between an inside 'i' and outside 'o'
    compartment, with a transmembrane voltage, V_mem acting across the membrane.

    Public Attributes
    -------------------
    F : int
        Planck's constant (F=96485 C/mol)

    F_s : Symbol
        Symbolic Planck's constant.

    Keqm_s : int
        Symbolic reaction equilibrium constant.

    R : float
        Ideal gas constant (R=8.314 J/(K mol))

    R_s : Symbol
        Symbolic ideal gas constant.

    T : int
        Reaction temperature (degrees K).

    T_s : Symbol
        Symbolic reaction temperature.

    V_s : Symbol
        Symbolic voltage.

    Vi_s : Symbol
        Symbolic voltage in interior compartment.

    Vm_s : Symbol
        Symbolic transmembrane voltage Vi_s - Vo_s.

    Vo_s : Symbol
        Symbolic voltage in exterior compartment.

    c_mem : float
        Membrane patch capacitance (F/m^2).

    c_norm : float
        Normalizing constant for the reaction, where concentrations are converted into mol/L
        from mol/m3 for thermodynamic purposes.

    c_norm_s : Symbol
        Symbolic concentration normalization constant for the system.

    chem_vals : list
        List of values of all chemicals.

    chem_vals_reduced : list
        List of values of all chemicals that have dynamic concentrations.

    chem_vect : list
        List of all chemicals in the reaction system.

    chem_vect_reduced : list
        List of all chemicals in the reaction system that have dynamic concentrations.

    cm_s : Symbol
        Symbolic representation of membrane capacitance.

    d_ecm : float
        Thickness of the extracellular space [m].

    d_mem : float
        Thickness of the membrane [m].

    div_i : float
        Divergence operator for the interior compartment.

    div_i_s : Symbol
        Symbolic divergence operator for the interior compartment.

    div_o : float
        Divergence operator for the exterior compartment.

    div_o_s : Symbol
        Symbolic divergence operator for the exterior compartment.

    divflux_f : function
        Numeric function computing the divergence of the flux for all reactions as a
        vector of fluxes.

    divflux_params : list
        Numerical parameter values of the divflux_f function.

    divflux_params_list : list
        String attibute names for the numerical parameter values of the divflux_f function.

    divflux_s : MutableDenseMatrix
        Symbolic expression representing the divergence of the flux for all reactions as a
        vector of fluxes. Used to create divflux_f.

    ghk_f : NoneType or function
        Function to compute Vmem using the GHK Voltage equation.

    ghk_params : NoneType or list
        Parameter values for the GHK Voltage function ghk_f.

    ghk_params_list : NoneType or list
        String names of parameter attributes for the GHK Voltage function ghk_f.

    ghk_s : NoneType or Add
        Symbolic expression to compute Vmem using the GHK Voltage equation.

    ion_currents_dict : dict
        Dictionary of ion currents from all reactions summed
        indexed to each transmembrane permeable ion.

    jc_arg_vals : list
        List of numerical constants values (static variables) used in jc_f.

    jc_args_list : list
        List of string names of constant attibutes (static variables) used in jc_f.

    jc_f : function
        Numerical function for computing the total current across the membrane [A].

    jc_params_list : list
        List of parameters (dynamic variables) used in jc_f.

    r_cell : float
        Radius of the cell (inner compartment) [m].

    r_o : float
        Radius of the outer compartment [m].

    react_vect : list
        List of all reactions of the system as ReactionABC objects.

    total_current_s : Add
        Analytic expression representing the total current across the membrane.


    Private Attributes
    ---------------------
    _ignore_chem_tag : list
        List of chemicals that will have concentration updates ignored (remain fixed concentration).

    _ignore_region_tag : str
        List of regions where all chemicals in the region will have fixed concentration.

    _quasi_static_approx : bool
        Use the quasi-static approximation for Vmem?

    _transmem_chem_names : list
        Names of chemicals that traverse the membrane.

    '''

    # @beartype
    def __init__(self,
                 p: ModelParams,
                 chem_vect: list[Chemical],
                 transmem_chem_names: list[str],
                 reaction_sys_vect: list[ReactionABC],
                 rate_base_names: list[tuple[str, Union[float, ndarray]]],
                 quasi_static_approx: bool=False,
                 ignore_region_updates: Optional[str] = None,
                 ignore_chem_updates: Optional[list[Chemical]] = None,
                 ediff_reactions_vector: Optional[list[TransportReaction]] = None
                 ):
        '''
        Initialize the bioelectric system.

        Parameters
        -----------

         p : ModelParams
            Instance of bioelectricity parameters (BioeParams) object.

         chem_vect : list[Chemical]
            List of all chemical entities involved in all reactions constituting the reaction system.

         transmem_chem_names : list[str]
            List of all chemicals that traverse the membrane.

         reaction_sys_vect : list[ReactionABC]
            List of all reaction objects constituting the reaction system.

         rate_base_names : list[tuple[str, Union[float, ndarray]]]
            String names of all reaction rate constants. In the same order as reaction_sys_vect.

         quasi_static_approx : bool
            Use the quasi-static approximation for Vmem?

         ignore_region_updates : Optional[str]
            Ignore updates to concentrations of all chemicals in a region
            indicated by 'i' for 'in' and 'o' for 'out'.

         ignore_chem_updates : Optional[list[Chemical]]
            Ignore updates to concentrations of any chemicals present in this list.

         ediff_reactions_vector : Optional[list[TransportReaction]]
            List of reactions that are strictly electrodiffusion Transport reactions.
        '''

        # Define constants and parameters
        self.R = p.R
        self.F = p.F
        self.T = p.T
        self.eo = p.e_o

        # Typeset labels for chemicals and voltages:
        self._labels = StringNames()

        # Normalization constant for substances for thermodynamic compatibility:
        self.c_norm = 1.0e3
        self.c_norm_s = sp.symbols('c_norm', real=True, positive=True)

        # Define symbolic versions for convenience:
        self.R_s = sp.symbols('R', real=True, positive=True)
        self.F_s = sp.symbols('F', real=True, positive=True)
        self.T_s = sp.symbols('T', real=True, positive=True)

        self.r_cell = p.r_cell
        self.e_mem = p.e_mem
        self.d_mem = p.d_mem
        self.d_ecm = p.d_ecm
        self.r_o = self.r_cell + self.d_ecm

        self.c_mem = (self.eo*self.e_mem)/self.d_mem # Patch capacitance of membrane

        self.div_i = 2/self.r_cell
        self.div_o = (2*self.r_o)/(self.r_o**2 - self.r_cell**2)

        self._quasi_static_approx = quasi_static_approx

        self._ignore_region_tag = ignore_region_updates

        self._ignore_chem_tag = ignore_chem_updates

        self.chem_vect = chem_vect
        self.react_vect = reaction_sys_vect
        self.ediff_vect = ediff_reactions_vector

        # Define the ion objects:
        # For each chemical in the chem_vect, add it as an attribute to this class and
        # initialize its value to what's stored in the chem object.
        for chm in self.chem_vect:
            setattr(self, chm.name, chm.c)
            setattr(self, chm.name + '_s', chm.symbol)

        # Assign the reaction rates as an attribute to this class:
        for rnm, rnv in rate_base_names:
            setattr(self, rnm, rnv)

        # Now we need to scrape the free energy constants from each reaction:
        for react in self.react_vect:
            delg_val = react.deltaGo # Numerical value of delGo for the reaction

            if type(react.deltaGo_s) is not int:
                delg_name = react.deltaGo_s.name
                 # String name of delGo for the reaction
                setattr(self, delg_name, delg_val) # set the attribute to the object

        # self._chem_vals_list = [ci.name for ci in self._chem_vect]
        self._transmem_chem_names = transmem_chem_names

        # Initialize the float values of chemicals to use in simulations:
        self.chem_vals = self._get_param_vals(self.chem_vect)

        self.chem_vect_reduced = []

        for chm in self.chem_vect:
            if ignore_region_updates is not None and self._ignore_chem_tag is None:
                if chm.name.endswith(ignore_region_updates):
                    pass
                else:
                    self.chem_vect_reduced.append(chm)

            elif ignore_region_updates is not None and self._ignore_chem_tag is not None:
                if chm.name.endswith(ignore_region_updates) or chm in self._ignore_chem_tag:
                    pass
                else:
                    self.chem_vect_reduced.append(chm)

            elif self._ignore_chem_tag is not None and ignore_region_updates is None:
                if chm in self._ignore_chem_tag:
                    pass
                else:
                    self.chem_vect_reduced.append(chm)

            else:
                self.chem_vect_reduced.append(chm)

        # Initialize the float values of chemicals to use in simulations:
        self.chem_vals_reduced = self._get_param_vals(self.chem_vect_reduced)

        # Create a Vmem symbol:
        self.V_mem_s = sp.symbols('V_mem', real=True)
        # Parameters to use for v in and vo (two point voltage edge)
        self.V_o_s = sp.symbols('V_o', real=True)
        self.V_i_s = sp.symbols('V_i', real=True)

        # Initialize a dictionary to hold total transmem current for each ion
        self.ion_currents_dict = {}
        self.ediff_currents_dict = {}

        for chm in self._transmem_chem_names:
            self.ion_currents_dict[chm] = 0
            self.ediff_currents_dict[chm] = 0

        # compute currents by adding new terms iteratively:
        for chm in self._transmem_chem_names:
            for react in self.react_vect:
                if chm in react.jc_dict:
                    self.ion_currents_dict[chm] += react.jc_dict[chm]  # Add in the current term for this reaction

            if self.ediff_vect is not None: # Make a current equation just from passive electrodiffusion to calc ghk
                for react in self.ediff_vect:
                    if chm in react.jc_dict:
                        self.ediff_currents_dict[chm] += react.jc_dict[chm]

        self.total_current_s = 0

        for k, v in self.ion_currents_dict.items():
            self.total_current_s += v

        self.ediff_current_s = 0

        if self.ediff_vect is not None:
            for k, v in self.ediff_currents_dict.items():
                self.ediff_current_s += v

            self.ghk_s = sp.solve(self.ediff_current_s, self.V_mem_s)[0]

            # Lambdify (get numerical function) for the GHK equation:
            self.ghk_params_list = list(self.ghk_s.free_symbols)
            self.ghk_f = sp.lambdify(self.ghk_params_list, self.ghk_s)
            # Initialize params:
            self.ghk_params = self._get_param_vals(self.ghk_params_list)

            # Estimate initial Vmem:
            vmem_o = self.ghk_f(*self.ghk_params)

        else:  # Don't create a GHK equation for the system
            self.ghk_params_list = None
            self.ghk_f = None
            self.ghk_s = None
            self.ghk_params = None
            vmem_o = 0.0

        # set the initial value of Vmem:
        self.V_mem = vmem_o
        self.V_o = 0.0
        self.V_i = vmem_o

        if self._quasi_static_approx:
            # Use the quasi-static approximation to calculate Vmem at zero conduction current:
            vsol = sp.solve(self.total_current_s, self.V_mem_s)

            if type(vsol) is list:
                if len(vsol):
                    self.v_ss_s = sp.solve(self.total_current_s, self.V_mem_s)[0]  # Get the analytical expression for steady-state voltage

                else:
                    print("No solution; setting quasi-static option to False.")
                    self._quasi_static_approx = False
                    self.v_ss_s = vsol

            else:
                self.v_ss_s = vsol

        else:
            self.v_ss_s = None

        self.divflux_s = []  # Initialize a master flux vector for the system
        # self.div_s = sp.symbols('div', real=True)
        self.div_i_s = sp.symbols('div_i', real=True) # divergence for inner region
        self.div_o_s = sp.symbols('div_o', real=True) # divergence for outer region
        self.cm_s = sp.symbols('c_mem', real=True, positive=True)

        for chm in self.chem_vect_reduced:  # For each unique chemical in the system...
            flux_chm = 0  # Initialize an object to hold all the flux contributions to this chemical
            for react in self.react_vect:  # For each reaction in the system
                if chm in react.flux_dict:  # Find out if the chemical is involved in the reaction
                    if chm.name.endswith('_i'):
                        flux_chm += self.c_norm_s*self.div_i_s * react.flux_dict[chm]
                    elif chm.name.endswith('_o'):
                        flux_chm += self.c_norm_s*self.div_o_s * react.flux_dict[chm]
            self.divflux_s.append(flux_chm)

        self.divflux_s.append((self.c_norm_s/self.cm_s)*self.total_current_s)  # If we're not using the quasi static approximation, append on the total current vector

        self.divflux_s = sp.Matrix(self.divflux_s).T  # Change into a Sympy matrix so we can lambdify the whole thing

        # Lambdify (get numerical function) for the V_flux vector:
        self.divflux_params_list = list(self.divflux_s.free_symbols)
        self.divflux_f = sp.lambdify(self.divflux_params_list, self.divflux_s)

        # Initialize the float values to use in function call to div_flux_f:
        self.divflux_params = self._get_param_vals(self.divflux_params_list)

        if quasi_static_approx:
            self.divflux_subs_s = self.divflux_s.subs({self.V_mem_s: self.v_ss_s})

            self.divflux_subs_params_list = list(self.divflux_subs_s.free_symbols)

            self.divflux_subs_f = sp.lambdify(self.divflux_subs_params_list, self.divflux_subs_s)

        else:
            self.divflux_subs_s = None
            self.divflux_subs_params_list = None
            self.divflux_subs_f = None

        # Lambdify the equation for total current:
        # Constant parameters that don't change
        self.jc_params_list = list(self.total_current_s.free_symbols)

        self.jc_args_list = np.setdiff1d(np.asarray(self.jc_params_list), np.asarray([self.V_mem_s]),
                                 assume_unique=True).tolist()
        # Value of constants:
        self.jc_arg_vals = self._get_param_vals(self.jc_args_list)

        self.jc_f = sp.lambdify((self.V_mem_s, self.jc_args_list), self.total_current_s)


        # Lambdify the steady-state transmembrane voltage (i.e. use the quasi-steady state
        # approximation of Jc=0) equation:
        if self._quasi_static_approx:
            # Lambdify the steady state voltage function:
            self.vss_params_list = list(self.v_ss_s.free_symbols)
            self.vss_f = sp.lambdify(self.vss_params_list, self.v_ss_s)

            # Initialize the param values to use in function call to Vss_f:
            self.vss_params = self._get_param_vals(self.vss_params_list)

            # Set the V_mem to the appropriate steady-state value as an initial condition:
            # self.V_mem = self.solve_ss_vmem()

        else:
            self.vss_params_list = None
            self.vss_f = None
            self.vss_params = None

        # Create a dictionary of all the permittivity symbols for all reactions indexed by
        # their string name:
        self.P_mem_dict = {}
        for react in self.react_vect:
            self.P_mem_dict[react.P_s.name] = react.P_s


        print("Completed initialization of bioelectrical system.")

    # @beartype
    def _get_param_vals(self, param_names_list: Union[list, ndarray]):
        '''
        Given a list of param names, this method harvests the parameter
        values and returns a list of the values.
        '''

        param_vals_list = []
        for pi in param_names_list:
            vali = getattr(self, pi.name, None)
            if vali is not None:
                # If we've got an attribute, append it to the list:
                param_vals_list.append(vali)
            else:
                raise Exception(f"Parameter {pi.name} not found in BioelectricSystem object!")

        return param_vals_list

    # @beartype
    def _set_param_vals(self,
                        param_names_list,
                        param_names_vals,
                        keep_positive: Optional[list[bool]]=None):
        '''
        Given a list of param names, this method set the parameter
        values to the appropriate variable.
        '''

        if keep_positive is not None:
            # print("Keeping positive values")
            for pi, vi, posi in zip(param_names_list, param_names_vals, keep_positive):
                if posi is True and vi < 0.0: # if less than zero but kept positive, zero the val
                    vi = 0.0
                    # print(f"Value of {pi} has been zeroed {vi}")

                # else:
                    # print(f"Value of {pi} is maintained {vi}")
                setattr(self, pi.name, vi)

        else:
            for pi, vi in zip(param_names_list, param_names_vals):
                setattr(self, pi.name, vi)

    # @beartype
    def _get_chem_vals(self, param_chem_list: Union[list, ndarray]):
        '''
        Given a list of param names, this method harvests the parameter
        values and returns a list of the values.
        '''

        param_vals_list = []
        for pi in param_chem_list:
            vali = getattr(self, pi.name, None)
            if vali is not None:
                # If we've got an attribute, append it to the list:
                param_vals_list.append(vali)
                pi.c = vali
            else:
                raise Exception(f"Parameter {pi.name} not found in BioelectricSystem object!")

        return param_vals_list

    def _set_chem_vals(self):
        '''
        Use values in the chem_vals list to update attributes on the main object.
        '''
        # For each chemical in the system, update it to the new value:
        for chm, chm_val in zip(self.chem_vect_reduced, self.chem_vals_reduced):
            setattr(self, chm.name, chm_val)

        # Synchronize the full chem vector:
        self.chem_vals = self._get_chem_vals(self.chem_vect)

    def opti_jc(self, vm, argv):
        return np.sqrt(self.jc_f(vm, argv) ** 2)

    # @beartype
    def solve_ss_vmem(self, method: str='TNC', force_opti: bool=False):
        '''
        Estimate the steady-state voltage by finding the zero of the current equation.

        Parameters
        -----------
        method : str
            Method to use in optimization

        force_opti : bool
            If True, use Scipy minimize to provide an optimization estimate.


        '''

        if self._quasi_static_approx is False or force_opti is True:
            # Get the present value of the parameters:
            self.jc_arg_vals = self._get_param_vals(self.jc_args_list)

            self.sol_vss = minimize(self.opti_jc, self.V_mem, args=self.jc_arg_vals, method=method)
            vm_ss = self.sol_vss.x[0]

        else:
            self.vss_params_vals = self._get_param_vals(self.vss_params_list)
            vm_ss = self.vss_f(*self.vss_params_vals)

        return vm_ss

    def get_reversal_potentials(self):
        '''
        Calculate the reversal potentials for transmembrane ions.
        '''
        self.v_rev_dict = {}
        self.v_ed_dict = {}

        for ion_base in self._transmem_chem_names:
            c_ion_out = getattr(self, ion_base + '_o', None)
            c_ion_in = getattr(self, ion_base + '_i', None)

            if c_ion_in is not None and c_ion_out is not None:
                # Find the charge of the ion:
                for chm in self.chem_vect:
                    if chm.name.startswith(ion_base):
                        ion_z = chm.z

                        if ion_z == 0 or ion_z == 0.0:
                            v_rev = 0.0
                            v_ed = 0.0

                        else:
                            v_rev = ((self.R * self.T) / (self.F*ion_z)) * np.log(c_ion_out / c_ion_in)
                            v_ed = self.V_mem - v_rev

                        self.v_rev_dict[ion_base] = v_rev
                        self.v_ed_dict[ion_base] = v_ed

                        # set the V_rev attribute to the main object for this:
                        setattr(self, 'V_' + ion_base, v_rev)

        return self.v_rev_dict, self.v_ed_dict

    def print_results(self):
        '''
        Print the value of chemicals and Vmem for the system state.
        '''
        for i, (chm, pval) in enumerate(zip(self.chem_vect, self.chem_vals)):
            print(chm.name, pval)
        print(self.V_mem_s.name, self.V_mem*1e3)

    # @beartype
    def get_pmem_vals(self, ion_perm_list: list[str], resting_pmem_dict: dict) -> dict:
        '''
        Create a dictionary with restime membrane permeability values for the system.
        '''
        # Get the specific permeability variables:
        for ion_name in ion_perm_list:

            pv = getattr(self, ion_name, None)

            if pv is not None:
                resting_pmem_dict[ion_name] = pv

            else:
                raise Exception('Ion membrane permeability not found in bioelectric system!')

        return resting_pmem_dict

    # @beartype
    def set_pmem_vals(self,
                      ion_perm_list: list[str],
                      perm_vals: Union[list, ndarray]):
        '''
        Set the list of permittivity values in self._perm_vals to the bes system.

        Parameters
        ----------
        ion_perm_list : list[str]
            List of attribute names for membrane permeabilities.

        perm_vals: Union[list, ndarray]
            List of values for membrane permeabilities, in same order as ion_perm_list.
        '''

        for pmem_name, pv in zip(ion_perm_list, perm_vals):
            setattr(self, pmem_name, pv)

    def get_jc_hh(self):
        '''
        Compute the Hodgkin-Huxley conduction current for the system.
        '''

        jc = 0.0

        vrev_dict, _ = self.get_reversal_potentials()

        for ion_base_name, ion_vrev in vrev_dict.items():
            Pi = getattr(self, 'P_' + ion_base_name, None)
            c_i = getattr(self, ion_base_name + '_i', None)
            c_o = getattr(self, ion_base_name + '_o', None)


            if Pi is not None and c_i is not None and c_o is not None:
                # Conversion to conductivity from P_mem: (1/d_mem)*c_mem where
                # c_mem is patch capacitance:
                cc = (c_i + c_o) / 2  # average concentration of the ion
                conv_Perm_to_g = (cc*(self.F**2))/(self.R*self.T)
                jc += conv_Perm_to_g*Pi*(self.V_mem - ion_vrev)

            else:
                raise Exception('Ion permeability not found in reaction system!')

        return jc

    def create_network(self, p: ModelParams) -> Dot:
        '''
        Create a plot of the mass transfer network for this reaction system. Returns a
        pydot graphviz object.

        Parameters
        ----------
        p : ModelParams
            An instance of the bioelectric parameters object.

        '''

        # FIXME: to display the pydot graph in Streamlit, we do:
        #  st.graphviz_chart(GG.to_string())

        # create a graph object
        GG = pydot.Dot(
            graph_type='digraph',
            concentrate=False,
            nodesep=0.1,
            ranksep=0.2,
            splines=True,
            strict=True,
            rankdir=p.net_layout,
        )

        # Add nodes to the graph
        for chm in self.chem_vect:
            if chm.z < 0:
                chm_col = p.anion_node_color

            elif chm.z > 0:
                chm_col = p.cation_node_color

            else:
                chm_col = p.neutral_node_color

            ndec = pydot.Node(
                chm.name,
                style='filled',
                color=chm_col,
                shape=p.conc_node_shape,
                fontcolor=p.conc_node_font_color,
                fontname=p.net_font_name,
                fontsize=p.node_font_size,
            )
            GG.add_node(ndec)

        for react in self.react_vect:
            nder = pydot.Node(
                react.reaction_name,
                style='filled',
                color=p.react_node_color,
                shape=p.react_node_shape,
                fontcolor=p.react_node_font_color,
                fontname=p.net_font_name,
                fontsize=p.node_font_size,
            )
            GG.add_node(nder)

        # Add directed edges to the graph:
        for react in self.react_vect:

            for rchm in react._react_list:
                GG.add_edge(
                    pydot.Edge(rchm.name,
                               react.reaction_name,
                               arrowhead='normal',
                               coeff=1.0,
                               penwidth=p.edge_width))

            for pchm in react._prod_list:
                GG.add_edge(
                    pydot.Edge(react.reaction_name,
                               pchm.name, arrowhead='normal',
                               coeff=1.0,
                               penwidth=p.edge_width))

            if react._coupled_reaction is not None:
                # add edges in for the coupled reagents and products:

                for creact in react._react_list_coup:
                    GG.add_edge(pydot.Edge(creact.name,
                                           react.reaction_name,
                                           arrowhead='normal', coeff=1.0,
                                           penwidth=p.edge_width))

                for cprod in react._prod_list_coup:
                    GG.add_edge(pydot.Edge(react.reaction_name,
                                           cprod.name, arrowhead='normal', coeff=1.0,
                                           penwidth=p.edge_width))

        return GG

    def return_elec_props_dict(self):
        '''
        Returns the electrical properties of the system as a Pandas Dataframe
         given a parameters and constants vector.

        '''
        n = self._labels

        vrevs, veds = self.get_reversal_potentials()

        Vmem = self.V_mem  # Vmem
        Na_rev = vrevs['Na'] # Na reversal
        K_rev = vrevs['K'] # K reversal
        Cl_rev = vrevs['Cl'] # Cl reversal
        Na_ed = Vmem - Na_rev # Electrochemical driving force Na
        K_ed = Vmem - K_rev # Electrochemical driving force K
        Cl_ed = Vmem - Cl_rev # Electrochemical driving force Cl

        unit_conv = 1.0e3 # unit conversion to mv

        # Electrical properties dictionary
        elec_props = {n.Vmem: np.round(Vmem*unit_conv, 1),
                      n.Vrev_Na: np.round(Na_rev*unit_conv, 1),
                      n.Vrev_K: np.round(K_rev*unit_conv, 1),
                      n.Vrev_Cl: np.round(Cl_rev*unit_conv, 1),
                      n.Ved_Na: np.round(Na_ed*unit_conv, 1),
                      n.Ved_K: np.round(K_ed*unit_conv, 1),
                      n.Ved_Cl: np.round(Cl_ed*unit_conv, 1)
        }

        elec_dataframe = pd.DataFrame.from_dict(elec_props,
                                                orient='index',
                                                columns=['Voltage [mV]'])

        return elec_dataframe

    def return_chem_props_dict(self, round_dec: int=2):
        '''
        Returns a Pandas dataframe of the ion concentrations inside and outside of the cell
        given a parameters and constants vector.
        '''
        n = self._labels

        steady_state_ions_dict = {
            n.Na_in: round(self.Na_i, round_dec),
            n.K_in: round(self.K_i, round_dec),
            n.Cl_in: round(self.Cl_i, round_dec),
            n.Na_out: round(self.Na_o, round_dec),
            n.K_out: round(self.K_o, round_dec),
            n.Cl_out: round(self.Cl_o, round_dec)
                 }

        ions_dataframe = pd.DataFrame.from_dict(steady_state_ions_dict,
                                                orient='index',
                                                columns=['Concentration [mM]'])

        return ions_dataframe