#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright (c) 2022-2023 Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
This module creates systems for working with ions and ionic reactions in IonSpire.

'''

from abc import ABCMeta, abstractmethod
from beartype import beartype
from beartype.typing import Optional, Union
import sympy as sp
from sympy.core.symbol import Symbol
import numpy as np
from numpy import ndarray
from calculion.science.chem_enum import ReactionClass, SystemDimension


@beartype
class Chemical(object):
    '''
    An object to store all the data associated with a typical ionic substance.

    Attributes
    ---------
    name: str
        Name of the chemical.
    z : int
        Charge of the chemical (ion).
    c: float or ndarray
        Concentration of the chemical as a single volume entity (float) or on an
        array of spatial locations.
    symbol: Symbol
        Sympy symbol object representing the chemical in analytic expressions.

    '''

    def __init__(self, name: str, z: int, c: Union[float, ndarray]):
        '''
        Initialize the ion object.

        Parameters
        ------------
        name: str
            Name of the chemical.
        z : int
            Charge of the chemical (ion).
        c: float or ndarray
            Concentration of the chemical as a single volume entity (float) or on an
            array of spatial locations.

        '''
        self.name = name
        self.symbol = sp.symbols(name, real=True, positive=True) # create a Sympy symbol for the name
        self.z = z # Charge of the ion
        self.c = c # Initial concentration of the ion
@beartype
class CoupledReaction(object):
    '''
    Class to define a reaction that will be coupled to another reaction, for example, this
    might define an ATP hydrolysis reaction that couples with a transporter reaction to
    define an ion pump.

    Public Attributes
    -----------------
    deltaGo : Union[int, float]
        The standard Gibbs Free Energy of the reaction in J/mol.

    deltaGo_s : Symbol
        A Sympy Symbol representing the standard Gibbs Free energy of the reaction, which is
        used in analytic sympy expressions.

    Private Attributes
    -----------------
    _react_list : list[Chemical]
        A list of the Chemical objects that are the reactants in the reaction.

    _react_stoic : list[Union[float, int]]
        A list of floats representing the stoichiometry of each reactant in the reaction.
        Must be in the same order as the react_list Chemicals.

    _prod_list : list[Chemical]
        A list of the Chemical objects that are the products in the reaction.

    _prod_stoic : list[Union[float, int]]
        A list of floats representing the stoichiometry of each product in the reaction.
        Must be in the same order as the prod_list Chemicals.

    '''

    def __init__(self,
                 react_list: list[Chemical],
                 react_stoic: list[Union[float, int]],
                 prod_list: list[Chemical],
                 prod_stoic: list[Union[float, int]],
                 deltaGo: Union[int, float] = 0.0,
                 deltaGo_base_name: Optional[str] = None):
        '''
        Initialize the reaction object.

        Parameters
        ----------
        react_list : list[Chemical]
            A list of the Chemical objects that are the reactants in the reaction.

        react_stoic : list[Union[float, int]]
            A list of floats representing the stoichiometry of each reactant in the reaction.
            Must be in the same order as the react_list Chemicals.

        prod_list : list[Chemical]
            A list of the Chemical objects that are the products in the reaction.

        prod_stoic : list[Union[float, int]]
            A list of floats representing the stoichiometry of each product in the reaction.
            Must be in the same order as the prod_list Chemicals.

        deltaGo : Union[int, float]
            The standard Gibbs Free Energy of the reaction in J/mol.

        deltaGo_base_name : Optional[str]
            A string name representing the standard Gibbs Free energy of the reaction, which is
            used in analytic sympy expressions.

        '''

        self._react_list = react_list
        self._react_stoic = react_stoic
        self._prod_list = prod_list
        self._prod_stoic = prod_stoic
        self.deltaGo = deltaGo

        if deltaGo_base_name is None:
            self.deltaGo_s = sp.symbols('deltaGo_coupled', real=True)

        else:
            self.deltaGo_s = sp.symbols(deltaGo_base_name, real=True)

@beartype
class ReactionABC(object, metaclass=ABCMeta):
    '''
    Base class for Reactions.

    Public Attributes
    -----------------
    deltaGo : Union[int, float]
        The standard Gibbs Free Energy of the reaction in J/mol.

    rate_const : Union[int, float]
        The rate constant for the reaction.

    reaction_name : Optional[str]
        The name of the reaction. This is used in network graphs.

    T : Union[int, float]
        The temperature for the reaction, in degrees K.

    R : float
        Ideal gas constant, R=8.314 J/(K mol)

    F : float
        Planck's constant, F=96485 C/mol

    P_s : Symbol
        Sympy Symbol for the rate constant. This is used as a symbol in analytic equations.

    R_s : Symbol
        Sympy Symbol for the ideal gas constant. This is used as a symbol in analytic equations.

    F_s : Symbol
        Sympy Symbol for Planck's constant. This is used as a symbol in analytic equations.

    T_s : Symbol
        Sympy Symbol for temperature. This is used as a symbol in analytic equations.

    V_s : Symbol
        Sympy Symbol for voltage. This is used as a symbol in analytic equations.

    Vo_s : Symbol
        Sympy Symbol for voltage inside cell. This is used as a symbol in analytic equations.

    Vi_s : Symbol
        Sympy Symbol for voltage outside cell. This is used as a symbol in analytic equations.

    Vm_s : Symbol
        Sympy Symbol for transmembrane voltage. This is used as a symbol in analytic equations.

    Keqm_s : Symbol
        Sympy Symbol for reaction eqm constant. This is used as a symbol in analytic equations.

    Keqm : Union[int, float]
        The equilibrium constant of the reaction.

    deltaGo_s : Symbol
        Sympy Symbol for reaction standard free energy. This is used as a symbol in analytic equations.


    Private Attributes
    -----------------
    _react_list : list[Chemical]
        A list of the Chemical objects that are the reactants in the reaction.

    _react_stoic : list[Union[float, int]]
        A list of floats representing the stoichiometry of each reactant in the reaction.
        Must be in the same order as the react_list Chemicals.

    _prod_list : list[Chemical]
        A list of the Chemical objects that are the products in the reaction.

    _prod_stoic : list[Union[float, int]]
        A list of floats representing the stoichiometry of each product in the reaction.
        Must be in the same order as the prod_list Chemicals.

    _sys_dim : SystemDimension
        Dimension of the system.

    _coupled_reaction : bool
        Coupled reaction used?

    _react_list_coup : list[Chemical]
        List of reactants in coupled reaction.

    _react_stoic_coup : list[Union[float, int]]
        Reaction stoichiometry of coupled reaction reactants.

    _prod_list_coup : list[Chemical]
        List of products in coupled reaction.

    _prod_stoic_coup : list[Union[float, int]]
        Reaction stoichiometry of coupled reaction products.

    '''

    def __init__(self,
                 react_list: list[Chemical],
                 react_stoic: list[Union[float, int]],
                 prod_list: list[Chemical],
                 prod_stoic: list[Union[float, int]],
                 deltaGo: Union[int, float] = 0.0,
                 rate_const: Union[int, float] = 1.0,
                 T: Union[int, float] = 310,
                 rate_base_name: Optional[str] = None,
                 system_dimension: SystemDimension = SystemDimension.d2,
                 reaction_name: Optional[str]=None):
        '''
        Initialize the Reaction object.

        Parameters
        ----------
        react_list : list[Chemical]
            A list of the Chemical objects that are the reactants in the reaction.

        react_stoic : list[Union[float, int]]
            A list of floats representing the stoichiometry of each reactant in the reaction.
            Must be in the same order as the react_list Chemicals.

        prod_list : list[Chemical]
            A list of the Chemical objects that are the products in the reaction.

        prod_stoic : list[Union[float, int]]
            A list of floats representing the stoichiometry of each product in the reaction.
            Must be in the same order as the prod_list Chemicals.

        deltaGo : Union[int, float]
            The standard Gibbs Free Energy of the reaction in J/mol.

        rate_const : Union[int, float]
            The rate constant for the reaction.

        T : Union[int, float]
            The temperature for the reaction, in degrees K.

        rate_base_name : Optional[str]
            String name for the rate constant. This is used as a symbol in analytic equations.

        system_dimension : SystemDimension
            Dimension of the system.

        reaction_name : Optional[str]
            The name of the reaction. This is used in network graphs.

        '''

        self._sys_dim = system_dimension # Remember the system dimension
        self.reaction_name = reaction_name # Remember name of the reaction

        # Symbolic variable of factor that normalizes concentrations:
        self._cnorm_s = sp.symbols('c_norm', real=True, positive=True)

        # Initialize various terms for the reaction:
        # Symbolic constants
        self.R_s = sp.symbols('R', real=True, positive=True)
        self.F_s = sp.symbols('F', real=True, positive=True)
        self.T_s = sp.symbols('T', real=True, positive=True)

        # Numerical constants
        self.R = 8.314
        self.F = 96485
        self.T = T

        self.rate_const = rate_const

        if len(react_list) != len(react_stoic):
            raise Exception("Length of reactants list doesn't equalreactant stoichiometry!")

        else:
            self._react_list = react_list
            self._react_stoic = react_stoic

        if len(prod_list) != len(prod_stoic):
            raise Exception("Length of products list doesn't equal product stoichiometry!")

        else:
            self._prod_list = prod_list
            self._prod_stoic = prod_stoic

        if rate_base_name is None: # If no symbolic name is given for rate base name
            self.P_s = sp.symbols('P_r', real=True, positive=True)
            setattr(self, 'P_r', rate_const)

        else:
            self.P_s = sp.symbols(rate_base_name, real=True, positive=True)
            setattr(self, rate_base_name, rate_const) # Set the rate constant to the object

        # Voltages:
        self.V_s = sp.symbols('V', real=True)  # symbolic voltage
        self.Vo_s = sp.symbols('V_o', real=True)  # symbolic voltage outside compartment
        self.Vi_s = sp.symbols('V_i', real=True)  # symbolic voltage inside compartment
        self.Vm_s = sp.symbols('V_mem', real=True)  # symbolic transmem voltage (Vi - Vo)

        if deltaGo == 0.0:
            self.Keqm_s = 1
            self.deltaGo = 0.0
            self.deltaGo_s = 0

        else:
            self.Keqm_s = sp.symbols('K_eqm', real=True, positive=True)
            self.Keqm = np.exp(-deltaGo / (self.R * self.T))
            self.deltaGo = deltaGo
            self.deltaGo_s = sp.symbols('deltaGo', real=True)

        # set coupled reaction attributes to None as an initialization:
        self._coupled_reaction = None
        self._react_list_coup = None
        self._react_stoic_coup = None
        self._prod_list_coup = None
        self._prod_stoic_coup = None

    @abstractmethod
    def _compute_flux_eqn(self):
        '''

        '''
        pass

    @abstractmethod
    def _write_reaction_G(self):
        '''

        '''
        pass

    # @beartype
    def get_numerical(self, sympy_express):
        '''
        Convert an analytic sympy expression to a numerical numpy function with variable names list and
        initial variables.
        '''
        # Lambdify (get numerical function) for the V_flux vector:
        express_params_list = list(sympy_express.free_symbols)
        numpy_express = sp.lambdify(express_params_list, sympy_express)

        return numpy_express, express_params_list

@beartype
class TransportReaction(ReactionABC):
    '''
    Define a chemical reaction characterized by movement of ions across a membrane
    or from one compartment to another. This type of reaction can be used to define
    an electrodiffusive flux, a transporter, and when coupled to an ATP hydrolysis
    reaction, an ion pump.

    Public Attributes
    ------------------
    F : int
        Planck's constant (F=96485 C/mol)

    F_s : Symbol
        Symbolic Planck's constant.

    Keqm_s : int
        Symbolic reaction equilibrium constant.

    P_s : Symbol
        Symbolic reaction rate constant.

    Q_denom_s : Mul
        Symbolic reaction quotient denominator.

    Q_numer_s : Mul
        Symbolic reaction quotient numerator.

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

    deltaG_s : Add
        Symbolic reaction Gibbs Free Energy change.

    deltaGo : float
        Reaction Gibbs Free Energy change (J/mol).

    deltaGo_s : Symbol
        Symbolic reaction Gibbs Free Energy change.

    flux_dict : dict
        Dictionary of analytic average flux expessions indexed by each reactant or product.

    flux_f_dict : dict
        Dictionary of analytic forward flux expessions indexed by each reactant or product.

    flux_r_dict : dict
        Dictionary of analytic reverse flux expessions indexed by each reactant or product.

    flux_tm_dict : dict
        Dictionary specifying analytic transmembrane average flux expression for each
    transmembrane ion.

    flux_tm_f_dict : dict
        Dictionary specifying analytic transmembrane forward flux expression for each
    transmembrane ion.

    flux_tm_r_dict : dict
        Dictionary specifying analytic transmembrane reverse flux expression for each
    transmembrane ion.

    jc_dict : dict
        Dictionary specifying analytic transmembrane current expression for each
    transmembrane ion.

    rate_const : float
        Value of the reaction rate constant.

    react_flx_f_s : Mul
        Symbolic reaction forward flux expression.

    react_flx_r_s : Mul
        Symbolic reaction reverse flux expression.

    react_flx_s : Add
        Symbolic reaction flux expression.

    reaction_name : str
        Name of the reaction.

    sys_agents : list
        A list of all chemicals (i.e. both reactants and products) involved in the reaction.

    term_exp : Mul
        Term of the reaciton raised to an exponential function.

    Private Attributes
    ---------------------
    _base_names : list
        List of the base name strings of ions, omitting the location tags
    (i.e. 'Na_i', 'Na_o' base name is 'Na')

    _cnorm_s : Symbol
        Normalizing constant for the reaction, where concentrations are converted into mol/L
    from mol/m3.

    _coupled_reaction : NoneType or CoupledReaction
        A Classic reaction that is energetically coupled to drive this transport reaction.

    _linearize_eq : bool
        Linearize the exponentials of expressions using exp(x) ~ 1 + x (assumes small x)?

    _prod_list : list
        List of product chemicals of the reaction with concentrations that cross the membrane.

    _prod_list_full : list
        List of all product chemicals of the reaction, including coupled reaction products from
        a Coupled Reaction.

    _prod_pos : list
        List of position tags for each product that transverses the membrane.

    _prod_stoic : list
        List of stoichiometry coefficients for each product that crosses the membrane.

    _prod_stoic_full : list
        List of all product stoichimetry coefficients, including coupled reaction products from
    a Coupled Reaction.

    _react_list : list
        List of reactant chemicals of the reaction with concentrations that cross the membrane.

    _react_list_full : list
        List of all reactant chemicals of the reaction, including coupled reaction products from
    a Coupled Reaction.

    _react_pos : list
        List of position tags for each reactant that transverses the membrane.

    _react_stoic : list
        List of stoichiometry coefficients for each reactant that crosses the membrane.

    _react_stoic_full : list
        List of all reactant stoichimetry coefficients, including coupled reaction products from
    a Coupled Reaction.

    _reaction_type : ReactionClass
        The class of the reaction (always TransportReaction).

    _sys_dim : SystemDimension
        Dimention of the system.

    _use_symmetric_flux : bool
        Use the flux expressions averaging forward and backward flux reactions (True)
        or only foward flux (False)?

    _write_as_vmem : bool
        Subsitute Vmem = Vi - Vo into all analytic expressions (True) or leave as-is (False)?

    '''
    def __init__(self,
                 react_list: list[Chemical],
                 react_stoic: list[Union[float, int]],
                 prod_list: list[Chemical],
                 prod_stoic: list[Union[float, int]],
                 react_pos: Optional[list[int]] = None,
                 prod_pos: Optional[list[int]] = None,
                 deltaGo: Union[int, float] = 0.0,
                 rate_const: Union[int, float] = 1.0,
                 T: Union[int, float]=310,
                 write_as_vmem: bool=True,
                 base_names: Optional[list[str]]= None,
                 rate_base_name: Optional[str] = None,
                 coupled_reaction: Optional[CoupledReaction] = None,
                 use_symmetric_flux: bool = True,
                 linearize_eq: bool = False,
                 reaction_name: Optional[str] = None
                 ):
        '''
        Initialize the TransportReaction object.

        Parameters
        ----------
        react_list : list[Chemical]
            A list of the Chemical objects that are the reactants in the reaction.

        react_stoic : list[Union[float, int]]
            A list of floats representing the stoichiometry of each reactant in the reaction.
        Must be in the same order as the react_list Chemicals.

        react_pos : list[int]
            A list of position tag integers specifying where each reactant is. Here
        0 indicates outside the cell, while 1 indicates inside the cell, thereby dictating
        the direction of the gradient and positive flux into the cell.
        Must be in the same order as the react_list Chemicals.

        prod_list : list[Chemical]
            A list of the Chemical objects that are the products in the reaction.

        prod_stoic : list[Union[float, int]]
            A list of floats representing the stoichiometry of each product in the reaction.
        Must be in the same order as the prod_list Chemicals.

        prod_pos : list[int]
            A list of position tag integers specifying where each product is. Here
        0 indicates outside the cell, while 1 indicates inside the cell, thereby dictating
        the direction of the gradient and positive flux into the cell.
        Must be in the same order as the prod_list Chemicals.

        deltaGo : Union[int, float]
            The standard Gibbs Free Energy of the reaction in J/mol.

        rate_const : Union[int, float]
            The rate constant for the reaction.

        T : Union[int, float]
            The temperature for the reaction, in degrees K.

        rate_base_name : Optional[str]
            String name for the rate constant. This is used as a symbol in analytic equations.

        reaction_name : Optional[str]
            The name of the reaction. This is used in network graphs.

        write_as_vmem: bool
            Re-express the interior and exterior voltage symbols as a Vmem = Vi - Vo?

        base_names: Optional[list[str]]
            Base name of ions transversing the membrane. For example, if 'Na_i' and 'Na_o' are
        the symbols for sodium transversing the membrane, then the base name 'Na' would be
        included in this list.

        coupled_reaction: Optional[CoupledReaction]
            Object instance of a classic reaction coupled to this transport reaction. For example,
        a coupled_reaction may be a ClassicReaction defining ATP hydrolysis.

        use_symmetric_flux : bool
            Incorporate forward and backward reaction expressions into the expression for flux?

        linearize_eq : bool
            Use the small-argument linearization for exponentials (i.e. exp(x) ~ 1 + x) to
            linearize the equations?

        '''

        self._reaction_type = ReactionClass.transport

        self._use_symmetric_flux = use_symmetric_flux

        self._linearize_eq = linearize_eq

        # Initialize the base class
        super().__init__(react_list,
                       react_stoic,
                       prod_list,
                       prod_stoic,
                       deltaGo=deltaGo,
                       rate_const=rate_const,
                       T=T,
                       rate_base_name=rate_base_name,
                       reaction_name=reaction_name)

        self._coupled_reaction = coupled_reaction

        # Save the position indicators use in transport reactions:
        self._react_pos = react_pos
        self._prod_pos = prod_pos

        self._write_as_vmem = write_as_vmem

        if base_names is None:
            self._base_names = [str(i) for i in range(len(react_list))]
        else:
            self._base_names = base_names

        # Write the symbolic delta_G expression for the reaction:
        self._write_reaction_G()
        self._compute_flux_eqn()
        self._compute_current_eqn()

    def _write_reaction_G(self):
        '''
        Calculate the Gibbs Free energy change (delta_G) of the reaction as an analytical
        Sympy expression.
        '''

        if self._coupled_reaction is not None:
            self._process_coupled_reaction()

        else:
            self.coupled_deltaG_s = 0
            self._react_list_full = self._react_list
            self._react_stoic_full = self._react_stoic
            self._prod_list_full = self._prod_list
            self._prod_stoic_full = self._prod_stoic
            self.coupled_Q_numer_s = 1
            self.coupled_Q_denom_s = 1

        # When working with a transport reaction, the voltages are positional:
        term_r = 0
        term_q_denominator = 0
        term_exp = 0
        for ci, pi, stoi in zip(self._react_list, self._react_pos, self._react_stoic):

            if pi == 0:  # If position tag is 0
                Vp = self.Vo_s  # Voltage is 'out'
            else:
                Vp = self.Vi_s  # Voltage is 'in'

            # Add to the reaction coefficient denominator:
            term_q_denominator += sp.log((ci.symbol/self._cnorm_s)**stoi)

            term_exp += stoi * ci.z * self.F_s * Vp

            term_r += stoi * self.R_s * self.T_s * sp.log(ci.symbol/self._cnorm_s) + stoi * ci.z * self.F_s * Vp

        term_p = 0
        term_q_numerator = 0
        for ci, pi, stoi in zip(self._prod_list, self._prod_pos, self._prod_stoic):

            if pi == 0:  # If position tag is 0
                Vp = self.Vo_s  # Voltage is 'out'
            else:
                Vp = self.Vi_s  # Voltage is 'in'

            # Add to the reaction coefficient numerator:
            term_q_numerator += sp.log((ci.symbol/self._cnorm_s)**stoi)
            term_p += stoi * self.R_s * self.T_s * sp.log(ci.symbol/self._cnorm_s) + stoi * ci.z * self.F_s * Vp

            term_exp += -stoi * ci.z * self.F_s * Vp

        term_eq = (term_p.simplify() -
                   term_r.simplify() +
                   self.deltaGo_s).simplify()

        term_exp += self.deltaGo_s # add on the standard free energy of reaction
        self.term_exp = term_exp

        self.deltaG_s = term_eq # assign the term to the equation
        self.Q_numer_s = sp.exp(sp.logcombine(term_q_numerator)) # Save the final Q numerator
        self.Q_denom_s = sp.exp(sp.logcombine(term_q_denominator)) # Save the final Q denominator

        if self._write_as_vmem: # If user wants things in terms of Vmem instead of Vi Vo
            self._as_vmem() # Substitute in Vmem = Vi - Vo

    def _as_vmem(self):
        '''
        process the symbolic delta G expression to be in terms of Vmem
        '''
        # Process the reaction term to be in Vmem:
        self.deltaG_s = self.deltaG_s.subs({self.Vi_s: self.Vm_s + self.Vo_s,
                                    self.Vo_s: -self.Vm_s + self.Vi_s}).simplify()
        self.term_exp = self.term_exp.subs({self.Vi_s: self.Vm_s + self.Vo_s,
                                            self.Vo_s: -self.Vm_s + self.Vi_s}).simplify()

    def _compute_flux_eqn(self):
        '''
        Compute analytic equations describing the transmembrane/transcompartment
        flux of each ion in a transport equation.
        '''

        # Write the flux equations for the transport equation of an ion or passive transporter:
        self.flux_dict = {}
        self.flux_f_dict = {}
        self.flux_r_dict = {}
        self.flux_tm_dict = {}
        self.flux_tm_f_dict = {}
        self.flux_tm_r_dict = {}

        self.sys_agents = [] # List storing the chemical entities involved for each flux

        # Forward reaction:
        if self._linearize_eq is False:
            exp_f = sp.exp(-(1/(self.R_s*self.T_s))*self.term_exp)
            exp_r = sp.exp((1/(self.R_s*self.T_s))*self.term_exp)

        else:
            exp_f = 1 -(1/(self.R_s*self.T_s))*self.term_exp
            exp_r = 1 + (1/(self.R_s*self.T_s))*self.term_exp

        self.react_flx_f_s = self.P_s*(self.Q_numer_s*self.coupled_Q_numer_s -
                                       self.Q_denom_s*self.coupled_Q_denom_s*exp_f)

        # Reverse reaction:
        self.react_flx_r_s = self.P_s*(self.Q_denom_s*self.coupled_Q_denom_s -
                                       self.Q_numer_s*self.coupled_Q_numer_s*exp_r)

        # True reaction rate is average of forward and revered reverse rates:
        self.react_flx_s = (1/2)*(self.react_flx_f_s - self.react_flx_r_s)

        # Compute the transmembrane fluxes:
        for ri, pi, rsi, psi, rpos, ppos, bn in zip(self._react_list,
                                  self._prod_list,
                                  self._react_stoic,
                                  self._prod_stoic,
                                  self._react_pos,
                                  self._prod_pos,
                                  self._base_names
                                  ):

            flx_f = rsi*self.react_flx_f_s
            flx_r = psi*self.react_flx_r_s
            flx = rsi*self.react_flx_s

            # Save the asymmetric fluxes:
            self.flux_f_dict[ri] = flx_f
            self.flux_f_dict[pi] = -flx_f

            self.flux_r_dict[ri] = -flx_r
            self.flux_r_dict[pi] = flx_r

            # Symmetric fluxes
            if self._use_symmetric_flux:
                self.flux_dict[ri] = flx
                self.flux_dict[pi] = -flx
            else:
                self.flux_dict[ri] = flx_f
                self.flux_dict[pi] = -flx_f

            self.sys_agents.append(ri)
            self.sys_agents.append(pi)

            # Calculate the sign change for movement from one compartment to another:
            pos_change = np.sign(ppos - rpos)
            self.flux_tm_dict[bn] = pos_change*flx
            self.flux_tm_f_dict[bn] = pos_change*flx_f
            self.flux_tm_r_dict[bn] = pos_change*flx_r

        # If there's a coupled reaction, add in the standard fluxes (non transmem) for
        # this whole:
        if self._coupled_reaction is not None:
            for ri, rsi in zip(self._react_list_coup,
                               self._react_stoic_coup):

                flx = rsi * self.react_flx_s

                self.flux_dict[ri] = flx
                self.sys_agents.append(ri)

            for pi, psi in zip(self._prod_list_coup, self._prod_stoic_coup):
                # flx = self.P_s*(psi * pi.symbol -
                #                 psi * sp.solve(self.deltaG_s, pi.symbol)[0].simplify())
                flx = -psi * self.react_flx_s

                self.flux_dict[pi] = flx
                self.sys_agents.append(pi)

    def _process_coupled_reaction(self):
        '''
        Process the coupled reaction.
        '''

        coupled_react_list = self._coupled_reaction._react_list
        coupled_react_stoic = self._coupled_reaction._react_stoic
        coupled_prod_list = self._coupled_reaction._prod_list
        coupled_prod_stoic = self._coupled_reaction._prod_stoic
        coupled_deltaGo_s = self._coupled_reaction.deltaGo_s
        coupled_deltaGo = self._coupled_reaction.deltaGo

        coupled_term_r = 0  # Free energy reactants term
        term_q_denominator = 0
        for rea, stoi in zip(coupled_react_list, coupled_react_stoic):
            coupled_term_r += stoi * self.R_s * self.T_s * sp.log(rea.symbol)
            term_q_denominator += sp.log(rea.symbol**stoi)

        coupled_term_p = 0  # Free energy products term
        term_q_numerator = 0
        for prod, stoi in zip(coupled_prod_list, coupled_prod_stoic):
            coupled_term_p += stoi * self.R_s * self.T_s * sp.log(prod.symbol)
            term_q_numerator += sp.log(prod.symbol ** stoi)

        coupled_term_eq = (coupled_term_p.simplify() -
                           coupled_term_r.simplify() +
                           coupled_deltaGo_s).simplify()

        self.coupled_Q_numer_s = sp.exp(sp.logcombine(term_q_numerator)) # Save the final Q numerator
        self.coupled_Q_denom_s = sp.exp(sp.logcombine(term_q_denominator)) # Save the

        self.coupled_deltaG_s = coupled_term_eq
        # self.deltaG_s = self.deltaG_s + self.coupled_deltaG_s
        self.deltaGo_s = coupled_deltaGo_s # override the original symbol with coupled
        self.deltaGo += coupled_deltaGo

        self.coupled_term_exp_s = sp.exp((1 / (self.R_s * self.T_s)) * self.deltaGo_s)

        # Define values for the full list of reactants and products and their
        # stoichiometry:
        self._react_list_full = self._react_list + coupled_react_list
        self._react_stoic_full = self._react_stoic + coupled_react_stoic
        self._prod_list_full = self._prod_list + coupled_prod_list
        self._prod_stoic_full = self._prod_stoic + coupled_prod_stoic

        self._react_list_coup = coupled_react_list
        self._react_stoic_coup = coupled_react_stoic
        self._prod_list_coup = coupled_prod_list
        self._prod_stoic_coup = coupled_prod_stoic

    def _compute_current_eqn(self):
        '''
        Compute analytic equations describing the transmembrane conduction
        current (jc) of each ion in a transport equation.
        '''

        self.jc_dict = {} # Transmembrane conduction current for each ion

        # Compute the transmembrane fluxes:
        for ri, pi, rsi, psi, bn in zip(self._react_list,
                                    self._prod_list,
                                    self._react_stoic,
                                    self._prod_stoic,
                                    self._base_names
                                    ):

            if self._use_symmetric_flux:
                jc = ri.z*self.F_s*self.flux_tm_dict[bn]

            else:
                jc = ri.z*self.F_s*self.flux_tm_f_dict[bn]

            # Symmetric currents wrt the ion base name:
            self.jc_dict[bn] = jc

@beartype
class ClassicReaction(ReactionABC):
    '''
    Define a chemical reaction characterized by consumption and production of
    chemical agents in one spatial location.

        Public Attributes
    ------------------
    F : int
        Planck's constant (F=96485 C/mol)

    F_s : Symbol
        Symbolic Planck's constant.

    Keqm_s : Add
        Symbolic reaction equilibrium constant.

    P_s : Symbol
        Symbolic reaction rate constant.

    Q_denom_s : Mul
        Symbolic reaction quotient denominator.

    Q_numer_s : Mul
        Symbolic reaction quotient numerator.

    R : float
        Ideal gas constant (R=8.314 J/(K mol))

    R_s : Symbol
        Symbolic ideal gas constant.

    T : int
        Reaction temperature (degrees K).

    T_s : Symbol
        Symbolic reaction temperature.

    deltaG_s : Add
        Symbolic reaction Gibbs Free Energy change.

    deltaGo : float
        Reaction Gibbs Free Energy change (J/mol).

    deltaGo_s : Symbol
        Symbolic reaction Gibbs Free Energy change.

    flux_dict : dict
        Dictionary of analytic average flux expessions indexed by each reactant or product.

    flux_f_dict : dict
        Dictionary of analytic forward flux expessions indexed by each reactant or product.

    flux_r_dict : dict
        Dictionary of analytic reverse flux expessions indexed by each reactant or product.

    rate_const : float
        Value of the reaction rate constant.

    react_flx_f_s : Mul
        Symbolic reaction forward flux expression.

    react_flx_r_s : Mul
        Symbolic reaction reverse flux expression.

    react_flx_s : Add
        Symbolic reaction flux expression.

    reaction_name : str
        Name of the reaction.

    sys_agents : list
        A list of all chemicals (i.e. both reactants and products) involved in the reaction.


    Private Attributes
    --------------------
    _cnorm_s : Symbol
        Normalizing constant for the reaction, where concentrations are converted into mol/L
    from mol/m3.

    _prod_list : list
        List of product chemicals of the reaction with concentrations that cross the membrane.

    _prod_stoic : list
        List of stoichiometry coefficients for each product that crosses the membrane.

    _react_list : list
        List of reactant chemicals of the reaction with concentrations that cross the membrane.

    _react_stoic : list
        List of stoichiometry coefficients for each reactant that crosses the membrane.

    _reaction_type : ReactionClass
        The class of the reaction (always TransportReaction).

    _sys_dim : SystemDimension
        Dimention of the system.

    _use_symmetric_flux : bool
        Use the flux expressions averaging forward and backward flux reactions (True)
        or only foward flux (False)?

    '''

    def __init__(self,
                 react_list: list[Chemical],
                 react_stoic: list[Union[float, int]],
                 prod_list: list[Chemical],
                 prod_stoic: list[Union[float, int]],
                 deltaGo: Union[int, float] = 0.0,
                 rate_const: Union[int, float] = 1.0,
                 T: Union[int, float]=310,
                 rate_base_name: Optional[str] = None,
                 use_symmetric_flux: bool = True,
                 reaction_name: Optional[str] = None
                ):

        '''
        Initialize the ClassicReaction object.

                Parameters
        ----------
        react_list : list[Chemical]
            A list of the Chemical objects that are the reactants in the reaction.

        react_stoic : list[Union[float, int]]
            A list of floats representing the stoichiometry of each reactant in the reaction.
            Must be in the same order as the react_list Chemicals.

        prod_list : list[Chemical]
            A list of the Chemical objects that are the products in the reaction.

        prod_stoic : list[Union[float, int]]
            A list of floats representing the stoichiometry of each product in the reaction.
            Must be in the same order as the prod_list Chemicals.

        deltaGo : Union[int, float]
            The standard Gibbs Free Energy of the reaction in J/mol.

        rate_const : Union[int, float]
            The rate constant for the reaction.

        T : Union[int, float]
            The temperature for the reaction, in degrees K.

        rate_base_name : Optional[str]
            String name for the rate constant. This is used as a symbol in analytic equations.

        reaction_name : Optional[str]
            The name of the reaction. This is used in network graphs.

        use_symmetric_flux : bool
            Incorporate forward and backward reaction expressions into the expression for flux?

        '''

        self._reaction_type = ReactionClass.classic

        # Initialize the base class
        super().__init__(react_list,
                       react_stoic,
                       prod_list,
                       prod_stoic,
                       deltaGo=deltaGo,
                       rate_const=rate_const,
                       T=T,
                       rate_base_name=rate_base_name,
                       reaction_name=reaction_name)

        self._use_symmetric_flux = use_symmetric_flux

        # Write the symbolic delta_G expression for the reaction:
        self._write_reaction_G()
        self._compute_flux_eqn()

        # A classic reaction won't have current as it's not transmembrane, but add an empty
        # currents dictionary to be consistent with transport reactions.:
        self.jc_dict = {}


    def _write_reaction_G(self):
        '''
        Calculate the Gibbs Free energy change (delta_G) of the reaction as an analytical
        Sympy expression.
        '''

        term_r = 0 # Free energy reactants term
        term_q_denominator = 0
        for rea, stoi in zip(self._react_list, self._react_stoic):
            term_r += stoi * self.R_s * self.T_s * sp.log(rea.symbol/self._cnorm_s) + stoi * rea.z * self.F_s * self.V_s

            # Add to the reaction coefficient numerator:
            term_q_denominator += sp.log((rea.symbol/self._cnorm_s)**stoi)

        term_p = 0 # Free energy products term
        term_q_numerator = 0
        for prod, stoi in zip(self._prod_list, self._prod_stoic):
            term_p += stoi * self.R_s * self.T_s * sp.log(prod.symbol/self._cnorm_s) + prod.z * self.F_s * self.V_s * stoi

            # Add to the reaction coefficient numerator:
            term_q_numerator += sp.log((prod.symbol/self._cnorm_s)**stoi)

        term_eq = (term_p.simplify() -
                   term_r.simplify() +
                   self.deltaGo_s).simplify()

        self.deltaG_s = term_eq
        self.Q_numer_s = sp.exp(sp.logcombine(term_q_numerator)) # Save the final Q numerator
        self.Q_denom_s = sp.exp(sp.logcombine(term_q_denominator)) # Save the

    def _compute_flux_eqn(self):
        '''
        Compute analytic equations describing the flux of each ion in a transport
        equation.

        delG_drive : Add
            A Sympy equation describing the Gibbs free energy change of a driving
            expression, such as ATP hydrolysis.
        '''

        # Write the flux equations for the classic reaction:
        self.flux_dict = {}
        self.sys_agents = []

        # Forward reaction:
        self.react_flx_f_s = self.P_s * (self.Q_numer_s -
                                       self.Q_denom_s * sp.exp(-(self.deltaGo_s / (self.R_s * self.T_s)))
                                       )
        # Reverse reaction:
        self.react_flx_r_s = self.P_s * (self.Q_denom_s -
                                       self.Q_numer_s * sp.exp((self.deltaGo_s / (self.R_s * self.T_s)))
                                       )

        # True reaction rate is average of forward and reversed reverse rates:
        self.react_flx_s = (1/2)*(self.react_flx_f_s - self.react_flx_r_s)

        for ri, rsi in zip(self._react_list,
                           self._react_stoic):
            # flx = -self.P_s * (rsi * ri.symbol -
            #                    rsi * sp.solve(self.deltaG_s, ri.symbol)[0].simplify())

            if self._use_symmetric_flux:
                flx = rsi*self.react_flx_s
            else:
                flx = rsi*self.react_flx_f_s

            self.flux_dict[ri] = flx
            self.sys_agents.append(ri)

        for pi, psi in zip(self._prod_list, self._prod_stoic):
            # flx = self.P_s * (psi * pi.symbol -
            #                   psi * sp.solve(self.deltaG_s, pi.symbol)[0].simplify())
            flx = -psi * self.react_flx_s
            self.flux_dict[pi] = -flx
            self.sys_agents.append(pi)

def get_numerical(sympy_express):
    '''
    Convert an analytic sympy expression to a numerical numpy function with variable names list and
    initial variables.
    '''
    # Lambdify (get numerical function) for the V_flux vector:
    express_params_list = list(sympy_express.free_symbols)
    numpy_express = sp.lambdify(express_params_list, sympy_express)

    return numpy_express, express_params_list