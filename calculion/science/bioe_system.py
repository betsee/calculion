#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright (c) 2022-2023 Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Create a bioelectric system for a single cell

'''
from beartype import beartype
from beartype.typing import Optional
from calculion.science.model_params import ModelParams
from calculion.science.chem_base import (Chemical,
                                            TransportReaction,
                                            ClassicReaction,
                                            CoupledReaction)
from calculion.science.reaction_system import ReactionSystem

class BioElectricSystem(object):
    '''
    Create a bioelectric system for modeling a single cell.

    Public Attributes
    ------------------
    c_norm : float
        Normalizing constant for the reaction, where concentrations are converted into mol/L
        from mol/m3 for thermodynamic purposes.

    c_norm_s : Symbol
        Symbolic normalizing constant for the reaction.

    chem_vals : list
        List of values of all Chemicals in the system.

    chem_vals_reduced : list
        List of values of Chemicals that have updated concentrations.

    chem_vect : list
        List of all Chemicals in the system.

    chem_vect_reduced : list
        List of Chemicals that have updated concentrations.

    divflux_f : function
        The 'change vector' function, computing all concentration fluxes and
        the conduction current.

    divflux_params : list
        Parameter values for the divflux_f function.

    divflux_params_list : list
        Parameter names for the divflux_f function.

    divflux_s : MutableDenseMatrix
        Symbolic 'change vector' expression, computing all concentration fluxes and
        the conduction current.

    ion_currents_dict : dict
        Dictionary of ion current expressions indexed to each transmembrane ion base name.

    jc_arg_vals : list
        Constant parameter values for the jc_f function.

    jc_args_list : list
        Constant parameter names for the jc_f function.

    jc_f : function
        The numerical conduction current function for the system.

    jc_params_list : list
        List of parameters used in the jc_f function.

    react_vect : list
        List of all ReactionABC objects in the system.

    total_current_s : Add
        Symbolic expression of total conduction current for the system.

    Private Attibutes
    ------------------
    _ignore_chem_tag : list
        List of chemicals to ignore updates.

    _ignore_region_tag : str
        List of regions where chemical concentrations are not updated.

    _quasi_static_approx : bool
        Use the quasi-static approximation for V_mem (True)?

    _transmem_chem_names : list
        Names of substances that traverse the membrane.

    '''

    # @beartype
    def __init__(self,
                 p: ModelParams,
                 irrev_fwd_only: bool = True):
        '''
        Initialize the Bioelectrical System.

        Parameters
        ------------
        p : Optional[ModelParams]
            An instance of BioeParams object.

        irrev_fwd_only : bool
            Do ion pumps and reactions with large negative standard free energy change
            (like ATP hydrolysis) only have a forward flux direction?

        '''

        if p is None:
            # Then make a default parameters object:
            p = ModelParams()

        self._p = p

        if irrev_fwd_only:
            self._flux_sym_irrev = False
        else:
            self._flux_sym_irrev = p.use_symmetric_flux

        self._build_system(p)

    def _build_system(self, p: ModelParams):
        '''
        Build the bioelectrical system.

        Parameters
        ------------
        p : Optional[ModelParams]
            An instance of BioeParams object.

        '''

        self._ion_base_names = ['Na', 'K', 'Cl']
        self.rate_name_opti_list = ['P_Na', 'P_K', 'P_Cl']

        # Define the ion objects:
        Na_i = Chemical('Na_i', 1, p.cNa_i)
        Na_o = Chemical('Na_o', 1, p.cNa_o)
        K_i = Chemical('K_i', 1, p.cK_i)
        K_o = Chemical('K_o', 1, p.cK_o)
        Cl_i = Chemical('Cl_i', -1, p.cCl_i)
        Cl_o = Chemical('Cl_o', -1, p.cCl_o)

        # And the metabolic objects:
        atp = Chemical('ATP', -4, p.cATP)
        adp = Chemical('ADP', -3, p.cADP)
        pi = Chemical('P', -1, p.cPi)

        # List of all chemical agents defined for this bioelectricity model
        self._all_chem = [Na_i, Na_o, K_i, K_o, Cl_i, Cl_o, atp, adp, pi]

        # Save main ions for use elsewhere:
        self._Na_i = Na_i
        self._Na_o = Na_o
        self._K_i = K_i
        self._K_o = K_o
        self._Cl_i = Cl_i
        self._Cl_o = Cl_o

        # Define ion passive transport reactions:
        ediff_Na = TransportReaction([Na_i], [1], [Na_o], [1],
                                         react_pos=[0], prod_pos=[1],
                                         deltaGo=0.0,
                                         rate_const=p.PNa,
                                         write_as_vmem=p.as_vmem,
                                         base_names=['Na'],
                                         rate_base_name='P_Na',
                                         T=p.T,
                                         use_symmetric_flux=p.use_symmetric_flux,
                                         linearize_eq=p.linearize_eq,
                                         reaction_name='Na Channel')

        self.ediff_Na = ediff_Na

        ediff_K = TransportReaction([K_i], [1], [K_o], [1],
                                        react_pos=[0], prod_pos=[1],
                                        deltaGo=0.0,
                                        rate_const=p.PK,
                                        write_as_vmem=p.as_vmem,
                                        base_names=['K'],
                                        rate_base_name='P_K',
                                        T=p.T,
                                        use_symmetric_flux=p.use_symmetric_flux,
                                        linearize_eq=p.linearize_eq,
                                        reaction_name='K Channel')


        ediff_Cl = TransportReaction([Cl_i], [1], [Cl_o], [1],
                                         react_pos=[0], prod_pos=[1],
                                         deltaGo=0.0,
                                         rate_const=p.PCl,
                                         write_as_vmem=p.as_vmem,
                                         base_names=['Cl'],
                                         rate_base_name='P_Cl',
                                         T=p.T,
                                         use_symmetric_flux=p.use_symmetric_flux,
                                         linearize_eq=p.linearize_eq,
                                         reaction_name='Cl Channel')

        # Initialize reactions vector:
        self._all_reactions = [ediff_Na, ediff_K, ediff_Cl]

        self._all_ediff_reactions = [ediff_Na, ediff_K, ediff_Cl] # create a separate vector only for electrodiffusion

        # Initialize permeability/rates vector:
        self._all_rates = [('P_Na', p.PNa), ('P_K', p.PK), ('P_Cl', p.PCl)]

        # Define the ATP hydrolysis reaction:
        # react_ATP_hydr = ion.ClassicReaction([atp], [1], [adp, pi], [1, 1], deltaGo = delG_ATP)
        self._react_ATP_hydr = CoupledReaction([atp], [1],
                                             [adp, pi], [1, 1],
                                             p.delG_ATP,
                                             deltaGo_base_name='delG_ATP')

        # For the sake of chemical reaction study, save a classic chemical reaction corresponding to ATP hydrolysis:
        self.ATP_hydrolysis_react = ClassicReaction([atp], [1], [adp, pi], [1, 1],
                                                      deltaGo = p.delG_ATP,
                                                      rate_const = p.PATP,
                                                      T=p.T,
                                                      rate_base_name='P_ATP',
                                                      use_symmetric_flux = p.use_symmetric_flux,
                                                      reaction_name='ATP Hydrolysis')

        # Inhibit updates to concentration of ATP, ADP, and P:
        if p.update_atp is False:
            self._ignore_chem_updates = [atp, adp, pi]
        else:
            self._ignore_chem_updates = None
            self.rate_name_opti_list.append('P_ATP')

        if p.update_extracellular is False:
            self._ignore_region_updates = '_o'
        else:
            self._ignore_region_updates = None

        if p.use_NaK_ATPase:

            # And the transport reaction for the Na, K ATPase pump's ion trasport coupled to an ATP hydrolysis reaction:
            transp_NaKpump = TransportReaction([Na_i, K_o], [3, 2],
                                                   [Na_o, K_i], [3, 2],
                                                   react_pos=[0, 1],
                                                   prod_pos=[1, 0],
                                                   deltaGo=0.0,
                                                   rate_const=p.PNaK_ATPase,
                                                   base_names=['Na', 'K'],
                                                   coupled_reaction=self._react_ATP_hydr,
                                                   write_as_vmem=True,
                                                   rate_base_name='P_NaKATP',
                                                   T=p.T,
                                                   # use_symmetric_flux=False,
                                                   use_symmetric_flux=self._flux_sym_irrev,
                                                   linearize_eq=p.linearize_eq,
                                                   reaction_name='Na/K-ATPase')

            self.transp_NaKpump = transp_NaKpump  # save for example case

            self._all_reactions.append(transp_NaKpump)
            self._all_rates.append(('P_NaKATP', p.PNaK_ATPase))
            self.rate_name_opti_list.append('P_NaKATP')


        if p.use_NaKCl:
            # Define an Na-K-2Cl cotransporter:
            # And the transport reaction for the Na, K ATPase pump's ion trasport coupled to an ATP hydrolysis reaction:
            transp_NaKCl = TransportReaction([Na_o, K_o, Cl_o], [1, 1, 2],
                                                 [Na_i, K_i, Cl_i], [1, 1, 2],
                                                 react_pos=[1, 1, 1],
                                                 prod_pos=[0, 0, 0],
                                                 deltaGo=0.0,
                                                 rate_const=p.PNaKCl,
                                                 base_names=['Na', 'K', 'Cl'],
                                                 coupled_reaction=None,
                                                 write_as_vmem=True,
                                                 rate_base_name='P_NaKCl',
                                                 T=p.T,
                                                 use_symmetric_flux=p.use_symmetric_flux,
                                                 linearize_eq=p.linearize_eq,
                                                 reaction_name='Na/K/Cl-cotrans')

            self.transp_NaKCl = transp_NaKCl # Save to study

            self._all_reactions.append(transp_NaKCl)
            self._all_rates.append(('P_NaKCl', p.PNaKCl))
            self.rate_name_opti_list.append('P_NaKCl')

        if p.use_KCl:

            transp_KCl = TransportReaction([K_i, Cl_i], [1, 1],
                                          [K_o, Cl_o], [1, 1],
                                          react_pos=[0, 0],
                                          prod_pos=[1, 1],
                                          deltaGo=0.0,
                                          rate_const=p.PKCl,
                                          base_names = ['K', 'Cl'],
                                          coupled_reaction = None,
                                          write_as_vmem=True,
                                          rate_base_name = 'P_KCl',
                                          T=p.T,
                                          use_symmetric_flux = p.use_symmetric_flux,
                                          linearize_eq=p.linearize_eq,
                                          reaction_name='K/Cl-Cotrans')

            self._all_reactions.append(transp_KCl)
            self._all_rates.append(('P_KCl', p.PKCl))
            self.rate_name_opti_list.append('P_KCl')

        # BUILD REACTION SYSTEM--------------------------------------------

        # Create a bioelectric system object based on the desired attributes:
        self.bes = ReactionSystem(p,
                                 self._all_chem,
                                 self._ion_base_names,
                                 self._all_reactions,
                                 self._all_rates,
                                 quasi_static_approx=False,
                                 ignore_region_updates=self._ignore_region_updates,
                                 ignore_chem_updates=self._ignore_chem_updates
                                 )
