#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright (c) 2021-2023 Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Input parameter definition and storage class for bioelectric systems modeling.

'''
from beartype import beartype
from beartype.typing import Optional
from dataclasses import dataclass

@beartype
@dataclass
class BioeParams(object):
    '''
    Default parameter initialization and storage class for use in
    bioelectrical system modeling.

    '''

    def __init__(self):

        self.R: float = 8.314  # Ideal gas constant [J/(K*mol)]
        self.F: float = 96485.0  # Faraday's Constant [C/mol]S
        self.e_o: float = 8.3145e-12  # Electric permittivity free space [F/m]
        self.T_C: float = 37.0 # Temp in degrees C

        # constant to convert to normalized concentrations in units mol/L for consistency with thermodynamic params:
        self.c_norm = 1e3

        # UI properties:
        self.use_iterative_solver = False
        self.solve_for_rate_constants = False

        # System properties:
        self.as_vmem: bool = True
        self.update_atp: bool = False
        self.update_extracellular: bool = False
        # Use symmetric flux equations (True) or forward only (False)
        self.use_symmetric_flux: bool = True

        # Linearize equations wrt to Vmem by expressing exp(x) as 1 + x?
        self.linearize_eq: bool = False

        # Auto-solve the constructed bioesystem to steady state?
        self.solve_for_ss: bool = True

        self.dt: float = 5.0e-5 # Time-step for iterative simulations

        # By default, simulations will include Na, K, Cl ions and a Na,K-ATPase ion pump.
        # The following are options that can be included:
        self.use_NaK_ATPase: bool = True
        self.use_NaKCl: bool = False
        self.use_KCl: bool = False

        # Iterative solver settings:
        self.use_hodgkin_huxley: bool = False
        self.clamp_vmem_at: Optional[float] = None

        # Initialize to human blood plasma:
        self.cNa_i: float = 12.0
        self.cNa_o: float = 140.0

        self.cK_i: float = 135.0
        self.cK_o: float = 4.0

        self.cCl_i: float = 4.0
        self.cCl_o: float = 116.0

        self.cATP: float = 4.0
        self.cADP: float = 0.1
        self.cPi: float = 0.1

        self.delG_ATP = -32e3  # Gibbs free energy for ATP hydrolysis in J

        # Membrane permittivities and reaction rate constants:
        self.base_pmem: float = 1.0e-9

        self.pump_unit_modifier = 1.0 # 5.0e-6

        self.base_PNa: float = 1.0
        self.base_PK: float = 25.0
        self.base_PCl: float = 5.0

        self.base_NaKpump: float = 1.0e4
        self.base_NaKCl: float = 2.0e4
        self.base_KCl: float = 1.0e4

        self.r_cell_um: float = 7.5
        self.d_mem: float = 5.0e-9
        self.d_ecm: float = 25e-9
        self.e_mem: float = 10.0

        self.vmemo = 0.0  # initial Vmem in mV

        # Network plotting options:--------------------------------------
        self.net_font_name = 'DejaVu Sans'
        self.node_font_size = 16
        self.tit_font_size = 24
        self.net_layout = 'TB'
        self.edge_width = 2.0

        self.conc_node_font_color = 'Black'
        self.conc_node_color = 'PaleTurquoise' #'LightCyan'
        self.anion_node_color = 'PaleTurquoise'
        self.cation_node_color = 'LightSalmon'
        self.neutral_node_color = 'Silver'
        self.conc_node_shape = 'ellipse'

        self.react_node_font_color = 'White'
        self.react_node_color = 'Gunmetal' # 'DarkSeaGreen'
        self.react_node_shape = 'rect'

        self.transp_node_font_color = 'White'
        self.transp_node_color = 'DarkSeaGreen'
        self.transp_node_shape = 'diamond'

        self.chan_node_font_color = 'White'
        self.chan_node_color = 'DarkSeaGreen'
        self.chan_node_shape = 'pentagon'


        self.update_parameters()

    def update_parameters(self):
        '''
        Recalculate calculated params after updating some settings.
        '''

        self.T = self.T_C + 273.15 # Temp in Kelvin
        self.cmem = (self.e_o * self.e_mem) / self.d_mem  # membrane capacitance
        self.r_cell = self.r_cell_um*1.0e-6

        # Electrodiffusion permittivity constants:
        self.PNa = self.base_pmem * self.base_PNa
        self.PK = self.base_pmem * self.base_PK
        self.PCl = self.base_pmem * self.base_PCl

        # Ion pump rate constants:
        self.PNaK_ATPase = self.base_pmem * self.base_NaKpump * self.pump_unit_modifier

        # Transporter rate constants:
        self.PNaKCl = self.base_pmem * self.base_NaKCl * self.pump_unit_modifier
        self.PKCl = self.base_pmem * self.base_KCl * self.pump_unit_modifier