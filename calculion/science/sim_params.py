#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright (c) 2022-2025 Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Input parameter definition and storage class for simulation settings of
bioelectric system models.
'''

from dataclasses import dataclass

@dataclass
class SimParams(object):
    '''
    Default parameter initialization and storage class for use in
    bioelectrical system modeling.
    '''

    def __init__(self):

        # UI properties:-----------------------------------------
        self.solve_for_rate_constants: bool = False
        self.solve_for_steady_state: bool = False
        self.use_iterative_solver: bool = False

        # Iterative solver default properties:
        self.delta_t: float = 1.0e-3  # Iterative solver time step
        # self.dt: float = 5.0e-5 # Time-step for iterative simulations
        self.start_time: float = 0.5  # Start time for plotting; always starts at 0.0
        self.end_time: float = 4.5

        # Membrane perturbations:
        self.perturb_PNa: bool = True
        self.perturb_PNa_start: float = 1.20
        self.perturb_PNa_end: float = 1.50
        self.perturb_PNa_multi: float = 10.0

        self.perturb_PK: bool = True
        self.perturb_PK_start: float = 2.20
        self.perturb_PK_end: float = 2.50
        self.perturb_PK_multi: float = 10.0

        self.perturb_PCl: bool = True
        self.perturb_PCl_start: float = 3.20
        self.perturb_PCl_end: float = 3.50
        self.perturb_PCl_multi: float = 10.0

        # Generate a comparison simulation with HH math?
        self.use_hh_math: bool = False

        self.update_parameters()


    def update_parameters(self):

        # Iterative solver calculated properties:
        self.starttime_plot_ind = int(self.start_time / self.delta_t)
        self.N_iter = int(self.end_time / self.delta_t)
