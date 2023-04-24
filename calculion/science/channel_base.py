#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright (c) 2022-2023 Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
This module creates different ion channel classes to use in Calculion.

'''

from abc import ABCMeta, abstractmethod
from beartype import beartype
from beartype.typing import Union, Optional
import numpy as np
from numpy import ndarray

class DynamicChannelABC(object, metaclass=ABCMeta):
    '''
    Abstract base class of all dynamic channel classes.

    Attributes
    ----------
    '''

    @beartype
    def __init__(self,
                 ion_permeability_list: list[str],
                 Pmax_chan: Union[float, ndarray],
                 vm: Union[float, ndarray]
                 ):
        '''

        '''

        # List of ion permeabilities affected by channel (as bes var name in BES)
        self.ion_perm_list = ion_permeability_list

        self.set_channel_params()

        self.Pmax_chan = Pmax_chan # scale ion channel perm

        self.init_channel(vm)

    @abstractmethod
    def set_channel_params(self):
        '''
        Sets the unique parameters for the ion channel.
        '''

        pass



    @beartype
    def init_channel(self, vm: Union[float, ndarray]):
        '''
        Runs the initialization sequence for a voltage gated ion channel to set
        the channel to steady-state at the requested vm.
        '''

        self.update_mh_gates(vm)

        self.m = self.m_inf
        self.h = self.h_inf

        self.P_chan = self.Pmax_chan*(self.m**self.n_m)*(self.h**self.n_h)

    @beartype
    def run_channel(self, vm: Union[float, ndarray], dt: float, t: float):
        '''
        Runs the voltage gated ion channel.
        '''

        self.update_mh_gates(vm)

        # Update channel state using semi-Implicit Euler method:-------------------
        dt = dt*self.time_unit

        self.m = (self.m_tau*self.m + dt*self.m_inf)/(self.m_tau + dt)
        self.h = (self.h_tau*self.h + dt*self.h_inf)/(self.h_tau + dt)

        self.P_chan = self.Pmax_chan*(self.m**self.n_m)*(self.h**self.n_h)

        return self.P_chan

    @abstractmethod
    @beartype
    def update_mh_gates(self, vm: Union[float, ndarray]):
        """
        Updates the 'm' and 'h' gating functions of the channel model for
        standard Hodgkin-Huxley formalism channel models.

        """
        pass

class StepFunctionChannel(DynamicChannelABC):

    def __init__(self,
                ion_permeability_list: list[str],
                Pmax_chan: float,
                T_chan: float,
                t_delay: float = 0.0,
                t_end: float = 1.0e15):

        dummy_vmem = 1.0 # We don't need Vmem for the step channel

        super().__init__(ion_permeability_list, Pmax_chan, dummy_vmem)

        self._T_chan = T_chan
        self._tdelay = t_delay
        self._tend = t_end

    def set_channel_params(self):
        self.pmem_scale = 1.0

    def init_channel(self, vm):
        pass

    def update_mh_gates(self, vm: Union[float, ndarray]):
        pass

    @beartype
    def run_channel(self,
                    vm: Union[float, ndarray],
                    dt: float,
                    t: float):
        '''
        Function determining change to membrane permeability.
        '''

        if t >= self._tdelay and t < self._tend:
            Pmem_mod = (self.Pmax_chan/2) * (1 + np.sign((np.sin((2 * np.pi * t) / self._T_chan))))
            # Pmem_mod = (self.Pmax_chan / 2) * (1 + (np.sin((2 * np.pi * t) / self._T_chan)))
        else:
            if type(vm) is ndarray:
                Pmem_mod = np.zeros(len(vm))

            else:
                Pmem_mod = 0.0

        return Pmem_mod

class PulseFunctionChannel(DynamicChannelABC):

    def __init__(self,
                ion_permeability_list: list[str],
                Pmax_chan: float,
                t_start: float,
                t_stop: float
                ):

        dummy_vmem = 1.0 # We don't need Vmem for the step channel

        super().__init__(ion_permeability_list, Pmax_chan, dummy_vmem)

        self._tstart = t_start
        self._tstop = t_stop

    def set_channel_params(self):
        self.pmem_scale = 1.0

    def init_channel(self, vm):
        pass

    def update_mh_gates(self, vm: Union[float, ndarray]):
        pass

    @beartype
    def run_channel(self,
                    vm: Union[float, ndarray],
                    dt: float,
                    t: float):
        '''
        Function determining change to membrane permeability.
        '''

        if t >= self._tstart and t < self._tstop:
            Pmem_mod = self.Pmax_chan

        else:
            if type(vm) is ndarray:
                Pmem_mod = np.zeros(len(vm))

            else:
                Pmem_mod = 0.0

        return Pmem_mod