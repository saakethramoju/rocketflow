from __future__ import annotations

from typing import TYPE_CHECKING

from System import Component, State
from Utilities import Fluid

if TYPE_CHECKING:
    from System import Network


class DensityfromPT(Component):

    def __init__(self, 
                 name: str,
                 network: Network,
                 fluid: str,
                 pressure: State,
                 temperature: State,
                 density: State):
        self.initialize_component(name, network)
        self.fluid = fluid
        self.pressure = pressure
        self.temperature = temperature
        self.density = density

    def pre_evaluation(self):
        self.evaluate_states()

    def evaluate_states(self) -> None:
        self.density.value = Fluid(self.fluid, P=self.pressure.value, T=self.temperature.value).density

    @property
    def iteration_variables(self) -> list[State]:
        return []
    
    @property
    def residuals(self) -> list[float]:
        return []
    
    '''
    @property
    def density(self) -> State:
        return self.rho
    
    @property
    def pressure(self) -> State:
        return self.p
    
    @property
    def temperature(self) -> State:
        return self.T
    '''