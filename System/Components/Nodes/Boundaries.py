from __future__ import annotations

from typing import TYPE_CHECKING

from System import Component, State
from Utilities import Fluid

if TYPE_CHECKING:
    from System import Network




class PressureBoundary(Component):

    def __init__(self, 
                 name: str,
                 network: Network,
                 pressure: State):
        self.initialize_component(name, network)
        self.pressure = pressure

    def evaluate_states(self) -> None:
        pass

    @property
    def iteration_variables(self) -> list[State]:
        return []
    
    @property
    def residuals(self) -> list[float]:
        return []
    

class IsothermalPressureBoundary(Component):

    def __init__(self, 
                 name: str,
                 network: Network,
                 fluid: str,
                 pressure: State,
                 density: State,
                 temperature: State):
        self.initialize_component(name, network)
        self.fluid = fluid
        self.pressure = pressure
        self.rho = density
        self.temperature = temperature

    def pre_evaluation(self) -> None:
        # This is needed since the branches need a density to evaluate their states
        self.rho.value = Fluid(self.fluid, P=self.pressure.value, T=self.temperature.value).density

    def evaluate_states(self) -> None:
        pass

    @property
    def iteration_variables(self) -> list[State]:
        return []
    
    @property
    def residuals(self) -> list[float]:
        return []