from __future__ import annotations

from typing import TYPE_CHECKING

from System import Component, State

if TYPE_CHECKING:
    from System import Network




class PressureNode(Component):

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