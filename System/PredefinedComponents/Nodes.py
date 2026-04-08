from __future__ import annotations

from typing import TYPE_CHECKING

from System import Component, Variable

if TYPE_CHECKING:
    from System import Network, State




class PressureNode(Component):

    def __init__(self, 
                 name: str,
                 network: Network,
                 pressure: State):
        self.initialize_component(name, network)
        self.pressure = Variable(pressure)

    def evaluate_states(self) -> None:
        pass

    @property
    def iteration_variables(self) -> list[Variable]:
        return []
    
    @property
    def residuals(self) -> list[float]:
        return []