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
                 pressure: State,
                 temperature: State,):
        self.initialize_component(name, network)
        self.p = pressure
        self.T = temperature

    def evaluate_states(self) -> None:
        pass

    @property
    def iteration_variables(self) -> list[State]:
        return []
    
    @property
    def residuals(self) -> list[float]:
        return []