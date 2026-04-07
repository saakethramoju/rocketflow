from __future__ import annotations

from typing import TYPE_CHECKING
from Components import Component, Variable


if TYPE_CHECKING:
    from System import State


class PressureNode(Component):

    def __init__(self, 
                 name: str,
                 pressure: State):
        super().__init__(name)
        self.pressure = Variable(pressure)

    def evaluate_states(self) -> None:
        pass

    @property
    def iteration_variables(self) -> list[Variable]:
        return []
    
    @property
    def residuals(self) -> list[float]:
        return []