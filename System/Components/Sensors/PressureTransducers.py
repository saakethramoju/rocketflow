from __future__ import annotations

from typing import TYPE_CHECKING

from System import Component

if TYPE_CHECKING:
    from System import Network, State



class PressureTransducer(Component):

    def __init__(self, 
                 name: str, 
                 network: Network,
                 target: State,
                 time_constant: float | None = None,
                 measurement: State | None = None):
        self.setup()

    def evaluate_states(self):
        
        self.measurement.value = self.target.value
