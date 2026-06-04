from __future__ import annotations

from typing import TYPE_CHECKING

from System import Component

if TYPE_CHECKING:
    from System import Network, State


class SolidConductor(Component):

    def __init__(self, 
                 name: str, 
                 network: Network,
                 upstream_temperature: State,
                 downstream_temperature: State,
                 thermal_conductivity: State,
                 length: float,
                 cross_sectional_area: float,
                 heat_rate: State | None = None):
        self.setup()

    def evaluate_states(self):
        k = self.thermal_conductivity.value
        A = self.cross_sectional_area.value
        l = self.length.value
        T1 = self.upstream_temperature.value
        T2 = self.downstream_temperature.value

        self.heat_rate.value = k*A/l * (T1 - T2)