from __future__ import annotations

import numpy as np
from typing import TYPE_CHECKING

from System import Component, State

if TYPE_CHECKING:
    from System import Network



class DischargeCoefficient(Component):

    def __init__(self,
                 name: str,
                 network: Network,
                 upstream_pressure: State,
                 downstream_pressure: State,
                 density: State,
                 discharge_coefficient: float,
                 cross_sectional_area: float,
                 mass_flow: State | None = None):
        self.setup()

    def evaluate_states(self) -> None:
        P1 = self.upstream_pressure.value
        P2 = self.downstream_pressure.value
        rho = self.density.value
        Cd = self.discharge_coefficient.value
        A = self.cross_sectional_area.value

        self.mass_flow.value = np.sign(P1 - P2) * Cd * A * np.sqrt(2.0 * rho * np.abs(P1 - P2))
