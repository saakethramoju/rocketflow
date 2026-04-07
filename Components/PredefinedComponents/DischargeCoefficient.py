from __future__ import annotations

import numpy as np
from typing import TYPE_CHECKING

from Components import Component, Variable

if TYPE_CHECKING:
    from System import State


class DischargeCoefficient(Component):

    def __init__(self,
                 name: str,
                 upstream_pressure: State | float,
                 downstream_pressure: State | float,
                 density: State | float,
                 discharge_coefficient: float,
                 cross_sectional_area: float,
                 mass_flow: State | float):
        super().__init__(name)

        self.upstream_pressure = Variable(upstream_pressure)
        self.downstream_pressure = Variable(downstream_pressure)
        self.density = Variable(density)
        self.mass_flow = Variable(mass_flow)
        self.Cd = Variable(discharge_coefficient)
        self.A = Variable(cross_sectional_area)

    def evaluate_states(self) -> None:
        P1 = self.upstream_pressure.value
        P2 = self.downstream_pressure.value
        rho = self.density.value
        Cd = self.Cd.value
        A = self.A.value

        self.mass_flow.value = (
            np.sign(P1 - P2) * Cd * A * np.sqrt(2.0 * rho * np.abs(P1 - P2))
        )

    @property
    def iteration_variables(self) -> list[Variable]:
        return []

    @property
    def residuals(self) -> list[float]:
        return []