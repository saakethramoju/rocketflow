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




class SeriesCdA(Component):

    def __init__(
        self,
        name: str,
        network: Network,
        effective_areas: list[State | float],
        effective_area: State | None = None,
    ):
        """
        Equivalent effective area for flow restrictions in series.

        Uses:
            1 / CdA_eq^2 = sum(1 / CdA_i^2)
        """
        self.setup()

    def evaluate_states(self):
        inverse_area_squared_sum = 0.0

        for effective_area in self.effective_areas:
            CdA = effective_area.value

            if abs(CdA) < 1e-12:
                self.effective_area.value = 0.0
                return

            inverse_area_squared_sum += 1.0 / CdA**2

        self.effective_area.value = 1.0 / inverse_area_squared_sum**0.5



class ParallelCdA(Component):

    def __init__(
        self,
        name: str,
        network: Network,
        effective_areas: list[State | float],
        effective_area: State | None = None,
    ):
        """
        Equivalent effective area for flow restrictions in parallel.

        Uses:
            CdA_eq = sum(CdA_i)
        """
        self.setup()

    def evaluate_states(self):
        self.effective_area.value = sum(
            effective_area.value
            for effective_area in self.effective_areas
        )
