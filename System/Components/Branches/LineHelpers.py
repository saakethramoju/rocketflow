from __future__ import annotations

from typing import TYPE_CHECKING

from System import Component, State

if TYPE_CHECKING:
    from System import Network


class GravityPressureChange(Component):

    def __init__(
        self,
        name: str,
        network: Network,
        upstream_pressure: State | float,
        density: State | float,
        elevation_change: State | float,
        gravitional_acceleration: State | float = 9.80665,
        downstream_pressure: State | None = None,
    ):
        """
        Elevation change is positive upwards
        """
        self.setup()

        if downstream_pressure is None:
            self.downstream_pressure = State()

    def evaluate_states(self) -> None:
        self.downstream_pressure.value = (
            self.upstream_pressure.value
            - self.density.value
            * self.gravitional_acceleration.value
            * self.elevation_change.value
        )