from __future__ import annotations

from typing import TYPE_CHECKING

from System import Component

if TYPE_CHECKING:
    from System import Network, State



class Solid(Component):

    def __init__(
        self,
        name: str,
        network: Network,
        temperature: State,
        mass: float | None = None,
        specific_heat: State | None = None,
        heat_rate: State | float = 0.0,
    ):
        self.setup()

    @property
    def iteration_variables(self) -> list[State]:
        return [self.temperature]

    @property
    def residuals(self) -> list[float]:
        return [self.heat_rate.value]