from __future__ import annotations

from typing import TYPE_CHECKING

from System import Component, State

if TYPE_CHECKING:
    from System import Network


class SimpleIncompressibleVolume(Component):

    def __init__(self,
                 name: str,
                 network: Network,
                 pressure: State ,
                 density: State,
                 volume: State,
                 mass_flow_in: State,
                 mass_flow_out: State):
        self.initialize_component(name, network)

        self.p = pressure
        self.rho = density
        self.V = volume
        self.mdot_in = mass_flow_in
        self.mdot_out = mass_flow_out

    def pre_evaluation(self) -> None:
        pass

    def evaluate_states(self) -> None:
        pass

    @property
    def iteration_variables(self) -> list[State]:
        return [self.p]

    @property
    def residuals(self) -> list[float]:
        return [self.mdot_in.value - self.mdot_out.value]


