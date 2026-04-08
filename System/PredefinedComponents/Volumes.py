from __future__ import annotations

from typing import TYPE_CHECKING

from System import Component, Variable

if TYPE_CHECKING:
    from System import Network, State


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

        self.p = Variable(pressure)
        self.rho = Variable(density)
        self.V = Variable(volume)
        self.mdot_in = Variable(mass_flow_in)
        self.mdot_out = Variable(mass_flow_out)

    def pre_evaluation(self) -> None:
        self.rho.value = 998 * (1 + (self.p.value - 101325) / (2.2e9))

    def evaluate_states(self) -> None:
        pass

    @property
    def iteration_variables(self) -> list[Variable]:
        return [self.p]

    @property
    def residuals(self) -> list[float]:
        return [self.mdot_in.value - self.mdot_out.value]


