from __future__ import annotations

from typing import TYPE_CHECKING

from System import Component, Variable

if TYPE_CHECKING:
    from System import Network, State


class RocketCEACombustionChamber(Component):

    def __init__(self, 
                 name: str, 
                 network: Network,
                 chamber_pressure: State,
                 oxidizer_mass_flow : State,
                 fuel_mass_flow: State, 
                 nozzle_mass_flow: State,):
        self.initialize_component(name, network)

        self.Pc = Variable(chamber_pressure)
        self.ox_mdot = Variable(oxidizer_mass_flow)
        self.fuel_mdot = Variable(fuel_mass_flow)
        self.nozzle_mdot = Variable(nozzle_mass_flow)

    @property
    def iteration_variables(self) -> list[Variable]:
        return [self.Pc]

    def evaluate_states(self) -> None:
        pass

    @property
    def residuals(self) -> list[float]:
        return [self.fuel_mdot.value + self.ox_mdot.value - self.nozzle_mdot.value]