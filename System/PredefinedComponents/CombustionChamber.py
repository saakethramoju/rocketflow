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
                 fuel: str,
                 oxidizer: str,
                 oxidizer_mass_flow : State,
                 fuel_mass_flow: State, 
                 nozzle_mass_flow: State,
                 characterstic_velocity_efficiency: float):
        self.initialize_component(name, network)

        self.Pc = Variable(chamber_pressure)
        self.fuel = fuel
        self.ox = oxidizer
        self.ox_mdot = Variable(oxidizer_mass_flow)
        self.fuel_mdot = Variable(fuel_mass_flow)
        self.nozzle_mdot = Variable(nozzle_mass_flow)
        self.eta_cstar = Variable(characterstic_velocity_efficiency)

    @property
    def iteration_variables(self) -> list[Variable]:
        return []

    def evaluate_states(self) -> None:
        pass

    @property
    def residuals(self) -> list[float]:
        return []