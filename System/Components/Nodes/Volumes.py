from __future__ import annotations

from typing import TYPE_CHECKING

from System import Component, State
from Utilities import Fluid

if TYPE_CHECKING:
    from System import Network

class IsothermalIncompressibleVolume(Component):

    def __init__(self,
                 name:str,
                 network: Network,
                 fluid: str,
                 pressure: State,
                 temperature: State,
                 density: State,
                 volume: float,
                 mass_flow_in: State,
                 mass_flow_out: State):
        self.initialize_component(name, network)

        self.fluid = fluid
        self.p = pressure
        self.T = temperature
        self.rho = density
        self.V = State(volume)
        self.mdot_in = mass_flow_in
        self.mdot_out = mass_flow_out

    def pre_evaluation(self) -> None:
        # This is needed since the branches need a density to evaluate their states
        self.evaluate_states()

    def evaluate_states(self) -> None:
        self.rho.value = Fluid(self.fluid, P=self.p.value, T=self.T.value).density

    @property
    def iteration_variables(self) -> list[State]:
        return [self.p]

    @property
    def residuals(self) -> list[float]:
        return [self.mdot_in.value - self.mdot_out.value]



class SimpleIncompressibleVolume(Component):

    def __init__(self,
                 name: str,
                 network: Network,
                 pressure: State ,
                 density: State,
                 volume: float,
                 mass_flow_in: State,
                 mass_flow_out: State):
        self.initialize_component(name, network)

        self.p = pressure
        self.rho = density
        self.V = State(volume)
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


