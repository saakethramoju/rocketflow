from __future__ import annotations

from typing import TYPE_CHECKING

from System import Component, State

if TYPE_CHECKING:
    from System import Network, Composition


class SimpleVolume(Component):

    def __init__(
        self,
        name: str,
        network: Network,
        pressure: State,
        volume: float,
        density: State | None = None,
        temperature: State | None = None,
        enthalpy: State | None = None,
        composition: Composition | None = None,
        composition_in: Composition | None = None,
        mass_flow_in: State | None = None,
        mass_flow_out: State | None = None,
    ):
        self.setup()

        if self.composition is not None:
            self.composition.constrain_species()

    @property
    def iteration_variables(self) -> list[State]:
        return [
            self.pressure,
            *(self.composition[species] for species in self.composition & self.composition_in)
        ]

    @property
    def residuals(self) -> list[float]:
        self.composition.enforce_constraint()

        return [
            self.mass_flow_in.value - self.mass_flow_out.value,
            *(self.mass_flow_in.value * self.composition_in[species].value
                - self.mass_flow_out.value * self.composition[species].value
                for species in self.composition & self.composition_in
            ),
        ]




class IsothermalVolume(Component):

    def __init__(self,
                 name:str,
                 network: Network,
                 pressure: State,
                 temperature: State,
                 volume: float,
                 density: State | None = None,
                 enthalpy: State | None = None,
                 composition: Composition | None = None,
                 mass_flow_in: State | None = None,
                 mass_flow_out: State | None = None):
        
        self.setup()
    
    @property
    def iteration_variables(self) -> list[State]:
        return [self.pressure]

    @property
    def residuals(self) -> list[float]:
        return [self.mass_flow_in.value - self.mass_flow_out.value]



class Volume(Component):
    def __init__(
        self,
        name: str,
        network: Network,
        pressure: State,
        enthalpy: State,
        volume: float,
        total_enthalpy_in: State,
        total_enthalpy_out: State | None = None,
        heat_rate: State | float | None = None,
        temperature: State | None = None,
        density: State | None = None,
        internal_energy: State | None = None,
        composition: Composition | None = None,
        mass_flow_in: State | None = None,
        mass_flow_out: State | None = None,
    ):
        self.setup()

        if not self.total_enthalpy_out.is_assigned:
            self.total_enthalpy_out = self.enthalpy

    @property
    def iteration_variables(self) -> list[State]:
        return [self.pressure, self.enthalpy]

    @property
    def residuals(self) -> list[float]:
        qdot = self.heat_rate.value if self.heat_rate.is_assigned else 0.0

        return [
            self.mass_flow_in.value - self.mass_flow_out.value,
            (
                self.mass_flow_in.value * self.total_enthalpy_in.value
                - self.mass_flow_out.value * self.total_enthalpy_out.value
                + qdot
            ),
        ]