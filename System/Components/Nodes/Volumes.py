from __future__ import annotations

from typing import TYPE_CHECKING

from System import Component, State

if TYPE_CHECKING:
    from System import Network, Composition



class SimpleVolume(Component):

    def __init__(self,
                 name: str,
                 network: Network,
                 pressure: State ,
                 volume: float,
                 density: State | None = None,
                 temperature: State | None = None,
                 enthalpy: State | None = None,
                 mass_flow_in: State | None = None,
                 mass_flow_out: State | None = None):

        self.setup()

    @property
    def iteration_variables(self) -> list[State]:
        return [self.pressure]

    @property
    def residuals(self) -> list[float]:
        return [self.mass_flow_in.value - self.mass_flow_out.value]


class IsothermalVolume(Component):

    def __init__(self,
                 name:str,
                 network: Network,
                 pressure: State,
                 temperature: State,
                 volume: float,
                 density: State | None = None,
                 enthalpy: State | None = None,
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




class FlowSplitter(Component):
    """
    Assume the flow entering and leaving carries total enthalpy.

    composition is the internal node composition.
    composition_in is the inlet stream composition.
    composition_out1 and composition_out2 are outlet stream compositions.

    If an outlet composition is not provided, it defaults to the node
    composition. Species conservation is only solved when at least one outlet
    composition defaults to the node composition.
    """

    def __init__(
        self,
        name: str,
        network: Network,
        pressure: State,
        enthalpy: State,
        volume: float,
        composition: Composition,
        total_enthalpy_in: State,
        total_enthalpy_out1: State | None = None,
        total_enthalpy_out2: State | None = None,
        composition_in: Composition | None = None,
        composition_out1: Composition | None = None,
        composition_out2: Composition | None = None,
        heat_rate: State | float | None = None,
        temperature: State | None = None,
        density: State | None = None,
        internal_energy: State | None = None,
        mass_flow_in: State | None = None,
        mass_flow_out1: State | None = None,
        mass_flow_out2: State | None = None,
    ):
        self.setup()

        if not self.total_enthalpy_out1.is_assigned:
            self.total_enthalpy_out1 = self.enthalpy

        if not self.total_enthalpy_out2.is_assigned:
            self.total_enthalpy_out2 = self.enthalpy

        if not self.composition_in.is_assigned:
            self.composition_in = self.composition

        out1_uses_node_composition = not self.composition_out1.is_assigned
        out2_uses_node_composition = not self.composition_out2.is_assigned

        if out1_uses_node_composition:
            self.composition_out1 = self.composition

        if out2_uses_node_composition:
            self.composition_out2 = self.composition

        self._solve_species = (
            out1_uses_node_composition
            or out2_uses_node_composition
        )

        if self._solve_species:
            self.composition.constrain_species()

    @property
    def iteration_variables(self) -> list[State]:
        variables = [
            self.pressure,
            self.enthalpy,
        ]

        if self._solve_species:
            variables.extend(self.composition.states[:-1])

        return variables

    @property
    def residuals(self) -> list[float]:
        qdot = self.heat_rate.value if self.heat_rate.is_assigned else 0.0

        energy_in = self.mass_flow_in.value * self.total_enthalpy_in.value

        energy_out = (
            self.mass_flow_out1.value * self.total_enthalpy_out1.value
            + self.mass_flow_out2.value * self.total_enthalpy_out2.value
        )

        residuals = [
            self.mass_flow_in.value
            - (self.mass_flow_out1.value + self.mass_flow_out2.value),

            energy_in + qdot - energy_out,
        ]

        if self._solve_species:
            self.composition.enforce_constraint()

            for species in self.composition.species[:-1]:
                species_in = (
                    self.mass_flow_in.value
                    * self.composition_in[species].value
                )

                species_out = (
                    self.mass_flow_out1.value
                    * self.composition_out1[species].value
                    + self.mass_flow_out2.value
                    * self.composition_out2[species].value
                )

                residuals.append(species_in - species_out)

        return residuals