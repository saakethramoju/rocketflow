from __future__ import annotations

from typing import TYPE_CHECKING

from System import Component, State, Composition

if TYPE_CHECKING:
    from System import Network




class FlowSplitter(Component):

    def __init__(
        self,
        name: str,
        network: Network,
        pressure: State,
        volume: float,
        enthalpy: State | None = None,
        temperature: State | None = None,
        density: State | None = None,
        internal_energy: State | None = None,
        heat_rate: State | float | None = None,
        composition_in: Composition = Composition(),
        composition_out1: Composition = Composition(),
        composition_out2: Composition = Composition(),
        total_enthalpy_in: State | None = None,
        total_enthalpy_out1: State | None = None,
        total_enthalpy_out2: State | None = None,
        mass_flow_in: State | None = None,
        mass_flow_out1: State | None = None,
        mass_flow_out2: State | None = None,
    ):
        self.setup()


        self._solve_energy = self.total_enthalpy_in.is_assigned and self.enthalpy.is_assigned

        if self.enthalpy.is_assigned and not self.total_enthalpy_out1.is_assigned:
            self.total_enthalpy_out1 = self.enthalpy

        if self.enthalpy.is_assigned and not self.total_enthalpy_out2.is_assigned:
            self.total_enthalpy_out2 = self.enthalpy


        self._solve_species = (
            self.composition_in.is_assigned
            and self.composition_out1.is_assigned
        )

        if not self._solve_species:
            if self.composition_in.is_assigned:
                self.composition_out1.copy_from(self.composition_in)
                self.composition_out2.copy_from(self.composition_in)
        else:
            in_species = set(self.composition_in.species)
            out1_species = set(self.composition_out1.species)

            extra_species = out1_species - in_species

            if extra_species:
                raise ValueError(
                    f"{self.name}: composition_out1 contains species not in "
                    f"composition_in: {extra_species}"
                )

            self.composition_out2.copy_from(self.composition_in)


    def evaluate_states(self):
        if not self._solve_species:
            return

        mdot_total = self.mass_flow_in.value
        mdot1 = self.mass_flow_out1.value
        mdot2 = self.mass_flow_out2.value

        if abs(mdot2) < 1e-12:
            raise ValueError(
                f"{self.name}: cannot compute composition_out2 because "
                f"mass_flow_out2 is zero."
            )

        for species, _ in self.composition_in:
            self.composition_out2[species].value = (
                mdot_total * self.composition_in[species].value
                - mdot1 * self.composition_out1[species].value
            ) / mdot2

        self.composition_out2.validate()


    @property
    def iteration_variables(self) -> list[State]:
        variables = [self.pressure]

        if self._solve_energy and self.enthalpy.is_assigned:
            variables.append(self.enthalpy)

        return variables


    @property
    def residuals(self) -> list[float]:
        residuals = [
            self.mass_flow_in.value
            - self.mass_flow_out1.value
            - self.mass_flow_out2.value
        ]

        if self._solve_energy:
            qdot = self.heat_rate.value if self.heat_rate.is_assigned else 0.0

            residuals.append(
                self.mass_flow_in.value * self.total_enthalpy_in.value
                - self.mass_flow_out1.value * self.total_enthalpy_out1.value
                - self.mass_flow_out2.value * self.total_enthalpy_out2.value
                + qdot
            )

        return residuals