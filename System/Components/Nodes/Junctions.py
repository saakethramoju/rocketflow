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
        composition: Composition = Composition(),
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

        # Energy is solved only when inlet enthalpy and internal enthalpy exist.
        self._solve_energy = (
            self.total_enthalpy_in.is_assigned
            and self.enthalpy.is_assigned
        )

        # Default outlet enthalpies to the internal enthalpy.
        if self.enthalpy.is_assigned and not self.total_enthalpy_out1.is_assigned:
            self.total_enthalpy_out1 = self.enthalpy

        if self.enthalpy.is_assigned and not self.total_enthalpy_out2.is_assigned:
            self.total_enthalpy_out2 = self.enthalpy

        # Add inlet species to the internal composition without changing values.
        if self.composition_in.is_assigned:
            self.composition.copy_from(
                self.composition_in,
                copy_values=False,
            )

        # Species splitting is active when internal composition and outlet 1 are known.
        self._solve_species = (
            self.composition.is_assigned
            and self.composition_out1.is_assigned
        )

        if self._solve_species:
            # Outlet 1 cannot contain species absent from the internal mixture.
            extra_species = (
                set(self.composition_out1.species)
                - set(self.composition.species)
            )

            if extra_species:
                raise ValueError(
                    f"{self.name}: composition_out1 contains species not in "
                    f"composition: {extra_species}"
                )

            # Give outlet 2 the same species basis as the internal composition.
            self.composition_out2.copy_from(
                self.composition,
                copy_values=False,
            )

        elif self.composition.is_assigned:
            # Normal splitter: both outlets inherit the internal composition.
            self.composition_out1.copy_from(self.composition)
            self.composition_out2.copy_from(self.composition)

    def evaluate_states(self) -> None:
        # Steady-state behavior: inlet composition overrides internal composition.
        if self.composition_in.is_assigned:
            self.composition.copy_from(
                self.composition_in,
                copy_values=True,
            )

        if not self._solve_species:
            return

        mdot1 = self.mass_flow_out1.value
        mdot2 = self.mass_flow_out2.value
        mdot_total = mdot1 + mdot2

        if abs(mdot2) < 1e-12:
            raise ValueError(
                f"{self.name}: cannot compute composition_out2 because "
                f"mass_flow_out2 is zero."
            )

        # Species conservation across the outlet split.
        for species, _ in self.composition:
            self.composition_out2[species].value = (
                mdot_total * self.composition[species].value
                - mdot1 * self.composition_out1[species].value
            ) / mdot2

        self.composition_out2.validate()

    @property
    def iteration_variables(self) -> list[State]:
        variables = [self.pressure]

        # Enthalpy is only an iteration variable when energy is solved.
        if self._solve_energy:
            variables.append(self.enthalpy)

        return variables

    @property
    def residuals(self) -> list[float]:
        residuals = [
            # Total mass conservation.
            self.mass_flow_in.value
            - self.mass_flow_out1.value
            - self.mass_flow_out2.value
        ]

        if self._solve_energy:
            qdot = self.heat_rate.value if self.heat_rate.is_assigned else 0.0

            residuals.append(
                # Steady-flow energy conservation.
                self.mass_flow_in.value * self.total_enthalpy_in.value
                - self.mass_flow_out1.value * self.total_enthalpy_out1.value
                - self.mass_flow_out2.value * self.total_enthalpy_out2.value
                + qdot
            )

        return residuals