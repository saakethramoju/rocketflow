from __future__ import annotations

from typing import TYPE_CHECKING

from System import Component, State, Composition

if TYPE_CHECKING:
    from System import Network



class FlowSplitter(Component):
    """
    Notes
    -----
    1) Species conservation solves for the composition at 
       outlet 2 only.
    2) If species conser
    """    
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

        # Default outlet enthalpies to internal enthalpy.
        if self.enthalpy.is_assigned and not self.total_enthalpy_out1.is_assigned:
            self.total_enthalpy_out1 = self.enthalpy

        if self.enthalpy.is_assigned and not self.total_enthalpy_out2.is_assigned:
            self.total_enthalpy_out2 = self.enthalpy

        # Steady-state: inlet composition overrides internal composition.
        if self.composition_in.is_assigned:
            self.composition.copy_from(
                self.composition_in,
                copy_values=True,
            )

        # Species splitting is only possible if internal composition and outlet1 exist.
        self._solve_species = (
            self.composition.is_assigned
            and self.composition_out1.is_assigned
        )

        if self._solve_species:
            # Outlet 1 cannot contain species absent from internal composition.
            extra_species = (
                set(self.composition_out1.species)
                - set(self.composition.species)
            )

            if extra_species:
                raise ValueError(
                    f"{self.name}: composition_out1 contains species not in "
                    f"composition: {extra_species}"
                )

            # Outlet 2 gets the same species basis, but values are solved later.
            self.composition_out2.copy_from(
                self.composition,
                copy_values=True,
            )

        elif self.composition.is_assigned:
            # Ordinary splitter: both outlets match internal composition.
            self.composition_out1.copy_from(
                self.composition,
                copy_values=True,
            )

            self.composition_out2.copy_from(
                self.composition,
                copy_values=True,
            )


    def evaluate_states(self) -> None:
        # Steady-state: inlet composition overrides internal composition.
        if self.composition_in.is_assigned:
            self.composition.copy_from(
                self.composition_in,
                copy_values=True,
            )

        if not self._solve_species:
            return

        # Outlet mdots may be placeholder States before the solver assigns them.
        if (
            not self.mass_flow_out1.is_assigned
            or not self.mass_flow_out2.is_assigned
        ):
            return

        mdot1 = self.mass_flow_out1.value
        mdot2 = self.mass_flow_out2.value
        mdot_total = mdot1 + mdot2

        if abs(mdot2) < 1e-12:
            raise ValueError(
                f"{self.name}: cannot compute composition_out2 because "
                f"mass_flow_out2 is zero."
            )

        # Solve outlet 2 composition from species conservation.
        for species, _ in self.composition:
            self.composition_out2[species].value = (
                mdot_total * self.composition[species].value
                - mdot1 * self.composition_out1[species].value
            ) / mdot2

        self.composition_out2.validate()



    @property
    def iteration_variables(self) -> list[State]:
        variables = [self.pressure]

        if self._solve_energy:
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