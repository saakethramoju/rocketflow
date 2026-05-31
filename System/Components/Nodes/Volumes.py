from __future__ import annotations

from typing import TYPE_CHECKING

from System import Component, State, Composition

if TYPE_CHECKING:
    from System import Network



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
        composition: Composition = Composition(),
        composition_in: Composition = Composition(),
        mass_flow_in: State | None = None,
        mass_flow_out: State | None = None,
    ):
        self.setup()
        self.composition.update()

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
        composition: Composition = Composition(),
        composition_in: Composition = Composition(),
        mass_flow_in: State | None = None,
        mass_flow_out: State | None = None
    ):
        self.setup()
        self.composition.update()
    
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
        volume: float,
        enthalpy: State | None = None,
        temperature: State | None = None,
        density: State | None = None,
        internal_energy: State | None = None,
        heat_rate: State | float | None = None,
        total_enthalpy_in: State | None = None,
        total_enthalpy_out: State | None = None,
        mass_flow_in: State | None = None,
        mass_flow_out: State | None = None,
    ):
        self.setup()
        self.composition.update()

        self._solve_energy = self.enthalpy.is_assigned and self.total_enthalpy_in.is_assigned

        if self._solve_energy and not self.total_enthalpy_out.is_assigned:
            self.total_enthalpy_out = self.enthalpy

    @property
    def iteration_variables(self) -> list[State]:
        variables = [self.pressure]

        if self._solve_energy:
            variables.append(self.enthalpy)

        return variables

    @property
    def residuals(self) -> list[float]:
        residuals = [
            self.mass_flow_in.value - self.mass_flow_out.value
        ]

        if self._solve_energy:
            qdot = self.heat_rate.value if self.heat_rate.is_assigned else 0.0

            residuals.append(
                self.mass_flow_in.value * self.total_enthalpy_in.value
                - self.mass_flow_out.value * self.total_enthalpy_out.value
                + qdot
            )

        return residuals
    

'''

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
        composition: Composition = Composition(),
        composition_in: Composition = Composition(),
        mass_flow_in: State | None = None,
        mass_flow_out: State | None = None,
    ):
        self.setup()
        self.composition.update()

    @property
    def iteration_variables(self) -> list[State]:
        return [
            self.pressure,
            *(self.composition[species] for species in self.composition & self.composition_in)
        ]

    @property
    def residuals(self) -> list[float]:
        self.composition.update()

        return [
            self.mass_flow_in.value - self.mass_flow_out.value,
            *(self.mass_flow_in.value * self.composition_in[species].value
                - self.mass_flow_out.value * self.composition[species].value
                for species in self.composition & self.composition_in
            ),
        ]
'''