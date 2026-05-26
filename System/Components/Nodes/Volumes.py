from __future__ import annotations

from typing import TYPE_CHECKING

from System import Component, State

if TYPE_CHECKING:
    from System import Network



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





class SimpleFlowSplitter(Component):
    """
    Assume the flow entering and leaving carries total enthalpy.
    Optional heat_rate is added directly to the splitter control volume.
    Positive heat_rate adds energy to the fluid.
    """
    def __init__(
        self,
        name: str,
        network: Network,
        pressure: State,
        enthalpy: State,
        volume: float,
        total_enthalpy_in: State,
        total_enthalpy_out1: State | None = None,
        total_enthalpy_out2: State | None = None,
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

    @property
    def iteration_variables(self) -> list[State]:
        return [self.pressure, self.enthalpy]

    @property
    def residuals(self) -> list[float]:
        qdot = self.heat_rate.value if self.heat_rate.is_assigned else 0.0

        energy_in = self.mass_flow_in.value * self.total_enthalpy_in.value

        energy_out = (
            self.mass_flow_out1.value * self.total_enthalpy_out1.value
            + self.mass_flow_out2.value * self.total_enthalpy_out2.value
        )

        return [
            self.mass_flow_in.value
            - (self.mass_flow_out1.value + self.mass_flow_out2.value),

            energy_in + qdot - energy_out,
        ]