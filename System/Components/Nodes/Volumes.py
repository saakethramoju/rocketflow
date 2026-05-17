from __future__ import annotations

from typing import TYPE_CHECKING

from System import Component, State

if TYPE_CHECKING:
    from System import Network



class SimpleIncompressibleVolume(Component):

    def __init__(self,
                 name: str,
                 network: Network,
                 pressure: State ,
                 volume: float,
                 density: State | None = None,
                 mass_flow_in: State | None = None,
                 mass_flow_out: State | None = None):

        self.setup()

    @property
    def iteration_variables(self) -> list[State]:
        return [self.pressure]

    @property
    def residuals(self) -> list[float]:
        return [self.mdot_in.value - self.mdot_out.value]


class IsothermalIncompressibleVolume(Component):

    def __init__(self,
                 name:str,
                 network: Network,
                 pressure: State,
                 temperature: State,
                 volume: float,
                 density: State | None = None,
                 mass_flow_in: State | None = None,
                 mass_flow_out: State | None = None):
        
        self.setup()
    
    @property
    def iteration_variables(self) -> list[State]:
        return [self.pressure]

    @property
    def residuals(self) -> list[float]:
        return [self.mass_flow_in.value - self.mass_flow_out.value]


class SimpleVolume(Component):

    def __init__(self,
                 name:str,
                 network: Network,
                 pressure: State,
                 enthalpy: State,
                 volume: float,
                 temperature: State | None = None,
                 density: State | None = None,
                 internal_energy: State | None = None,
                 mass_flow_in: State | None = None,
                 mass_flow_out: State | None = None,
                 enthalpy_in: State | None = None,):
        
        self.setup()
    
    @property
    def iteration_variables(self) -> list[State]:
        return [self.pressure, self.enthalpy]

    @property
    def residuals(self) -> list[float]:
        return [self.mass_flow_in.value - self.mass_flow_out.value,
                (self.mass_flow_in.value * self.enthalpy_in.value) - (self.mass_flow_out.value * self.enthalpy.value)]
    '''
    @property
    def residual_scalar(self) -> list[float]:

        mdot_scale = max(
            abs(self.mass_flow_in.value),
            abs(self.mass_flow_out.value),
            1e-9,
        )

        energy_scale = max(
            abs(self.mass_flow_in.value * self.enthalpy_in.value),
            abs(self.mass_flow_out.value * self.enthalpy_out.value),
            1e-9,
        )

        return [
            mdot_scale,
            energy_scale,
        ]
    '''


class SimpleFlowSplitter(Component):

    def __init__(self,
                 name:str,
                 network: Network,
                 pressure: State,
                 enthalpy: State,
                 volume: float,
                 temperature: State | None = None,
                 density: State | None = None,
                 internal_energy: State | None = None,
                 mass_flow_in: State | None = None,
                 mass_flow_out1: State | None = None,
                 mass_flow_out2: State | None = None,
                 enthalpy_in: State | None = None,):
        
        self.setup()
    
    @property
    def iteration_variables(self) -> list[State]:
        return [self.pressure, self.enthalpy]

    @property
    def residuals(self) -> list[float]:
        energy_out = (self.mass_flow_out1.value * self.enthalpy.value) + (self.mass_flow_out2.value * self.enthalpy.value)
        return [self.mass_flow_in.value - (self.mass_flow_out1.value + self.mass_flow_out2.value),
                (self.mass_flow_in.value * self.enthalpy_in.value) - energy_out]
    