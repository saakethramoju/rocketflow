from __future__ import annotations

from typing import TYPE_CHECKING

from System import Component, State
from Utilities import Fluid

if TYPE_CHECKING:
    from System import Network


class Helium(Component):pass


class IdealGas(Component):

    def __init__(self, 
                 name: str, 
                 network: Network,
                 molar_mass: float,
                 degrees_of_freedom: State | float = 3,
                 pressure: State | None = None,
                 temperature: State | None = None,
                 universal_gas_constant: float = 8.314,
                 specific_gas_constant: float | None = None,
                 density: float | None = None,
                 constant_pressure_specific_heat: float | None = None,
                 constant_volume_specific_heat: float | None = None,
                 specific_heat_ratio: float | None = None):
        self.setup()