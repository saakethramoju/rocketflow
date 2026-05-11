from __future__ import annotations

from typing import TYPE_CHECKING

from System import Component, State
from Utilities import Fluid

if TYPE_CHECKING:
    from System import Network


class PressureBoundary(Component):

    def __init__(self, 
                 name: str,
                 network: Network,
                 pressure: State):
        self.setup()
    

class IsothermalPressureBoundary(Component):

    def __init__(self, 
                 name: str,
                 network: Network,
                 pressure: State,
                 density: State,
                 temperature: State):
        
        self.setup()
