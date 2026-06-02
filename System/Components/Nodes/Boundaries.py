from __future__ import annotations

from typing import TYPE_CHECKING

from System import Component, State

if TYPE_CHECKING:
    from System import Network


class Boundary(Component):
    
    def __init__(self, 
                 name: str,
                 network: Network,
                 pressure: State,
                 temperature: State):
        self.setup()

class PressureBoundary(Component):

    def __init__(self, 
                 name: str,
                 network: Network,
                 pressure: State):
        self.setup()