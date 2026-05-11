from __future__ import annotations

import numpy as np
from typing import TYPE_CHECKING

from System import Component, State

if TYPE_CHECKING:
    from System import Network


class DarcyWeisbach(Component):
    
    def __init__(self, 
                 name: str, 
                 network: Network,
                 upstream_pressure: State,
                 downstream_pressure: State,
                 darcy_friction_factor: State,
                 ):
        pass