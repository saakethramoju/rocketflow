from __future__ import annotations

import numpy as np
from typing import TYPE_CHECKING

from System import Component, State

if TYPE_CHECKING:
    from System import Network, Composition


class FlowSplitter(Component):
    '''
    If you don't assign a composition to an outlet,
    it's assumed to be the same as the inlet.
    '''

    def __init__(
        self,
        name: str,
        network: Network,
        mass_flow_in: State,
        mass_flow_out1: State,
        composition_in: Composition,
        composition_out1: Composition,
        mass_flow_out2: State | float = 10,
        composition_out2: Composition = 1,
    ):
        self.setup()

        if not self.composition_out1.is_assigned:
            self.composition_out1 = self.composition_in

        if not self.composition_out1 <= self.composition_in:
            raise ValueError(
                f"{self.name}: composition_out1 contains species not present "
                f"in composition_in."
            )

    def evaluate_states(self): pass