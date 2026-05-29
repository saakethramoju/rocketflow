from __future__ import annotations

import numpy as np
from typing import TYPE_CHECKING

from System import Component, State

if TYPE_CHECKING:
    from System import Network


class FlowSplitter(Component): pass

'''
If you don't assign a composition to an outlet,
it's assumed to be the same as the inlet. If you
do, that will be set to composition_out, and the
node will take that as composition_in (so add 
composition_in to the nodes). 
'''