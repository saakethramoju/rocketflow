from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from System import Variable, Network

class Balance:

    def __init__(self, name: str,
                 network: Network,
                 iteration_variable: Variable,
                 ):
        pass