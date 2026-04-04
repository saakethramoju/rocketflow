from __future__ import annotations

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from .Component import Branch, Node

class Port:
    def __init__(self, 
                 branch: Branch,
                 node: Node):
        self.branch = branch
        self.node = node