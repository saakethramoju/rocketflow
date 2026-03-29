from typing import Optional
from Global import SIUnits


class Component:
    """
    Base class for every component
    """

    def __init__(self, name: Optional[str] = None):
        self.name = name

        self.system_variables = {}
        self.configuration_variables = {}
        self.iteration_variables = {}
