from __future__ import annotations

from typing import TYPE_CHECKING
from abc import ABC, abstractmethod


from System import Variable

if TYPE_CHECKING:
    from System import Network

class Component(ABC):

    def __init__(self, 
                 name: str,
                 network: Network):

        self.initialize_component(name, network)


    def initialize_component(self, name: str, network: Network) -> None:
        self.name = name
        self.network = network
        self.network.add_component(component=self)


    #@abstractmethod
    def pre_evaluation(self) -> None:
        pass

    # never put a derived State in the solver’s iteration-variable list!!!
    @property
    @abstractmethod
    def iteration_variables(self) -> list[Variable]:
        return []

    @abstractmethod
    def evaluate_states(self) -> None:
        pass

    @property
    @abstractmethod
    def residuals(self) -> list[float]:
        return []

    def __str__(self):
        return f"Component ({self.__class__.__name__}: {self.name})"