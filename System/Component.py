from __future__ import annotations

from typing import TYPE_CHECKING
from abc import ABC, abstractmethod

from System.State import State

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
    def iteration_variables(self) -> list[State]:
        return []

    @abstractmethod
    def evaluate_states(self) -> None:
        pass

    @property
    @abstractmethod
    def residuals(self) -> list[float]:
        return []
    
    @property
    def timestep_variables(self) -> list[State]:
        return []
    
    @property
    def time_derivative(self) -> list[State]:
        return []

    def __repr__(self):
        return f"Component ({self.__class__.__name__}: {self.name})"
        
    def __str__(self):

        lines = [f"Component {self.name} ({self.__class__.__name__})"]

        skip_attrs = {"network"}

        for attr, value in self.__dict__.items():

            if attr in skip_attrs:
                continue

            if isinstance(value, State):
                formatted_value = value.value

            elif isinstance(value, list):
                formatted_value = [
                    item.value if isinstance(item, State) else item
                    for item in value
                ]

            elif isinstance(value, tuple):
                formatted_value = tuple(
                    item.value if isinstance(item, State) else item
                    for item in value
                )

            elif isinstance(value, dict):
                formatted_value = {
                    k: (v.value if isinstance(v, State) else v)
                    for k, v in value.items()
                }

            else:
                formatted_value = value

            lines.append(f"    {attr}: {formatted_value}")

        return "\n".join(lines)