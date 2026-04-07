from abc import ABC, abstractmethod
from Components import Variable

class Component(ABC):

    def __init__(self, 
                 name: str):
        self.name = name

    @property
    @abstractmethod
    def iteration_variables(self) -> list[Variable]:
        pass

    #@abstractmethod
    def pre_state_evaluation(self) -> None:
        pass

    @abstractmethod
    def evaluate_states(self) -> None:
        pass

    @property
    @abstractmethod
    def residuals(self) -> list[float]:
        pass

    def __str__(self):
        return f"Component ({self.__class__.__name__}: {self.name})"