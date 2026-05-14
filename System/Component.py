from __future__ import annotations

from typing import TYPE_CHECKING
from abc import ABC, abstractmethod
import inspect

from System.State import State

if TYPE_CHECKING:
    from System import Network

class Component(ABC):

    def __init__(self, 
                 name: str,
                 network: Network):
        
        self.setup()

    def setup(self) -> None:
        """
        Automatically initialize component and constructor attributes.

        Intended to be called inside subclass __init__():

            self.setup()

        This will:
        - initialize the component using `name` and `network`
        - automatically create instance attributes from constructor locals
        - pass values through initialize_attribute()
        """

        # Get caller frame (__init__ frame).
        frame = inspect.currentframe().f_back

        if frame is None:
            raise RuntimeError("Could not access caller frame.")

        local_vars = frame.f_locals

        # Required constructor arguments.
        name = local_vars["name"]
        network = local_vars["network"]

        # Initialize component.
        self.initialize_component(name, network)

        # Automatically assign remaining constructor args.
        for attr_name, value in local_vars.items():

            # Skip internal/setup variables.
            if attr_name in {"self", "name", "network"}:
                continue

            setattr(
                self,
                attr_name,
                self.initialize_attribute(value),
            )


    def initialize_attribute(self, value: State | float | int | str | None = None) -> State:
        """
        Normalize input into a State object.

        - State -> returned directly
        - float/int -> wrapped in State
        - None -> empty placeholder State
        - Anything else -> returned directly
        """

        if isinstance(value, State):
            return value

        if value is None:
            return State()

        if isinstance(value, (float, int)):
            return State(float(value))
        
        return value

        '''
        raise TypeError(
            f"Expected State, float, int, or None; "
            f"got {type(value).__name__}."
        )'''


    def initialize_component(self, name: str, network: Network) -> None:
        self.name = name
        self.network = network
        self.network.add_component(component=self)

    #@abstractmethod
    def pre_evaluation(self) -> None:
        pass

    # never put a derived State in the solver’s iteration-variable list!!!
    @property
    #@abstractmethod
    def iteration_variables(self) -> list[State]:
        return []

    #@abstractmethod
    def evaluate_states(self) -> None:
        pass

    @property
    #@abstractmethod
    def residuals(self) -> list[float]:
        return []
    
    @property
    def residual_scalar(self) -> list[float]:
        return [1.0] * len(self.residuals)

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
        
    @property
    def ignored_export_attributes(self) -> set[str]:
        return set()
    
    # ---- transient stuff ----#
    @property
    def timestep_variables(self) -> list[State]:
        return []
    
    @property
    def time_derivative(self) -> list[State]:
        return []